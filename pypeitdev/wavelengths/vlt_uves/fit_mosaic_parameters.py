import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pypeit import wavecalib
from pypeit import edgetrace

from scipy.optimize import minimize
from scipy.optimize import brute
from numpy.polynomial.polynomial import polyval2d


# ---------------------------------------------------------------------------
# 2D polynomial wavelength model (no PypeIt dependency)
# ---------------------------------------------------------------------------

def fit_wave2d(spec, spat, wave, order_spec=6, order_spat=4):
    """
    Fit a 2D polynomial to the wavelength solution across all orders.

    Inputs are normalised internally to [-1, 1] for numerical stability.

    Parameters
    ----------
    spec : np.ndarray, shape (N,)
        Spectral pixel coordinates (dispersion direction).
    spat : np.ndarray, shape (N,)
        Spatial pixel coordinates (cross-dispersion direction).
    wave : np.ndarray, shape (N,)
        Wavelength values corresponding to each (spec, spat) point.
    order_spec : int
        Polynomial degree along the spectral direction.
    order_spat : int
        Polynomial degree along the spatial direction.

    Returns
    -------
    coeffs : np.ndarray, shape (order_spec+1, order_spat+1)
        2D polynomial coefficients.
    norm_params : dict
        Normalisation parameters needed to evaluate the fit on new data.
    residuals : np.ndarray, shape (N,)
        Fit residuals (wave - wave_fit) for each input point.
    rms : float
        RMS of residuals in wavelength units.
    """
    spec_min, spec_max = spec.min(), spec.max()
    spat_min, spat_max = spat.min(), spat.max()

    spec_norm = 2.0 * (spec - spec_min) / (spec_max - spec_min) - 1.0
    spat_norm = 2.0 * (spat - spat_min) / (spat_max - spat_min) - 1.0

    norm_params = dict(spec_min=spec_min, spec_max=spec_max,
                       spat_min=spat_min, spat_max=spat_max)

    n_terms = (order_spec + 1) * (order_spat + 1)
    A = np.zeros((len(spec), n_terms))
    col = 0
    for i in range(order_spec + 1):
        for j in range(order_spat + 1):
            A[:, col] = spec_norm**i * spat_norm**j
            col += 1

    result = np.linalg.lstsq(A, wave, rcond=None)
    coeffs_flat = result[0]
    coeffs = coeffs_flat.reshape(order_spec + 1, order_spat + 1)

    wave_fit = A @ coeffs_flat
    residuals = wave - wave_fit
    rms = np.sqrt(np.mean(residuals**2))

    return coeffs, norm_params, residuals, rms


def eval_wave2d(spec, spat, coeffs, norm_params):
    """
    Evaluate a 2D polynomial wavelength model at arbitrary (spec, spat) positions.

    Parameters
    ----------
    spec, spat : np.ndarray
        Coordinates in original pixel units.
    coeffs : np.ndarray, shape (order_spec+1, order_spat+1)
        Coefficients from fit_wave2d.
    norm_params : dict
        Normalisation parameters from fit_wave2d.

    Returns
    -------
    wave : np.ndarray
        Predicted wavelengths.
    """
    spec_norm = 2.0 * (spec - norm_params['spec_min']) / (norm_params['spec_max'] - norm_params['spec_min']) - 1.0
    spat_norm = 2.0 * (spat - norm_params['spat_min']) / (norm_params['spat_max'] - norm_params['spat_min']) - 1.0
    return polyval2d(spec_norm, spat_norm, coeffs)


# ---------------------------------------------------------------------------
# Rigid-body transform: detector 2 -> detector 1 frame
# ---------------------------------------------------------------------------

def transform_det2(spec2, spat2, dx, dy, theta):
    """
    Apply a rigid-body transform to detector 2 coordinates.

    Rotates about the centroid of detector 2's points, then translates.

    Parameters
    ----------
    spec2, spat2 : np.ndarray
        Spectral and spatial coordinates in detector 2's local frame.
    dx : float
        Shift in the spatial (cross-dispersion) direction, in pixels.
    dy : float
        Shift in the spectral (dispersion) direction, in pixels.
    theta : float
        Rotation angle in radians (positive = counter-clockwise).

    Returns
    -------
    spec2_t, spat2_t : np.ndarray
        Transformed coordinates in detector 1's frame.
    """
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    spat2_c = spat2 - np.mean(spat2)
    spec2_c = spec2 - np.mean(spec2)

    spat2_t = cos_t * spat2_c - sin_t * spec2_c + np.mean(spat2) + dx
    spec2_t = sin_t * spat2_c + cos_t * spec2_c + np.mean(spec2) + dy

    return spec2_t, spat2_t


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def unpack_data(alldata, allorders):
    """
    Unpack alldata and allorders into flat arrays for each detector,
    applying the per-pixel validity mask stored in channel 3.

    Parameters
    ----------
    alldata : list of list of np.ndarray, shape (nsetups, 2, nspec, norders, 4)
        alldata[setup][detector] has shape (nspec, norders, 4):
          channel 0: wavelength
          channel 1: spec pixel
          channel 2: spat pixel (offset so det2 follows det1)
          channel 3: validity mask (bool/float, nonzero = good)
    allorders : list of list of np.ndarray
        allorders[setup][detector] contains the echelle order numbers.

    Returns
    -------
    spec1, spat1, wave1 : np.ndarray
        Masked, flattened arrays for detector 1.
    spec2, spat2, wave2 : np.ndarray
        Masked, flattened arrays for detector 2.
    """
    spec1_list, spat1_list, wave1_list = [], [], []
    spec2_list, spat2_list, wave2_list = [], [], []

    for ii in range(len(alldata)):
        for det, (spec_list, spat_list, wave_list) in enumerate(
                [(spec1_list, spat1_list, wave1_list),
                 (spec2_list, spat2_list, wave2_list)]):
            data = alldata[ii][det]          # shape (nspec, norders, 4)
            nspec, norders, _ = data.shape
            for ord_idx in range(norders):
                ww = np.where(data[:, ord_idx, 3])[0]
                wave_list.append(data[ww, ord_idx, 0])
                spec_list.append(data[ww, ord_idx, 1])
                spat_list.append(data[ww, ord_idx, 2])

    return (np.concatenate(spec1_list), np.concatenate(spat1_list), np.concatenate(wave1_list),
            np.concatenate(spec2_list), np.concatenate(spat2_list), np.concatenate(wave2_list))


# ---------------------------------------------------------------------------
# Cost function (pooled fit across both detectors)
# ---------------------------------------------------------------------------

def cost_function(params, spec1, spat1, wave1, spec2, spat2, wave2,
                  order_spec=6, order_spat=4):
    """
    Cost function for the rigid-body alignment optimisation.

    Applies the transform to detector 2, pools all points with detector 1,
    fits a single 2D polynomial wavelength model to the combined data,
    and returns the RMS residual.

    Parameters
    ----------
    params : array-like, length 3
        (dx, dy, theta) — spatial shift, spectral shift, rotation in radians.
    spec1, spat1, wave1 : np.ndarray
        Detector 1 data (all orders and setups, masked and flattened).
    spec2, spat2, wave2 : np.ndarray
        Detector 2 data (all orders and setups, masked and flattened).
    order_spec, order_spat : int
        Polynomial degrees for the 2D wavelength fit.

    Returns
    -------
    rms : float
        RMS wavelength residual of the global 2D fit, in wavelength units.
    """
    dx, dy, theta = params

    # Transform detector 2 into detector 1's frame
    spec2_t, spat2_t = transform_det2(spec2, spat2, dx, dy, theta)

    # Pool all data and fit a single 2D polynomial
    spec_all = np.concatenate([spec1, spec2_t])
    spat_all = np.concatenate([spat1, spat2_t])
    wave_all = np.concatenate([wave1, wave2])

    _, _, _, rms = fit_wave2d(spec_all, spat_all, wave_all,
                              order_spec=order_spec, order_spat=order_spat)
    return rms


# ---------------------------------------------------------------------------
# Bootstrap uncertainty estimation
# ---------------------------------------------------------------------------

def bootstrap_errors(best_params, spec1, spat1, wave1, spec2, spat2, wave2,
                     order_spec=6, order_spat=4,
                     n_bootstrap=200, pixel_sigma=0.1,
                     bounds=None, verbose=True):
    """
    Estimate parameter uncertainties by adding Gaussian pixel-coordinate
    noise and re-optimising from the best-fit solution.

    For each bootstrap iteration:
      1. Add Gaussian noise (sigma = pixel_sigma) to all pixel coordinates.
      2. Re-run the fine L-BFGS-B optimiser starting from the known solution.
      3. Collect the resulting (dx, dy, theta) triple.

    The standard deviation across all iterations gives 1-sigma uncertainties.

    Parameters
    ----------
    best_params : array-like, length 3
        Best-fit (dx, dy, theta) from the main solver.
    spec1, spat1, wave1 : np.ndarray
        Detector 1 data (masked, flattened).
    spec2, spat2, wave2 : np.ndarray
        Detector 2 data (masked, flattened).
    order_spec, order_spat : int
        Polynomial degrees (must match those used in the main fit).
    n_bootstrap : int
        Number of bootstrap iterations.
    pixel_sigma : float
        Standard deviation of the Gaussian pixel-coordinate noise (pixels).
    bounds : list of (min, max) or None
        Parameter bounds passed to the optimiser.
    verbose : bool
        Print progress.

    Returns
    -------
    errors : dict with keys 'dx', 'dy', 'theta'
        1-sigma uncertainties on each parameter.
    samples : np.ndarray, shape (n_bootstrap, 3)
        Full bootstrap sample array for (dx, dy, theta).
    """
    rng = np.random.default_rng(seed=42)
    samples = np.zeros((n_bootstrap, 3))

    if verbose:
        print(f"\nRunning {n_bootstrap} bootstrap iterations (pixel_sigma={pixel_sigma} px)...")

    for i in range(n_bootstrap):
        if verbose and (i + 1) % 50 == 0:
            print(f"  Bootstrap iteration {i + 1}/{n_bootstrap}")

        # Perturb pixel coordinates with Gaussian noise
        s1 = spec1 + rng.normal(0, pixel_sigma, size=spec1.shape)
        p1 = spat1 + rng.normal(0, pixel_sigma, size=spat1.shape)
        s2 = spec2 + rng.normal(0, pixel_sigma, size=spec2.shape)
        p2 = spat2 + rng.normal(0, pixel_sigma, size=spat2.shape)

        args = (s1, p1, wave1, s2, p2, wave2, order_spec, order_spat)

        res = minimize(cost_function, x0=best_params, args=args,
                       method='L-BFGS-B', bounds=bounds,
                       options={'ftol': 1e-10, 'gtol': 1e-6, 'maxiter': 500})
        samples[i] = res.x

    errors = {
        'dx':    samples[:, 0].std(),
        'dy':    samples[:, 1].std(),
        'theta': samples[:, 2].std(),
    }

    if verbose:
        print(f"Bootstrap errors:  dx={errors['dx']:.4f} pix,  "
              f"dy={errors['dy']:.4f} pix,  "
              f"theta={np.degrees(errors['theta']):.6f} deg")

    return errors, samples


# ---------------------------------------------------------------------------
# Diagnostic plots
# ---------------------------------------------------------------------------

def plot_mosaic(alldata, allorders, result=None,
                det_nspec=4096, det_nspat=2042):
    """
    Plot the echelle order traces for both detectors with detector outlines.

    Two panels are produced:
      Left  — raw data as loaded (detector 2 offset in spat, spec local).
      Right — after applying the best-fit transform to detector 2
               (only drawn if `result` is provided).

    Parameters
    ----------
    alldata : list of list of np.ndarray
        As returned by load_data(), shape (nsetups, 2, nspec, norders, 4).
    allorders : list of list of np.ndarray
        As returned by load_data().
    result : dict or None
        Output of fit_mosaic_parameters(). If provided, the right panel
        shows the aligned coordinates.
    det_nspec : int
        Number of pixels in the spectral direction (default 4096).
    det_nspat : int
        Number of pixels in the spatial direction per detector (default 2048).
    """
    spec1, spat1, wave1, spec2, spat2, wave2 = unpack_data(alldata, allorders)

    n_panels = 2 if result is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6),
                             sharex=False, sharey=False)
    if n_panels == 1:
        axes = [axes]

    def draw_panel(ax, s1, p1, s2, p2, title, offx=0, offy=0):
        det1_rect = mpatches.Rectangle((0, 0), det_nspat, det_nspec,
                                       linewidth=1.5, edgecolor='steelblue',
                                       facecolor='none', label='Detector 1')
        det2_rect = mpatches.Rectangle((det_nspat + offx, offy), det_nspat, det_nspec,
                                       linewidth=1.5, edgecolor='darkorange',
                                       facecolor='none', label='Detector 2')
        ax.add_patch(det1_rect)
        ax.add_patch(det2_rect)

        sc = ax.scatter(p1, s1, c=wave1, cmap='viridis', s=0.5, alpha=0.4, rasterized=True)
        ax.scatter(p2, s2, c=wave2, cmap='plasma',  s=0.5, alpha=0.4, rasterized=True)

        ax.set_xlabel('Spatial pixel (spat)')
        ax.set_ylabel('Spectral pixel (spec)')
        ax.set_title(title)
        ax.set_aspect('equal')
        ax.legend(handles=[det1_rect, det2_rect], loc='upper right', fontsize=8)
        plt.colorbar(sc, ax=ax, label='Wavelength det 1 (Å/order)')

        margin = 100
        ax.set_xlim(-margin, 2 * det_nspat + margin)
        ax.set_ylim(-margin, det_nspec + margin)

    draw_panel(axes[0], spec1, spat1, spec2, spat2, 'Raw (pre-alignment)')

    if result is not None:
        spec2_t, spat2_t = transform_det2(spec2, spat2,
                                          result['dx'], result['dy'], result['theta'])
        draw_panel(axes[1], spec1, spat1, spec2_t, spat2_t, 'After alignment', offx=result['dx'], offy=result['dy'])

    plt.tight_layout()
    plt.savefig('mosaic_alignment.pdf', dpi=150, bbox_inches='tight')
    plt.show()
    print("Plot saved to mosaic_alignment.pdf")


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

def fit_mosaic_parameters(alldata, allorders,
                           dx_range=(-200, 200), dy_range=(-50, 50),
                           theta_range=(-0.01, 0.01),
                           coarse_steps=(21, 15, 11),
                           order_spec=6, order_spat=4,
                           n_bootstrap=200, pixel_sigma=0.1,
                           verbose=True):
    """
    Find the optimal rigid-body alignment (dx, dy, theta) of detector 2
    relative to detector 1 by minimising the RMS of a global 2D polynomial
    wavelength fit to both detectors pooled together.

    Uses a coarse brute-force grid search followed by a fine L-BFGS-B
    optimisation. Parameter uncertainties are estimated by bootstrap
    resampling with Gaussian pixel-coordinate perturbations.

    Parameters
    ----------
    alldata, allorders : as returned by load_data().
    dx_range : (min, max) search range for the spatial shift in pixels.
    dy_range : (min, max) search range for the spectral shift in pixels.
    theta_range : (min, max) search range for the rotation in radians.
    coarse_steps : (ndx, ndy, ntheta) grid points for brute-force search.
    order_spec, order_spat : polynomial degrees for the 2D wavelength fit.
    n_bootstrap : number of bootstrap iterations for uncertainty estimation.
    pixel_sigma : pixel-coordinate noise std dev for bootstrap (pixels).
    verbose : print progress and results.

    Returns
    -------
    result : dict with keys
        'dx', 'dy', 'theta'      — best-fit parameters
        'dx_err', 'dy_err',
        'theta_err'              — 1-sigma bootstrap uncertainties
        'rms'                    — RMS wavelength residual at optimum
        'bootstrap_samples'      — shape (n_bootstrap, 3) array
        'coarse_result'          — scipy.optimize.brute output
        'fine_result'            — scipy.optimize.minimize output
    """
    spec1, spat1, wave1, spec2, spat2, wave2 = unpack_data(alldata, allorders)

    if verbose:
        print(f"Loaded {len(spec1)} points from detector 1, {len(spec2)} from detector 2.")
        print(f"Det1 spat range: {spat1.min():.1f} – {spat1.max():.1f}")
        print(f"Det2 spat range: {spat2.min():.1f} – {spat2.max():.1f}")
        print(f"Det1 wave range: {wave1.min():.1f} – {wave1.max():.1f}")
        print(f"Det2 wave range: {wave2.min():.1f} – {wave2.max():.1f}")

    args = (spec1, spat1, wave1, spec2, spat2, wave2, order_spec, order_spat)
    bounds = [dx_range, dy_range, theta_range]

    # --- Stage 1: coarse brute-force search ---
    if verbose:
        print("\nRunning coarse grid search...")

    ranges = (slice(dx_range[0],    dx_range[1],    complex(coarse_steps[0])),
              slice(dy_range[0],    dy_range[1],    complex(coarse_steps[1])),
              slice(theta_range[0], theta_range[1], complex(coarse_steps[2])))

    coarse_result = brute(cost_function, ranges, args=args, finish=None)
    dx0, dy0, theta0 = coarse_result

    if verbose:
        rms0 = cost_function(coarse_result, *args)
        print(f"Coarse result: dx={dx0:.2f} pix, dy={dy0:.2f} pix, "
              f"theta={np.degrees(theta0):.4f} deg, RMS={rms0:.6f}")

    # --- Stage 2: fine gradient-based optimisation ---
    if verbose:
        print("Running fine optimisation...")

    fine_result = minimize(cost_function, x0=[dx0, dy0, theta0], args=args,
                           method='L-BFGS-B', bounds=bounds,
                           options={'ftol': 1e-12, 'gtol': 1e-8, 'maxiter': 1000})

    dx_fit, dy_fit, theta_fit = fine_result.x
    rms_fit = fine_result.fun

    if verbose:
        print(f"Fine result:   dx={dx_fit:.4f} pix, dy={dy_fit:.4f} pix, "
              f"theta={np.degrees(theta_fit):.6f} deg, RMS={rms_fit:.6f}")
        print(f"Optimisation {'succeeded' if fine_result.success else 'FAILED: ' + fine_result.message}")

    # --- Stage 3: bootstrap uncertainty estimation ---
    errors, samples = bootstrap_errors(
        best_params=[dx_fit, dy_fit, theta_fit],
        spec1=spec1, spat1=spat1, wave1=wave1,
        spec2=spec2, spat2=spat2, wave2=wave2,
        order_spec=order_spec, order_spat=order_spat,
        n_bootstrap=n_bootstrap, pixel_sigma=pixel_sigma,
        bounds=bounds, verbose=verbose)

    return {
        'dx':                dx_fit,
        'dy':                dy_fit,
        'theta':             theta_fit,
        'dx_err':            errors['dx'],
        'dy_err':            errors['dy'],
        'theta_err':         errors['theta'],
        'rms':               rms_fit,
        'bootstrap_samples': samples,
        'coarse_result':     coarse_result,
        'fine_result':       fine_result,
    }


# ---------------------------------------------------------------------------
# Data loader
# ---------------------------------------------------------------------------

def load_data():
    dirc = 'PypeIt_Templates/'
    # setups = [['564_l', '564_u'],
    #           ['580_l', '580_u']]
    setups = [['580_l', '580_u']]
    alldata   = [[None for _ in range(len(setups[0]))] for _ in setups]
    allorders = [[None for _ in range(len(setups[0]))] for _ in setups]

    for ii in range(len(setups)):
        for dd in range(2):
            grat, chip = setups[ii][dd].split('_')
            edges = edgetrace.EdgeTraceSet.from_file(f'{grat}/Edges_{chip}.fits.gz')
            offset = 0.0
            if dd == 0:
                next_offset = edges.traceimg.image.shape[1]
            else:
                offset = next_offset

            edge_fit = edges.edge_fit
            slitcen = np.mean(
                edge_fit.reshape((edge_fit.shape[0], edge_fit.shape[1] // 2, 2)), axis=2)
            nspec, norders = slitcen.shape

            wvcal_file = f"{dirc}/WaveCalib_{grat}{chip}.fits"
            wvcal = wavecalib.WaveCalib.from_file(wvcal_file)
            orders = wvcal.ech_orders
            assert orders.size == norders

            thisdata = np.zeros((nspec, norders, 4))
            for ord in range(norders):
                wave = (wvcal.wv_fit2d[0].eval(
                    np.linspace(0.0, 1.0, nspec),
                    x2=np.full(nspec, orders[ord])) / orders[ord])
                spec = np.arange(nspec, dtype=float)
                spat = slitcen[:, ord] + offset
                mask = ((slitcen[:, ord] >= 0) &
                        (slitcen[:, ord] < edges.traceimg.image.shape[1]))
                thisdata[:, ord, 0] = wave
                thisdata[:, ord, 1] = spec
                thisdata[:, ord, 2] = spat
                thisdata[:, ord, 3] = mask

            alldata[ii][dd]   = thisdata.copy()
            allorders[ii][dd] = orders.copy()

    return alldata, allorders


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    alldata, allorders = load_data()

    result = fit_mosaic_parameters(
        alldata, allorders,
        dx_range=(-200, 200),
        dy_range=(-50, 50),
        theta_range=(0.0, 0.0),
        coarse_steps=(21, 15, 1),
        n_bootstrap=200,
        pixel_sigma=0.1)

    print("\nFinal result:")
    print(f"  dx    = {result['dx']:.4f} ± {result['dx_err']:.4f} pixels")
    print(f"  dy    = {result['dy']:.4f} ± {result['dy_err']:.4f} pixels")
    print(f"  theta = {np.degrees(result['theta']):.6f} ± "
          f"{np.degrees(result['theta_err']):.6f} degrees")
    print(f"  RMS   = {result['rms']:.6f} (wavelength units)")

    plot_mosaic(alldata, allorders, result=result,
                det_nspec=4096, det_nspat=2042)
