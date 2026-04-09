import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pypeit import wavecalib
from pypeit import edgetrace
from pypeit.core import fitting
from pypeit.core.arc import fit2darc, fit2darc_global_qa, fit2darc_orders_qa

from scipy.optimize import minimize
from scipy.optimize import brute

from IPython import embed

det_nspec=4096
det_nspat=2042

# ---------------------------------------------------------------------------
# Rigid-body transform: detector 2 -> detector 1 frame
# ---------------------------------------------------------------------------

def transform_det2(all_specpix, all_spatpix, all_detector, dx, dy, theta):
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

    out_specpix, out_spatpix = all_specpix.copy(), all_spatpix.copy()
    d2 = np.where(all_detector == 1)

    spat2_c = all_spatpix[d2] - det_nspat/2
    spec2_c = all_specpix[d2] - det_nspec/2

    out_spatpix[d2] = cos_t * spat2_c - sin_t * spec2_c + det_nspat/2 + dx
    out_specpix[d2] = sin_t * spat2_c + cos_t * spec2_c + det_nspec/2 + dy

    return out_specpix, out_spatpix

# ---------------------------------------------------------------------------
# Cost function (pooled fit across both detectors)
# ---------------------------------------------------------------------------

def cost_function(params, all_specpix, all_spatpix, all_wave,
                  all_orders, all_detector, all_setup, all_slitcen, all_ordrcen,
                  nspec_coeff=6, norder_coeff=4):
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
    nspec_coeff, norder_coeff : int
        Polynomial degrees for the 2D wavelength fit.

    Returns
    -------
    rms : float
        RMS wavelength residual of the global 2D fit, in wavelength units.
    """
    dx, dy, theta = params

    # Transform detector 2 into detector 1's frame
    all_specpix_t, all_spatpix_t = transform_det2(all_specpix, all_spatpix, all_detector, dx, dy, theta)

    # Analyse each setup separately
    unq_setup = np.unique(all_setup)
    all_wavefit = np.zeros_like(all_wave)
    all_resid = np.zeros_like(all_wave)
    for ss in range(unq_setup.size):
        wgd = np.where(all_setup == unq_setup[ss])
        fit2d = fit2darc(all_wave[wgd], all_specpix_t[wgd], all_orders[wgd], det_nspec,
                         nspec_coeff=nspec_coeff, norder_coeff=norder_coeff)

        # Use this 2D fit to get the wavelength solution for each order
        waveord_fit = fit2d.eval(all_specpix_t[wgd]/(det_nspec-1), x2=all_orders[wgd])
        waveord_fit_p1 = fit2d.eval((1+all_specpix_t[wgd])/(det_nspec-1), x2=all_orders[wgd])
        dwav_dpix = np.abs(waveord_fit_p1 - waveord_fit) / all_orders[wgd]
        all_wavefit[wgd] = waveord_fit / all_orders[wgd]
        # Save some QA
        fit2darc_global_qa(fit2d, det_nspec, outfile=f"global_qa_{unq_setup[ss]}.pdf")
        fit2darc_orders_qa(fit2d, det_nspec, outfile=f"orders_qa_{unq_setup[ss]}.pdf")

        # Finally, calculate the spectral residual in pixels
        spec_resid = (all_wave[wgd] - all_wavefit[wgd]) / dwav_dpix

        ##########################
        # Now do the spatial fit as well, and combine the residuals in quadrature
        ##########################
        all_slitcen_spat = np.array([])
        all_slitcen_spec = np.array([])
        all_slitcen_ordr = np.array([])
        all_slitcen_det = np.array([])
        for dd in range(len(all_slitcen[ss])):
            all_speccen = np.arange(det_nspec).reshape(-1, 1).repeat(all_slitcen[ss][dd].shape[1], axis=1)
            all_ordrcen_t = np.tile(all_ordrcen[ss][dd].reshape(1, -1), (det_nspec, 1))
            # Mask of good values that are actually on the detector
            if dd == 0:
                mask = ((all_slitcen[ss][dd] >= 0) & (all_slitcen[ss][dd] < det_nspat))
            else:
                mask = ((all_slitcen[ss][dd] >= det_nspat) & (all_slitcen[ss][dd] < 2*det_nspat))
            all_slitcen_spec = np.concatenate([all_slitcen_spec, all_speccen[mask]])
            all_slitcen_spat = np.concatenate([all_slitcen_spat, all_slitcen[ss][dd][mask]])
            all_slitcen_ordr = np.concatenate([all_slitcen_ordr, all_ordrcen_t[mask]])
            all_slitcen_det = np.concatenate([all_slitcen_det, np.full_like(all_ordrcen_t[mask], dd, dtype=int)])
        # Transform the slitcen coordinates for det2 as well
        all_slitcen_spec_t, all_slitcen_spat_t = transform_det2(all_slitcen_spec, all_slitcen_spat, all_slitcen_det, dx, dy, theta)

        # Perform a 2D fit to the spatial coordinates as well
        pypeitFit = fitting.robust_fit(all_slitcen_spec_t/(det_nspec-1), all_slitcen_spat_t/(2*det_nspat), (nspec_coeff, norder_coeff), x2=all_slitcen_ordr/200,
                                       function="legendre2d", maxiter=100, lower=10, upper=10, minx=0, maxx=1,
                                       minx2=0, maxx2=1, use_mad=True, sticky=False)
        # Evaluate the spatial fit at the slitcen positions, and estimate the spatial residuals next the locations of the arc lines
        slitcen_spat_fit = (2*det_nspat) * pypeitFit.eval(all_specpix_t[wgd]/(det_nspec-1), x2=all_orders[wgd]/200)
        spat_resid = all_spatpix_t[wgd] - slitcen_spat_fit
        all_resid[wgd] = np.sqrt(spec_resid**2 + spat_resid**2)
    # residuals = all_wave - all_wavefit
    rms = np.sqrt(np.mean(all_resid**2))
    # print("Fit RMS = ", rms)
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

def fit_mosaic_parameters(all_specpix, all_spatpix, all_wave,
                          all_orders, all_detector, all_setup, all_slitcen, all_ordrcen,
                           dx_range=(-200, 200), dy_range=(-50, 50),
                           theta_range=(-0.01, 0.01),
                           coarse_steps=(21, 15, 11),
                           n_bootstrap=200, pixel_sigma=0.1,
                           do_coarse=False, verbose=True):
    """
    Find the optimal rigid-body alignment (dx, dy, theta) of detector 2
    relative to detector 1 by minimising the RMS of a global 2D polynomial
    wavelength fit to both detectors pooled together.

    Uses a coarse brute-force grid search followed by a fine L-BFGS-B
    optimisation. Parameter uncertainties are estimated by bootstrap
    resampling with Gaussian pixel-coordinate perturbations.

    Parameters
    ----------
    all_specpix : np.ndarray
        Spectral pixel coordinates for all points (both detectors, all setups).
    all_spatpix : np.ndarray
        Spatial pixel coordinates for all points (both detectors, all setups).
    all_wave : np.ndarray
        Wavelengths for all points (both detectors, all setups).
    all_orders : np.ndarray
        Echelle order numbers for all points (both detectors, all setups).
    all_detector : np.ndarray
        Detector index (0 or 1) for all points.
    all_setup : np.ndarray
        Setup index for all points (0, 1, 2, or 3).
    all_slitcen : list of list of np.ndarray
        Slit center coordinates for each setup and detector, as returned by load_data().
    all_ordrcen : list of list of np.ndarray
        Order center numbers for each setup and detector, as returned by load_data().
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
    w1 = (all_detector == 0)
    w2 = (all_detector == 1)
    if verbose:
        print(f"Loaded {len(all_specpix[w1])} points from detector 1, {len(all_specpix[w2])} from detector 2.")
        print(f"Det1 spat range: {all_spatpix[w1].min():.1f} – {all_spatpix[w1].max():.1f}")
        print(f"Det2 spat range: {all_spatpix[w2].min():.1f} – {all_spatpix[w2].max():.1f}")
        print(f"Det1 wave range: {all_wave[w1].min():.1f} – {all_wave[w1].max():.1f}")
        print(f"Det2 wave range: {all_wave[w2].min():.1f} – {all_wave[w2].max():.1f}")

    args = (all_specpix, all_spatpix, all_wave, all_orders, all_detector, all_setup, all_slitcen, all_ordrcen)
    bounds = [dx_range, dy_range, theta_range]

    # --- Stage 1: coarse brute-force search ---
    if do_coarse:
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

    else:
        dx0, dy0, theta0 = 70.0, 0.0, 0.0

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
    if False:
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
        # 'dx_err':            errors['dx'],
        # 'dy_err':            errors['dy'],
        # 'theta_err':         errors['theta'],
        'rms':               rms_fit,
        # 'bootstrap_samples': samples,
        # 'coarse_result':     coarse_result,
        'fine_result':       fine_result,
    }


# ---------------------------------------------------------------------------
# Data loader
# ---------------------------------------------------------------------------

def load_data(idx=None):
    dirc = 'PypeIt_Templates/'
    all_setups = [['564_l', '564_u'],
                  ['580_l', '580_u'],
                  ['760_l', '760_u'],
                  ['860_l', '860_u']]
    if idx is None:
        setups = all_setups
        setup_idx = np.arange(len(all_setups))
    elif 0 <= idx < len(all_setups):
        setups = [all_setups[idx]]
        setup_idx = np.array([idx])
    else:
        raise ValueError(f"Invalid idx={idx}. Must be in range [0, {len(all_setups)-1}] or None.")

    all_slitcen   = [[None for _ in range(len(setups[0]))] for _ in setups]
    all_ordrcen = [[None for _ in range(len(setups[0]))] for _ in setups]
    all_specpix = np.array([])
    all_spatpix = np.array([])
    all_wave = np.array([])
    all_orders = np.array([])
    all_detector = np.array([], dtype=int)
    all_setup = np.array([], dtype=int)

    for ii in range(len(setups)):
        for dd in range(2):
            grat, chip = setups[ii][dd].split('_')
            edges = edgetrace.EdgeTraceSet.from_file(f'{grat}/Edges_{chip}.fits.gz')
            offset = 0.0
            if dd == 0:
                next_offset = edges.traceimg.image.shape[1]
                assert next_offset == det_nspat, f"Expected {det_nspat} spatial pixels per detector, but got {next_offset} from the edges frame."
            else:
                offset = next_offset

            edge_fit = edges.edge_fit
            slitcen = np.mean(edge_fit.reshape((edge_fit.shape[0], edge_fit.shape[1] // 2, 2)), axis=2)
            nspec, norders = slitcen.shape

            wvcal_file = f"{dirc}/WaveCalib_{grat}{chip}.fits"
            wvcal = wavecalib.WaveCalib.from_file(wvcal_file)
            # The order numbers for the 580u and 860u setups were not stored in the WaveCalib files,
            # so we hardcode them here based on the values in the vlt_uves.py spectrograph class.
            # For all other setups, we use the orders from the WaveCalib file.
            if (grat== '580') and (chip == 'u'):
                orders = np.array([105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90], dtype=int)
            elif (grat == '860') and (chip == 'u'):
                orders = np.array([70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58], dtype=int)
            else:
                orders = wvcal.ech_orders
            # Check that the number of wavelength solution orders matches the number of orders in the edges frame
            assert orders.size == norders

            all_slitcen[ii][dd] = slitcen.copy() + offset
            all_ordrcen[ii][dd] = orders.copy()

            # For each wavelength detection, we need to store the wavelength, spec pixel, spat pixel, and order number.
            # We also apply the offset to the spat pixel so that det2 coordinates are in the same frame as det1.
            det_specpix = wvcal.wv_fit2d[0].xval * (nspec-1)
            det_spatpix = np.zeros_like(det_specpix)
            det_wave = wvcal.wv_fit2d[0].yval / wvcal.wv_fit2d[0].x2
            det_orders = wvcal.wv_fit2d[0].x2
            det_detector = np.full_like(det_specpix, dd, dtype=int)
            det_setup = np.full_like(det_specpix, setup_idx[ii], dtype=int)

            # Now determine the spatial pixel locations
            for dd in range(det_specpix.size):
                specpix = det_specpix[dd]
                order = det_orders[dd]
                # Find the closest order in the edges frame
                ord_idx = np.argmin(np.abs(orders - order))
                # Get the spatial pixel location from the slitcen array
                det_spatpix[dd] = offset + np.interp(specpix, np.arange(nspec), slitcen[:, ord_idx])

            # Store the data for this setup and detector
            all_specpix = np.concatenate([all_specpix, det_specpix])
            all_spatpix = np.concatenate([all_spatpix, det_spatpix])
            all_wave = np.concatenate([all_wave, det_wave])
            all_orders = np.concatenate([all_orders, det_orders])
            all_detector = np.concatenate([all_detector, det_detector])
            all_setup = np.concatenate([all_setup, det_setup])

            if False:
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

            # alldata[ii][dd]   = thisdata.copy()
            # allorders[ii][dd] = orders.copy()

    # return alldata, allorders
    return all_specpix, all_spatpix, all_wave, all_orders, all_detector, all_setup, all_slitcen, all_ordrcen

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    # Fit the mosaic parameters for all setups combined (idx=None) or a single setup (idx=0, 1, 2, or 3)
    setups = [0, 1, 2, 3, None]  # Run for each setup and the combined case
    for idx in setups:
        idxstr = f"idx={idx}" if idx is not None else "all setups combined"
        print(f"\n=== Fitting mosaic parameters for setup: {idxstr} ===")
        all_specpix, all_spatpix, all_wave, all_orders, all_detector, all_setup, all_slitcen, all_ordrcen = load_data(idx=idx)

        result = fit_mosaic_parameters(
            all_specpix, all_spatpix, all_wave, all_orders, all_detector, all_setup, all_slitcen, all_ordrcen,
            dx_range=(50, 150),
            dy_range=(-100, 0),
            theta_range=(0.0, 0.0),
            coarse_steps=(21, 21, 1),
            n_bootstrap=200,
            pixel_sigma=0.1)

        print("\nFinal result:")
        print(f"  dx    = {result['dx']:.4f} pixels")
        print(f"  dy    = {result['dy']:.4f} pixels")
        print(f"  theta = {np.degrees(result['theta']):.6f} degrees")
        # print(f"  dx    = {result['dx']:.4f} ± {result['dx_err']:.4f} pixels")
        # print(f"  dy    = {result['dy']:.4f} ± {result['dy_err']:.4f} pixels")
        # print(f"  theta = {np.degrees(result['theta']):.6f} ± "
        #       f"{np.degrees(result['theta_err']):.6f} degrees")
        print(f"  RMS   = {result['rms']:.6f} (pixels)")

        # if idx is not None:
        #     plot_mosaic(alldata, allorders, result=result)
