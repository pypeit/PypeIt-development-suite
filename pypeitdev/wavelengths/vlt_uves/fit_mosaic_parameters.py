import numpy as np

from pypeit.core.arc import fit2darc
from pypeit import wavecalib
from pypeit import edgetrace
from IPython import embed

from scipy.optimize import minimize
from scipy.optimize import brute
from numpy.polynomial.polynomial import polyval2d


# ---------------------------------------------------------------------------
# 2D polynomial wavelength model (no PypeIt dependency)
# ---------------------------------------------------------------------------

def fit_wave2d(spec, spat, wave, order_spec=4, order_spat=4):
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
    # Normalise coordinates to [-1, 1]
    spec_min, spec_max = spec.min(), spec.max()
    spat_min, spat_max = spat.min(), spat.max()

    spec_norm = 2.0 * (spec - spec_min) / (spec_max - spec_min) - 1.0
    spat_norm = 2.0 * (spat - spat_min) / (spat_max - spat_min) - 1.0

    norm_params = dict(spec_min=spec_min, spec_max=spec_max,
                       spat_min=spat_min, spat_max=spat_max)

    # Build Vandermonde matrix for least-squares fit
    # Each row: [spec_norm^i * spat_norm^j] for all (i,j) with i<=order_spec, j<=order_spat
    n_terms = (order_spec + 1) * (order_spat + 1)
    A = np.zeros((len(spec), n_terms))
    col = 0
    for i in range(order_spec + 1):
        for j in range(order_spat + 1):
            A[:, col] = spec_norm**i * spat_norm**j
            col += 1

    # Least-squares solve
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
        Coordinates to evaluate at (in original pixel units).
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

    The transform is: rotate about the centre of detector 2, then translate.
    Since theta is small, the rotation centre choice has a minor effect on dx/dy
    but it is important to be consistent between the coarse and fine stages.

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

    # Rotate about the centroid of detector 2 points
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
    Unpack alldata and allorders into flat arrays for each detector.

    Parameters
    ----------
    alldata : list of list of np.ndarray, shape (nsetups, 2, nspec, norders, 3)
        alldata[setup][detector] is an array of shape (nspec, norders, 3)
        where axis 2 contains (wavelength, spec, spat).
    allorders : list of list of np.ndarray
        allorders[setup][detector] contains the echelle order numbers.

    Returns
    -------
    spec1, spat1, wave1 : np.ndarray
        Flattened arrays for detector 1, pooled across all setups and orders.
    spec2, spat2, wave2 : np.ndarray
        Flattened arrays for detector 2, pooled across all setups and orders.
    """
    spec1_list, spat1_list, wave1_list = [], [], []
    spec2_list, spat2_list, wave2_list = [], [], []

    for ii in range(len(alldata)):
        for det, (spec_list, spat_list, wave_list) in enumerate(
                [(spec1_list, spat1_list, wave1_list),
                 (spec2_list, spat2_list, wave2_list)]):
            data = alldata[ii][det]          # shape (nspec, norders, 3)
            nspec, norders, _ = data.shape
            for ord_idx in range(norders):
                wave_list.append(data[:, ord_idx, 0])
                spec_list.append(data[:, ord_idx, 1])
                spat_list.append(data[:, ord_idx, 2])

    return (np.concatenate(spec1_list), np.concatenate(spat1_list), np.concatenate(wave1_list),
            np.concatenate(spec2_list), np.concatenate(spat2_list), np.concatenate(wave2_list))


# ---------------------------------------------------------------------------
# Cost function
# ---------------------------------------------------------------------------

def cost_function(params, spec1, spat1, wave1, spec2, spat2, wave2,
                  order_spec=6, order_spat=4):
    """
    Cost function for the rigid-body alignment optimisation.

    Applies the transform to detector 2, pools all points with detector 1,
    fits a 2D polynomial wavelength model, and returns the RMS residual.

    Parameters
    ----------
    params : array-like, length 3
        (dx, dy, theta) — spatial shift, spectral shift, rotation in radians.
    spec1, spat1, wave1 : np.ndarray
        Detector 1 data (all orders, all setups, flattened).
    spec2, spat2, wave2 : np.ndarray
        Detector 2 data (all orders, all setups, flattened).
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

    # Pool all data
    spec_all = np.concatenate([spec1, spec2_t])
    spat_all = np.concatenate([spat1, spat2_t])
    wave_all = np.concatenate([wave1, wave2])

    # Fit 2D polynomial and return RMS
    _, _, _, rms = fit_wave2d(spec_all, spat_all, wave_all,
                              order_spec=order_spec, order_spat=order_spat)
    return rms


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------

def fit_mosaic_parameters(alldata, allorders,
                           dx_range=(-200, 200), dy_range=(-10, 10),
                           theta_range=(-0.01, 0.01),
                           coarse_steps=(21, 11, 11),
                           order_spec=4, order_spat=4,
                           verbose=True):
    """
    Find the optimal rigid-body alignment (dx, dy, theta) of detector 2
    relative to detector 1 by minimising the RMS of a global 2D wavelength fit.

    Uses a coarse brute-force grid search followed by a fine gradient-based
    optimisation (L-BFGS-B).

    Parameters
    ----------
    alldata : list of list of np.ndarray
        As returned by load_data(), with shape (nsetups, 2, nspec, norders, 3).
    allorders : list of list of np.ndarray
        As returned by load_data().
    dx_range : tuple
        (min, max) search range for the spatial shift in pixels.
    dy_range : tuple
        (min, max) search range for the spectral shift in pixels.
    theta_range : tuple
        (min, max) search range for the rotation in radians.
    coarse_steps : tuple of int
        Number of grid points for (dx, dy, theta) in the brute-force search.
    order_spec, order_spat : int
        Polynomial degrees for the 2D wavelength fit.
    verbose : bool
        Print progress and results.

    Returns
    -------
    result : dict with keys:
        'dx'    : float, optimal spatial shift in pixels
        'dy'    : float, optimal spectral shift in pixels
        'theta' : float, optimal rotation in radians
        'rms'   : float, RMS wavelength residual at the optimum
        'coarse_result' : output of scipy.optimize.brute
        'fine_result'   : output of scipy.optimize.minimize
    """
    # Unpack data into flat arrays
    spec1, spat1, wave1, spec2, spat2, wave2 = unpack_data(alldata, allorders)

    if verbose:
        print(f"Loaded {len(spec1)} points from detector 1, {len(spec2)} from detector 2.")

    args = (spec1, spat1, wave1, spec2, spat2, wave2, order_spec, order_spat)

    # --- Stage 1: coarse brute-force search ---
    if verbose:
        print("Running coarse grid search...")

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
                           method='L-BFGS-B',
                           bounds=[dx_range, dy_range, theta_range],
                           options={'ftol': 1e-12, 'gtol': 1e-8, 'maxiter': 1000})

    dx_fit, dy_fit, theta_fit = fine_result.x
    rms_fit = fine_result.fun

    if verbose:
        print(f"Fine result:   dx={dx_fit:.4f} pix, dy={dy_fit:.4f} pix, "
              f"theta={np.degrees(theta_fit):.6f} deg, RMS={rms_fit:.6f}")
        print(f"Optimisation {'succeeded' if fine_result.success else 'FAILED: ' + fine_result.message}")

    return {
        'dx': dx_fit,
        'dy': dy_fit,
        'theta': theta_fit,
        'rms': rms_fit,
        'coarse_result': coarse_result,
        'fine_result': fine_result,
    }

def load_data():
    dirc = 'PypeIt_Templates/'
    # setups = [['564_l', '564_u'],
    #          ['580_l', '580_u']]
    setups = [['580_l', '580_u']]
    alldata = [[None for dd in range(len(setups[0]))] for all in setups]
    allorders = [[None for dd in range(len(setups[0]))] for all in setups]
    for ii in range(len(setups)):
        # Load the data for each detector
        for dd in range(2):
            grat, chip = setups[ii][dd].split('_')
            # Load the edges
            edges = edgetrace.EdgeTraceSet.from_file(f'{grat}/Edges_{chip}.fits.gz')
            offset = 0.0
            if dd == 0:
                next_offset = edges.traceimg.image.shape[1]
            else:
                offset = next_offset
            edge_fit = edges.edge_fit
            slitcen = np.mean(edge_fit.reshape((edge_fit.shape[0], edge_fit.shape[1] // 2, 2)), axis=2)
            nspec, norders = slitcen.shape
            wvcal_file = f"{dirc}/WaveCalib_{grat}{chip}.fits"
            wvcal = wavecalib.WaveCalib.from_file(wvcal_file)
            orders = wvcal.ech_orders
            assert orders.size==norders
            thisdata = np.zeros((nspec, norders, 3))
            for ord in range(norders):
                wave = wvcal.wv_fit2d[0].eval(np.linspace(0.0, 1.0, nspec), x2=np.full(nspec, orders[ord]))/orders[ord]
                spec = np.arange(nspec)
                spat = slitcen[:,ord] + offset
                thisdata[:, ord, 0] = wave
                thisdata[:, ord, 1] = spec
                thisdata[:, ord, 2] = spat
            alldata[ii][dd] = thisdata.copy()
            allorders[ii][dd] = orders.copy()
    return alldata, allorders


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    # Load the data
    alldata, allorders = load_data()

    # spec1, spat1, wave1, spec2, spat2, wave2 = unpack_data(alldata, allorders)
    # print("Det1 spat range:", spat1.min(), spat1.max())
    # print("Det2 spat range:", spat2.min(), spat2.max())
    # print("Det1 wave range:", wave1.min(), wave1.max())
    # print("Det2 wave range:", wave2.min(), wave2.max())

    result = fit_mosaic_parameters(alldata, allorders,
                                   dx_range=(-200, 200),
                                   dy_range=(-10, 10),
                                   theta_range=(-0.01, 0.01),
                                   coarse_steps=(21, 11, 11))

    print("\nFinal result:")
    print(f"  dx    = {result['dx']:.4f} pixels")
    print(f"  dy    = {result['dy']:.4f} pixels")
    print(f"  theta = {np.degrees(result['theta']):.6f} degrees")
    print(f"  RMS   = {result['rms']:.6f} (wavelength units)")
