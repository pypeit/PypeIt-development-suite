"""
Draft of qa_fit_profile_refactor.

Layout uses fig.add_axes() for fine control over every panel position.

Image panels (left half of figure):
  - Shared y-axis; horizontal colorbars placed above each image with
    tick marks and numeric labels on the top face of the colorbar.
  - Axes dimensions computed analytically for square image pixels:
      img_w = img_h * (FIG_H / FIG_W) * (nspat / nspec)
  - imshow called with aspect='auto' to fill the pre-sized axes box.

Diagnostic panels (right half): three panels stacked vertically.
"""

import os
import sys
from math import gcd
sys.path.insert(0, '/Users/westfall/Work/packages/pypeit-main/pypeit')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pypeit.core import spatialprofile_refactor
from pypeit.core.spatialprofile_refactor import _fit_spectrum_and_normalize

# ---------------------------------------------------------------------------
# Helpers replicated from test_spatialprofile.py
# ---------------------------------------------------------------------------

def make_fourier_spectrum(nspec, nmodes, mean=0.0, std=1.0, seed=42):
    rng = np.random.default_rng(seed)
    k_min = int(np.ceil(nspec / 3.0))
    k_max = nspec // 2
    modes = rng.choice(np.arange(k_min, k_max + 1), size=nmodes, replace=False)
    coeffs = np.zeros(nspec, dtype=complex)
    for k in modes:
        if nspec % 2 == 0 and k == nspec // 2:
            coeffs[k] = rng.standard_normal()
        else:
            a = rng.standard_normal() + 1j * rng.standard_normal()
            coeffs[k] = a
            coeffs[nspec - k] = np.conj(a)
    spectrum = np.real(np.fft.ifft(coeffs))
    spectrum -= spectrum.mean()
    rms = spectrum.std()
    if rms > 0.0:
        spectrum /= rms
    return mean + std * spectrum


def make_profile_inputs(nspec=200, nspat=100, fwhm=4.0, flux_level=1000.0,
                        sn_ratio=20.0, trace_offset=0.0, seed=42):
    rng = np.random.default_rng(seed)
    wave = np.linspace(4000.0, 5000.0, nspec)
    waveimg = wave[:, np.newaxis] * np.ones((nspec, nspat))
    spat_img = np.ones((nspec, 1)) * np.arange(nspat, dtype=float)
    trace_in = np.full(nspec, float(nspat // 2))
    trace_offset_arr = np.broadcast_to(
        np.asarray(trace_offset, dtype=float), (nspec,)
    ).copy()
    true_trace = trace_in + trace_offset_arr
    flux_level_arr = np.broadcast_to(
        np.asarray(flux_level, dtype=float), (nspec,)
    ).copy()
    sigma = fwhm / 2.3548
    dspat_true = spat_img - true_trace[:, np.newaxis]
    profile_gauss = np.exp(-0.5 * (dspat_true / sigma) ** 2)
    profile_norm = profile_gauss / profile_gauss.sum(axis=1, keepdims=True)
    true_image = flux_level_arr[:, None] * profile_norm
    noise_sigma = flux_level_arr / (sn_ratio * np.sqrt(nspat))
    image = true_image + rng.normal(scale=noise_sigma[:, None], size=(nspec, nspat))
    ivar = np.broadcast_to((1.0 / noise_sigma ** 2)[:, None], (nspec, nspat)).copy()
    thismask = np.ones((nspec, nspat), dtype=bool)
    inmask = thismask.copy()
    flux = flux_level_arr.copy()
    noise_sigma_1d = flux_level_arr / sn_ratio
    fluxivar = 1.0 / noise_sigma_1d ** 2
    return (image, ivar, waveimg, thismask, spat_img, trace_in,
            wave, flux, fluxivar, inmask, true_trace)


# ---------------------------------------------------------------------------
# QA function
# ---------------------------------------------------------------------------

def _bin_profile(sigma_x, norm_obj, model_vals, xmin=-7.0, xmax=7.0, nsamp=80):
    """Bin normalised-object values and profile model in sigma_x space."""
    bins = np.linspace(xmin, xmax, nsamp + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    y20 = np.full(nsamp, np.nan)
    y50 = np.full(nsamp, np.nan)
    y80 = np.full(nsamp, np.nan)
    ym  = np.full(nsamp, np.nan)
    for i in range(nsamp):
        mask = (sigma_x >= bins[i]) & (sigma_x < bins[i+1])
        n = mask.sum()
        if n >= 4:
            ys = np.sort(norm_obj[mask])
            y20[i] = ys[int(0.20 * (n - 1))]
            y50[i] = ys[int(0.50 * (n - 1))]
            y80[i] = ys[int(0.80 * (n - 1))]
            ym[i]  = np.median(model_vals[mask])
    return centers, y20, y50, y80, ym


def qa_fit_profile_refactor(
    image, ivar, waveimg, thismask, spat_img, trace_in,
    wave, flux, fluxivar, inmask,
    profile_model, xnew, fwhmfit, med_sn2,
    spec_img=None, percentile_sn2=70.0, thisfwhm=4.0,
    obj_string='', outfile=None, l_limit=0.0, r_limit=0.0,
):
    r"""
    QA diagnostic plot for :func:`fit_profile_refactor`.

    Produces a figure with two sections:

    **Left half** — three 2-D image panels (Data, Model, normalised
    Residual) sharing a common y-axis.  Each panel has a horizontal
    colorbar placed above it with tick marks and numeric labels on the
    top face of the bar.  Image pixels are rendered with a 1:1 aspect
    ratio (square pixels).

    **Right half** — three 1-D diagnostic panels stacked vertically:

    1. *Spectrum* — input ``flux`` vs. S/N-optimal extraction using the
       fitted profile; discrepancies flag profile errors.
    2. *Trace & FWHM* — ``xnew − trace_in`` (trace correction) and
       ``fwhmfit`` vs. spectral row.
    3. *Spatial profile* — ``norm_obj_x`` binned in
       :math:`\sigma_x` space with the profile-model overplotted.

    Parameters
    ----------
    image : :class:`numpy.ndarray`
        Sky-subtracted science image, shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    ivar : :class:`numpy.ndarray`
        Inverse variance of ``image``, same shape.
    waveimg : :class:`numpy.ndarray`
        Wavelength image, same shape.
    thismask : :class:`numpy.ndarray`
        Boolean slit mask, same shape.
    spat_img : :class:`numpy.ndarray`
        Spatial-coordinate image, same shape.
    trace_in : :class:`numpy.ndarray`
        Input object trace, shape :math:`(N_{\rm spec},)`.
    wave : :class:`numpy.ndarray`
        Extracted wavelength array, shape :math:`(N_{\rm spec},)`.
    flux : :class:`numpy.ndarray`
        Extracted flux array, shape :math:`(N_{\rm spec},)`.
    fluxivar : :class:`numpy.ndarray`
        Inverse variance of ``flux``, shape :math:`(N_{\rm spec},)`.
    inmask : :class:`numpy.ndarray`, optional
        Additional boolean pixel mask, same shape as ``image``.
    profile_model : :class:`numpy.ndarray`
        Normalised spatial profile returned by :func:`fit_profile_refactor`,
        shape :math:`(N_{\rm spec}, N_{\rm spat})`.
    xnew : :class:`numpy.ndarray`
        Corrected trace returned by :func:`fit_profile_refactor`, shape
        :math:`(N_{\rm spec},)`.
    fwhmfit : :class:`numpy.ndarray`
        Fitted FWHM per spectral pixel, shape :math:`(N_{\rm spec},)`.
    med_sn2 : :obj:`float`
        Median (S/N)^2 returned by :func:`fit_profile_refactor`.
    spec_img : :class:`numpy.ndarray`, optional
        Integer spectral-row image forwarded to
        :func:`_fit_spectrum_and_normalize`.  Default is ``None``.
    percentile_sn2 : :obj:`float`, optional
        Upper percentile used to estimate ``med_sn2``.  Default is 70.
    thisfwhm : :obj:`float`, optional
        Initial FWHM estimate.  Default is 4.0.
    obj_string : :obj:`str`, optional
        Object label for the figure title.  Default is ``''``.
    outfile : :obj:`str`, optional
        If provided, save the figure to this path instead of displaying it.
        Default is ``None``.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        The matplotlib figure object.
    """
    # -----------------------------------------------------------------------
    # Setup
    # -----------------------------------------------------------------------
    totmask = (ivar > 0.) & thismask
    if inmask is not None:
        totmask &= inmask
    nspec, nspat = image.shape
    row_idx = np.arange(nspec)
    sig2fwhm = np.sqrt(8.0 * np.log(2.0))
    sigma = fwhmfit / sig2fwhm

    # Profile x-axis limits in sigma_x space
    limits_set = (l_limit != 0.0) or (r_limit != 0.0)
    if limits_set:
        prof_center = 0.5 * (l_limit + r_limit)
        prof_half_w = 0.75 * (r_limit - l_limit)
    else:
        prof_center = 0.0
        prof_half_w = 2.5 * sig2fwhm   # total width = 5 × FWHM in sigma_x
    prof_xlim_lo = prof_center - prof_half_w
    prof_xlim_hi = prof_center + prof_half_w

    # -----------------------------------------------------------------------
    # Re-run spectral fitting to obtain norm_obj_x and sigma_x
    # -----------------------------------------------------------------------
    success, _, norm_obj_x, norm_ivar_x, _, _, spec_x = _fit_spectrum_and_normalize(
        wave=wave, flux=flux, fluxivar=fluxivar,
        waveimg=waveimg, image=image, ivar=ivar, totmask=totmask,
        spec_img=spec_img, percentile_sn2=percentile_sn2, fwhm=thisfwhm,
    )

    # -----------------------------------------------------------------------
    # S/N-optimal extraction using the fitted profile
    # -----------------------------------------------------------------------
    opt_num  = np.sum(image * profile_model * ivar, axis=1)
    opt_den  = np.sum(profile_model ** 2 * ivar, axis=1)
    opt_flux = np.where(opt_den > 0, opt_num / opt_den, 0.0)

    # -----------------------------------------------------------------------
    # 2D display arrays
    # -----------------------------------------------------------------------
    model_2d    = profile_model * opt_flux[:, None]
    display_img = np.where(thismask, image,           np.nan)
    display_mod = np.where(thismask, model_2d,        np.nan)
    resid_2d    = np.where(thismask, image - model_2d, np.nan)

    vlo  = np.nanpercentile(display_img,  5)
    vhi  = np.nanpercentile(display_img, 95)
    vres = np.nanpercentile(np.abs(resid_2d[np.isfinite(resid_2d)]), 97)

    # -----------------------------------------------------------------------
    # Spatial coordinate and collapsed profile
    # -----------------------------------------------------------------------
    if success:
        good_x    = norm_ivar_x > 0
        dspat_x   = spat_img[totmask] - trace_in[spec_x]
        sigma_x   = dspat_x / sigma[spec_x]
        profile_x = profile_model[totmask]
        cx, y20, y50, y80, ym = _bin_profile(
            sigma_x[good_x], norm_obj_x[good_x], profile_x[good_x],
            xmin=prof_xlim_lo, xmax=prof_xlim_hi, nsamp=80,
        )

    # -----------------------------------------------------------------------
    # Figure geometry (figure-fraction coordinates [0, 1])
    # -----------------------------------------------------------------------
    FIG_W, FIG_H = 16.0, 8.0
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    sn_str   = f'S/N = {np.sqrt(med_sn2):.1f}'
    fwhm_str = f'FWHM = {np.median(fwhmfit):.2f} px'
    dtrc_str = f'|Δtrace|max = {np.max(np.abs(xnew - trace_in)):.2f} px'

    # --- Image panel aspect ratio ---
    # Choose an integer panel ratio R_panel (height : width of axes box) in [1, 10].
    # R_panel = nspec/nspat gives square pixels; larger R_panel compresses the spectral axis.
    R_natural = nspec / nspat
    R_panel   = int(np.clip(round(R_natural), 1, 10))

    # Pixel display ratio x:y (spatial : spectral) at the chosen R_panel
    _g    = gcd(nspec, R_panel * nspat)
    pix_p = nspec // _g          # x (spatial) component
    pix_q = (R_panel * nspat) // _g  # y (spectral) component
    pix_str = f'pix x:y = {pix_p}:{pix_q}'

    # --- Horizontal layout: right edge of residual image panel fixed at img_right_max ---
    img_left_start = 0.06
    img_right_max  = 0.45
    img_w          = (img_right_max - img_left_start) / 3

    # --- Vertical layout: image height from aspect ratio, centered in available space ---
    cbar_h     = 0.010   # 1/3 of original
    cbar_pad   = 0.008
    bot_margin = 0.055   # room for x-axis labels
    top_margin = 0.085   # room for colorbar + label + suptitle
    avail_h    = 1.0 - bot_margin - cbar_pad - cbar_h - top_margin
    img_h      = min(img_w * (FIG_W / FIG_H) * R_panel, avail_h)
    img_bot    = bot_margin + (avail_h - img_h) / 2
    cbar_bot   = img_bot + img_h + cbar_pad

    img_l = [img_left_start + i * img_w for i in range(3)]

    # Colorbar geometry — 75 % of image width, centered above each panel
    cbar_w   = 0.75 * img_w
    cbar_off = 0.125 * img_w

    # Spatial profile panel — right third, aligned with image panels
    # Profile + residual panels together span the full image height (img_bot to img_bot + img_h).
    prof_r         = 0.97
    prof_w         = 0.29 * 0.75   # reduced by 25 %
    prof_l         = prof_r - prof_w
    resid_prof_h   = 0.2 * img_h
    prof_h         = img_h - resid_prof_h
    resid_prof_bot = img_bot
    prof_bot       = img_bot + resid_prof_h

    # Spectrum + spectrum-residual panels — rotated (row on y), same height as images
    # spec_w + spec_resid_w + spec_gap = prof_l - spec_l
    # spec_resid_w = spec_w / 2  →  spec_w * 3/2 + spec_gap = prof_l - spec_l
    img_right    = img_l[2] + img_w
    spec_l       = img_right
    spec_gap     = 0.045                           # gap to profile panel
    spec_w       = (prof_l - spec_gap - spec_l) * 2.0 / 3.0
    spec_resid_l = spec_l + spec_w
    spec_resid_w = spec_w / 2.0

    # Object string: centered between spectrum left and profile right,
    # halfway between the top of the panels and the top of the figure.
    title_x = 0.5 * (spec_l + prof_r)
    title_y = cbar_bot + cbar_h / 2 + 0.04
    fig.text(title_x, title_y,
             f'{obj_string}     [{sn_str},  {fwhm_str},  {dtrc_str}]',
             ha='center', va='center', fontsize=10)

    # -----------------------------------------------------------------------
    # Image axes (shared y-axis)
    # -----------------------------------------------------------------------
    ax_data  = fig.add_axes([img_l[0], img_bot, img_w, img_h])
    ax_model = fig.add_axes([img_l[1], img_bot, img_w, img_h], sharey=ax_data)
    ax_resid = fig.add_axes([img_l[2], img_bot, img_w, img_h], sharey=ax_data)
    ax_model.tick_params(labelleft=False)
    ax_resid.tick_params(labelleft=False)

    # Colorbar axes
    cax_data  = fig.add_axes([img_l[0] + cbar_off, cbar_bot, cbar_w, cbar_h])
    cax_model = fig.add_axes([img_l[1] + cbar_off, cbar_bot, cbar_w, cbar_h])
    cax_resid = fig.add_axes([img_l[2] + cbar_off, cbar_bot, cbar_w, cbar_h])

    # Images — viridis for data/model, RdBu_r for residual
    extent   = [0, nspat, 0, nspec]
    im_data  = ax_data.imshow(display_img, origin='lower', aspect='auto',
                               vmin=vlo,   vmax=vhi,  cmap='viridis', extent=extent)
    im_model = ax_model.imshow(display_mod, origin='lower', aspect='auto',
                                vmin=vlo,   vmax=vhi,  cmap='viridis', extent=extent)
    im_resid = ax_resid.imshow(resid_2d,    origin='lower', aspect='auto',
                                vmin=-vres, vmax=vres, cmap='RdBu_r',  extent=extent)

    # Pixel aspect-ratio label: text box with white background inside model panel
    ax_model.text(
        0.04, 0.03, pix_str,
        transform=ax_model.transAxes, ha='left', va='bottom',
        fontsize=7, color='0.20',
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.80, pad=1.5),
    )

    # Colorbars: ticks and labels on top face
    for im, cax, label in [
        (im_data,  cax_data,  'Data  [Counts]'),
        (im_model, cax_model, 'Model  [Counts]'),
        (im_resid, cax_resid, 'Residual  [Counts]'),
    ]:
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        cbar.set_label(label, fontsize=8)

    # Axis labels and trace overlays
    for ax in [ax_data, ax_model, ax_resid]:
        ax.set_xlabel('Spatial pixel', fontsize=8)
    ax_data.set_ylabel('Spectral row', fontsize=8)

    ax_data.plot(trace_in, row_idx, '--', color='lime', lw=1.0, alpha=0.85)
    ax_data.plot(xnew,     row_idx, '-',  color='red',  lw=1.5, alpha=0.90)
    ax_model.plot(xnew, row_idx, '-', color='red', lw=1.5, alpha=0.90)
    ax_resid.plot(xnew, row_idx, '-', color='k',   lw=0.8, alpha=0.50)

    # FWHM width lines on model panel (trace ± half-FWHM)
    ax_model.plot(xnew - fwhmfit / 2, row_idx, '--', color='black', lw=0.8, alpha=0.7)
    ax_model.plot(xnew + fwhmfit / 2, row_idx, '--', color='black', lw=0.8, alpha=0.7)

    # -----------------------------------------------------------------------
    # Spectrum panel — rotated: flux on x-axis, spectral row on y-axis
    # -----------------------------------------------------------------------
    ax_spec = fig.add_axes([spec_l, img_bot, spec_w, img_h], sharey=ax_data)
    ax_spec.step(flux,     row_idx, color='0.60',   lw=0.8, where='mid',
                 label='input flux',    zorder=2)
    ax_spec.step(opt_flux, row_idx, color='tomato', lw=1.5, where='mid',
                 label='optimal extr.', zorder=3)
    ax_spec.set_xlabel('Counts / row', fontsize=8)
    ax_spec.set_title(f'Spectrum  ({sn_str})', fontsize=9)
    ax_spec.tick_params(labelsize=7, labelleft=False)
    ax_spec.set_axisbelow(True)
    ax_spec.grid(True, color='lightgray', lw=0.5)

    # -----------------------------------------------------------------------
    # Spectrum residual panel — adjoins spectrum panel, step function
    # -----------------------------------------------------------------------
    ax_spec_res = fig.add_axes([spec_resid_l, img_bot, spec_resid_w, img_h],
                               sharey=ax_data)
    spec_res = flux - opt_flux
    ax_spec_res.step(spec_res, row_idx, color='k', lw=0.8, where='mid')
    ax_spec_res.axvline(0, color='0.50', lw=0.5, ls='--')
    ax_spec_res.set_xlabel('Residual', fontsize=8)
    ax_spec_res.tick_params(labelsize=7, labelleft=False)
    ax_spec_res.set_axisbelow(True)
    ax_spec_res.grid(True, color='lightgray', lw=0.5)

    # -----------------------------------------------------------------------
    # Spatial profile panel — right third, ±10 σ
    # -----------------------------------------------------------------------
    ax_prof = fig.add_axes([prof_l, prof_bot, prof_w, prof_h])

    if success:
        rng_sub  = np.random.default_rng(0)
        in10     = good_x & (sigma_x >= prof_xlim_lo) & (sigma_x <= prof_xlim_hi)
        idx_in10 = np.where(in10)[0]
        sub      = rng_sub.choice(idx_in10, size=min(len(idx_in10), 8000), replace=False)
        ax_prof.scatter(sigma_x[sub], norm_obj_x[sub],
                        s=0.4, c='k', alpha=0.12, rasterized=True, zorder=2)
        finite = np.isfinite(y50)
        if finite.any():
            ax_prof.vlines(cx[finite], y20[finite], y80[finite],
                           color='orange', lw=1.0, alpha=0.6, zorder=3)
            ax_prof.plot(cx[finite], y50[finite], 'o', color='lime',
                         ms=2.5, zorder=4, label='data 20/50/80 %ile')
            ax_prof.plot(cx[finite], ym[finite], '-', color='red',
                         lw=1.8, zorder=5, label='profile model')

    ax_prof.axhline(0, color='0.50', lw=0.5)
    if limits_set:
        ax_prof.axvline(l_limit, color='0.40', lw=1.0, ls='--', alpha=0.8)
        ax_prof.axvline(r_limit, color='0.40', lw=1.0, ls='--', alpha=0.8)
    ax_prof.set_xlim(prof_xlim_lo, prof_xlim_hi)
    ax_prof.set_ylabel('Normalised flux', fontsize=9)
    ax_prof.set_title('Spatial profile (spectrally collapsed)', fontsize=9)
    ax_prof.tick_params(labelsize=8, labelbottom=False)
    ax_prof.set_axisbelow(True)
    ax_prof.grid(True, color='lightgray', lw=0.5)
    if success:
        ax_prof.legend(fontsize=8)

    # -----------------------------------------------------------------------
    # Profile residual panel — below spatial profile, shared x-axis, no gap
    # -----------------------------------------------------------------------
    ax_pres = fig.add_axes([prof_l, resid_prof_bot, prof_w, resid_prof_h],
                           sharex=ax_prof)
    if success:
        res_all = norm_obj_x[idx_in10] - profile_x[idx_in10]
        res_sub = norm_obj_x[sub]       - profile_x[sub]
        ax_pres.scatter(sigma_x[sub], res_sub,
                        s=0.4, c='k', alpha=0.12, rasterized=True, zorder=2)
        ax_pres.axhline(0, color='0.50', lw=0.5, ls='--')
        if len(res_all) >= 2:
            p005   = np.nanpercentile(res_all, 0.5)
            p995   = np.nanpercentile(res_all, 99.5)
            half   = p995 - p005
            center = 0.5 * (p005 + p995)
            ax_pres.set_ylim(center - half, center + half)
    ax_pres.set_xlabel(r'$x / \sigma$', fontsize=9)
    ax_pres.set_ylabel('Resid.', fontsize=7)
    ax_pres.tick_params(labelsize=7)
    ax_pres.set_axisbelow(True)
    ax_pres.grid(True, color='lightgray', lw=0.5)

    # -----------------------------------------------------------------------
    # Global: all tick marks point inward
    # -----------------------------------------------------------------------
    for ax in fig.get_axes():
        ax.tick_params(direction='in', which='both')

    # -----------------------------------------------------------------------
    # Save or display
    # -----------------------------------------------------------------------
    if outfile is not None:
        fig.canvas.print_figure(outfile, dpi=140, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# Example — test_fourier_flux_and_curved_trace dataset
# ---------------------------------------------------------------------------

nspec, nspat = 200, 100
fwhm     = 4.0
sn_ratio = 20.0

flux_level   = make_fourier_spectrum(nspec, nmodes=5, mean=1000.0, std=50.0)
t            = np.linspace(-1.0, 1.0, nspec)
trace_offset = 10.0 * (2.0 * t ** 2 - 1.0)

image, ivar, waveimg, thismask, spat_img, _, wave, flux, fluxivar, inmask, true_trace = \
    make_profile_inputs(
        nspec=nspec, nspat=nspat, fwhm=fwhm, sn_ratio=sn_ratio,
        flux_level=flux_level, trace_offset=trace_offset,
    )

print("Running fit_profile_refactor ...")
profile_model, xnew, fwhmfit, med_sn2 = spatialprofile_refactor.fit_profile_refactor(
    image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
    spat_img=spat_img, trace_in=true_trace, wave=wave, flux=flux,
    fluxivar=fluxivar, inmask=inmask, thisfwhm=fwhm, sn_gauss=4.0,
)
print(f"  med_sn2          = {med_sn2:.1f}")
print(f"  median(fwhmfit)  = {np.median(fwhmfit):.3f} px")
print(f"  max|xnew-trace|  = {np.max(np.abs(xnew - true_trace)):.3f} px")

outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'qa_spatprof.png')
print(f"\nGenerating QA figure → {outfile}")
qa_fit_profile_refactor(
    image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
    spat_img=spat_img, trace_in=true_trace, wave=wave, flux=flux,
    fluxivar=fluxivar, inmask=inmask,
    profile_model=profile_model, xnew=xnew, fwhmfit=fwhmfit, med_sn2=med_sn2,
    thisfwhm=fwhm, obj_string='Fourier flux + curved trace',
    outfile=outfile,
)
print("Done.")
