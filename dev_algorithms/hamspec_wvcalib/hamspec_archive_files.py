
# imports
import os
import numpy as np
from importlib.resources import files

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.io import fits
from astropy import units

from pypeit.core.wavecal import templates
from pypeit.core.fitting import robust_fit
from pypeit.core import wave as core_wave 

from IPython import embed


# NIRSPEC files
nirspec_angle_fits_file = os.path.join(files('pypeit'), 'data', 'arc_lines', 'reid_arxiv', 'keck_nirspec_j_preupgrade_angle_fits.fits')
nirspec_composite_arc_file = os.path.join(files('pypeit'), 'data', 'arc_lines', 'reid_arxiv', 'keck_nirspec_j_preupgrade_composite_arc.fits')

# Load XIDL solution
hamspec_idl_file = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms', 'hamspec_wvcalib', 'Arc_01_fit.idl')
order_vec, wave, spec, NORDs = templates.xidl_esihires(hamspec_idl_file, ret_NORD=True)
norders = order_vec.size

# Check
if True:
    #idx = 2  # Order 71, which is legit
    # This is the same as -16 in my extracted arc.
    idx = 20

    print(f'Order: {order_vec[idx]}')
    fig = plt.figure(figsize=(12, 8))
    ax = plt.gca()
    ax.plot(core_wave.vactoair(wave[idx, :]*units.AA), spec[idx, :])
    plt.show()
    embed(header='hamspec_archive_files.py: Check 35')

#  Angles file
hdu = fits.open(nirspec_angle_fits_file)
hamspec_angle_params = Table(hdu[1].data)
hamspec_xd_angle_coeffs = hdu[3].data.copy()


## Fill in XD coeffs
hamspec_xd_angle_coeffs[:] = [order_vec[0], 0., 0.]

## Fit for ECH coeffs
n_final = 4
coeff_fit_order_max = 0
hamspec_ech_angle_coeffs = np.zeros((norders, n_final + 1, coeff_fit_order_max + 1))
sigrej = 3.
maxrej = 1
ech_min, ech_max = -1., 1.

nspec = wave.shape[1]
xnspecmin1 = float(nspec - 1)
xnspec = np.arange(nspec)/xnspecmin1

for iord in range(norders):
    # Fit the wavelengths
    pypeitFit = robust_fit(xnspec, wave[iord, :], n_final,
        function='legendre',
        minx=hamspec_angle_params['wave_xmin'], 
        maxx=hamspec_angle_params['wave_xmax'], maxiter=25,
        lower=sigrej, upper=sigrej, maxrej=maxrej, sticky=True, use_mad=True,
        weights=None)

    # Fill in
    hamspec_ech_angle_coeffs[iord, :, 0] = pypeitFit.fitc
    #embed(header='59 of hamspec_archive_files.py')

        #pypeitFit = fitting.robust_fit(
        #    [0.], coeff_this_order[:, ic], coeff_fit_order_vec[ic], function=func,
        #    minx=ech_min, maxx=ech_max, maxiter=25,
        #    lower=sigrej, upper=sigrej, maxrej=maxrej, sticky=True, use_mad=True,
        #    weights=None)
#embed(header='68 of hamspec_archive_files.py')

## Fill in Table
hamspec_angle_params['norders'] = norders
hamspec_angle_params['order_min'] = order_vec[0]
hamspec_angle_params['order_max'] = order_vec[-1]
hamspec_angle_params.remove_column('xdisp_vec')
hamspec_angle_params['xdisp_vec'] = 'Hamspec'
hamspec_angle_params['ech_n_final'] = n_final
hamspec_angle_params.remove_column('ech_coeff_fit_order')

t = np.ones_like(order_vec,dtype=int)*(coeff_fit_order_max+1)
t = t.tolist()

## Remake it because of painful formatting; thanks JFH
items = []
names = []
for key in hamspec_angle_params.keys():
    items.append([hamspec_angle_params[key][0]])
    names.append(key)
items.append([t])
names.append('ech_coeff_fit_order')
hamspec_angle_params=Table(items, names=names)


## Write
outfile = os.path.join(files('pypeit'), 'data', 'arc_lines', 'reid_arxiv', 'lick_hamspec_angle_fits.fits')
hdulist = fits.HDUList()
hdulist.append(fits.BinTableHDU(hamspec_angle_params))  # hdu = 1
hdulist.append(fits.ImageHDU(np.array(hamspec_ech_angle_coeffs)))  # hdu = 2
hdulist.append(fits.ImageHDU(np.array(hamspec_xd_angle_coeffs)))  # hdu = 3
hdulist.writeto(outfile, overwrite=True)
print(f'Wrote: {outfile}')


# Composite arc

hdu = fits.open(nirspec_composite_arc_file)
hamspec_composite_arc_params = Table(hdu[1].data)
hamspec_wave_composite = wave.T
hamspec_arc_composite = spec.T
hamspec_gpm_composite = hamspec_arc_composite.astype(bool)

outfile2 = os.path.join(files('pypeit'), 'data', 'arc_lines', 'reid_arxiv', 'lick_hamspec_composite_arc.fits')
hdulist = fits.HDUList()
hdulist.append(fits.BinTableHDU(hamspec_composite_arc_params)) # hdu = 1
hdulist.append(fits.ImageHDU(np.array(hamspec_wave_composite)))  # hdu = 2
hdulist.append(fits.ImageHDU(np.array(hamspec_arc_composite)))  # hdu = 3
hdulist.append(fits.ImageHDU(np.array(hamspec_gpm_composite.astype(float))))  # hdu = 3
hdulist.writeto(outfile2, overwrite=True)
print(f'Wrote: {outfile2}')