import numpy as np
from astropy import table, units
from astropy.io import fits
from pypeit.core.wave import airtovac
from matplotlib import pyplot as plt
import argparse
from IPython import embed

parser = argparse.ArgumentParser(description='Plot Murpy ThAr calibrated spectrum with lines. '
                                             'It is recommended to set the wavelength range with '
                                             '--wavemin and --wavemax.')
parser.add_argument('--wavemin', default=None, type=float, help='minimum wavelength to plot')
parser.add_argument('--wavemax', default=None, type=float, help='maximum wavelength to plot')
args = parser.parse_args()

linelist_file = 'hires_thar.lst'
calib_spec_file = 'thar_spec_MM201006.fits'

# load Murpy ThAr line list(from xidl/Spec/Arcs/Lists/hires_thar.lst)
linelist = table.Table.read(linelist_file, format='ascii')
linelist.rename_columns(('col1', 'col3'), ('wave', 'ion'))
# convert line list air wavelength to vacuum
linelist['wave'] = airtovac(linelist['wave'] * units.AA).value

# load calibrated spectrum (from xidl/Keck/HIRES/Redux/pro/Arcs/ThAr/thar_spec_MM201006.fits)
calib_spec = fits.open(calib_spec_file)
hdu = fits.open(calib_spec_file)
hdr = hdu[0].header
flux_all = hdu[0].data[0]
wave_all_air = 10**(hdr['CRVAL1'] + (hdr['CDELT1']*np.arange(flux_all.size)))
# convert air wavelength to vacuum
wave_all = airtovac(wave_all_air * units.AA).value

# cut down the wavelength range
gpm = np.ones(flux_all.size, dtype=bool)
if args.wavemin is not None:
    gpm &= wave_all > args.wavemin
if args.wavemax is not None:
    gpm &= wave_all < args.wavemax

if not gpm.any():
    print('No data in the wavelength range specified.')
    exit()

wave = wave_all[gpm]
flux = flux_all[gpm]

# plot
fig = plt.figure(figsize=(23, 6.))
ax = plt.subplot()
ax.minorticks_on()
# spectrum
ax.plot(wave, flux, lw=1, zorder=1)
# lines
line_num = 0
# text offset
toff = np.median(np.diff(wave)) * 2
for i in range(len(linelist)):
    if wave.min() < linelist['wave'][i] < wave.max():
        line_num += 1
        yannot = 0.83 if (line_num % 2) == 0 else 0.65
        ax.axvline(linelist['wave'][i], color='Gray', zorder=-1, alpha=0.5, lw=0.8)
        ax.annotate('{}  {}'.format(linelist['ion'][i], round(linelist['wave'][i],3)),
                    xy=(linelist['wave'][i], 1),
                    xytext=(linelist['wave'][i]+toff, yannot),
                    xycoords=('data', 'axes fraction'),
                    arrowprops=dict(facecolor='None', edgecolor='None', headwidth=0., headlength=0, width=0, shrink=0.),
                    annotation_clip=True, horizontalalignment='center', color='GRAY', fontsize=9, rotation=-90)
ax.set_xlim(wave.min(), wave.max())
ax.set_ylim(-1000, flux.max()*1.1)
ax.set_xlabel('Wavelength  (Angstrom)')

plt.tight_layout()

plt.show()


