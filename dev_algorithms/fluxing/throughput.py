
import numpy as np
from astropy import units as u
from astropy import constants as const
from pypeit import sensfunc
from pypeit.spectrographs import util
from pypeit import utils
from matplotlib import pyplot as plt


def sensfunc_to_thru(sensfile):

    wave, sensfunction, meta_table, out_table, header_sens = sensfunc.SensFunc.load(sensfile)
    spectrograph = util.load_spectrograph(header_sens['PYP_SPEC'])
    sensfunc_units = 1e-17*u.erg/u.cm**2

    sens_gpm = (utils.inverse(sensfunction) > np.finfo(float).tiny) &\
               (utils.inverse(sensfunction) < np.finfo(float).max) & (wave > 1.0)
    inv_wave = utils.inverse(wave[sens_gpm])/u.angstrom
    inv_sensfunc = utils.inverse(sensfunction[sens_gpm])/sensfunc_units
    eff_aperture = spectrograph.telescope['eff_aperture']*u.m**2
    throughput = np.zeros_like(sensfunction)

    thru = ((const.h*const.c)*inv_wave/eff_aperture*inv_sensfunc).decompose()
    throughput[sens_gpm] = thru

    return wave, throughput


sensfile = '/Users/joe/nires_1220/redux/sens_s201226_0027-Feige110_NIRES_2020Dec26T044943.096.fits'
wave, throughput  = sensfunc_to_thru(sensfile)

norders = wave.shape[1]
for iorder in range(norders):
    wv_mask = wave[:, iorder] > 1.0
    plt.plot(wave[wv_mask, iorder],sensfunction[wv_mask, iorder])
#    plt.plot(wave[wv_mask, iorder],throughput[wv_mask, iorder])

plt.ylim((0.0, 100))
plt.ylim((0.0, 0.5))