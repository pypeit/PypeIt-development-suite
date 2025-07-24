from astropy import units as u
from pypeit.coadd3d import DARcorrection
import numpy as np
from matplotlib import pyplot as plt


wave_max = 13000.0 # acquisition in J-band
parangle = 135.0 # degrees
temperature = 0.0 # Celsius
pressure = 600 # mbar 
humidity = 30.0 # percent
dec = 70.0
cosdec = np.cos(np.radians(dec))
waves = np.linspace(9500.0, wave_max, 100)

airmass = [1.7, 1.8, 1.9, 2.0]
# Make this a multi-panel plot
fig = plt.figure(figsize=(10, 10))
ax = fig.subplots(2, 1)

for ix, (wave_ref, band) in enumerate(zip([10500.0, 12500.0], ['Y', 'J'])):
    for am in airmass:
        darcorr = DARcorrection(am, parangle, pressure, temperature, humidity, cosdec, wave_ref=wave_ref)
        dispersion = darcorr.calculate_dispersion(waves)
        disp_arcsec = (dispersion*u.deg).to(u.arcsec).value
        ax[ix].plot(waves, disp_arcsec, label='airmass={:5.3f}'.format(am))
        
    ax[ix].set_xlabel('Wavelength (Angstroms)')
    ax[ix].set_ylabel('DAR offset (arcsec)')
    ax[ix].axvline(wave_ref, color='k', linestyle=':', label='wave_ref')
    ax[ix].set_title('Acquisition in {:s}-band at wave_ref={:5.1f}'.format(band, wave_ref))
    ax[ix].legend()

plt.show()



