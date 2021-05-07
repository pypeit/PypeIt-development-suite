from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np

# Task for you is to read in the standard star reduction from Extract directory and make some plots of the
# of the counts vs wavelength for each order. The relevant tags are wave, fx, var (sqrt var is sigma)

# Extra credit (you probably won't get this)

# Try to create a SpecObjs in echelle format that contains the data from that standard file. Look at echelle in Dev Suite forPypeit for examples. 

def plot_cts(fitsfile):
    # fitsfile can be 'fits.gz'

    data = fits.open(fitsfile)[1].data
    wave = data['wave']
    fx = data['fx']
    var = data['var']
    order = data['order']

    plt.figure(figsize=(10,7))
    for i in range(len(order)):
        iwant = np.argwhere(wave[i] != 0)
        plt.plot(wave[i][iwant], fx[i][iwant], label=order[i])

    plt.xlabel('Wavelength (A)')
    plt.ylabel('Counts')
    plt.legend()
    plt.ylim([-100, np.max(fx)])
    plt.tight_layout()
    plt.show()