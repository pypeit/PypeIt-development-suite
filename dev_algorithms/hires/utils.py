from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np


def plot_cts(fitsfile):
    # fitsfile can be 'fits.gz'

    data = fits.open(fitsfile)[1].data
    wave = data['wave']
    fx = data['fx']
    var = data['var']
    order = data['order']

    plt.figure(figsize=(10 ,7))
    for i in range(len(order)):
        iwant = np.argwhere(wave[i] != 0)
        plt.plot(wave[i][iwant], fx[i][iwant], label=order[i])

    plt.xlabel('Wavelength (A)')
    plt.ylabel('Counts')
    plt.legend()
    plt.ylim([-100, np.max(fx)])
    plt.tight_layout()
    plt.show()



