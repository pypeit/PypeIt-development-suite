""" Code for building the reid arxiv"""
import glob

from astropy.table import Table

from pypeit import wavecalib

from IPython import embed

def build_reid_arxiv(outfile:str='vlt_xshooter_uvb1x1.fits'):
    order_files = glob.glob('order*_wvcalib.fits')
    order_files.sort()

    orders = []
    waves = []
    fluxes = []
    # Loop on em
    for ifile in order_files[::-1]:
        orders.append(int(ifile[5:7]))
        # Load
        wvCalib = wavecalib.WaveCalib.from_file(ifile)
        # Save
        waves.append(wvCalib.wv_fits[0].wave_soln) 
        fluxes.append(wvCalib.arc_spectra[:,0])

    # Write
    tbl = Table()
    tbl['order'] = orders
    tbl['wave'] = waves
    tbl['flux'] = fluxes
    tbl.write(outfile)
    print(f"Wrote: {outfile}")

if __name__ == "__main__":

    build_reid_arxiv()
