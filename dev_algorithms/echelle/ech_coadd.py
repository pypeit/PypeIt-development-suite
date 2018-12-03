
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from pypeit.core import coadd
from pypeit.core import load
from pypeit import msgs

from linetools.spectra.utils import collate



def ech_load_spec(files,objid=None,norder=None,extract='OPT',flux=True):
    """
    files: A list of file names
    objid:
    norder:
    extract:
    flux:
    """

    nfiles = len(files)
    if objid is None:
        objid = ['OBJ0000'] * nfiles
    elif len(objid) == 1:
        objid = objid * nfiles
    elif len(objid) != nfiles:
        msgs.error('The length of objid should be either 1 or equal to the number of spectra files.')

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):
        if norder is None:
            ext_final = fits.getheader(fname,-1)
            norder = ext_final['ORDER']+1
            msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
        elif norder <=1:
            msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

        # read extension name into a list
        primary_header = fits.getheader(fname,0)
        nspec = primary_header['NSPEC']
        extnames = [primary_header['EXT0001']] * nspec
        for kk in range(nspec):
            extnames[kk] = primary_header['EXT'+'{0:04}'.format(kk+1)]
        extnameroot = extnames[0]

        for iord in range(norder):

            # Figure out which extension is the required data
            ordername = '{0:04}'.format(iord)
            extname = extnameroot.replace('OBJ0000',objid[ii])
            extname = extname.replace('ORDER0000','ORDER'+ordername)
            try:
                exten = extnames.index(extname) + 1
                msgs.info("Loading extension {:s} of spectrum {:s}".format(extname, fname))
            except:
                msgs.error("Spectrum {:s} does not contain {:s} extension".format(fname,extname))

            spectrum = load.load_1dspec(fname, exten=exten, extract=extract, flux=flux)
            # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
            #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
            bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                               np.abs(spectrum.flux) > 1e30,
                               spectrum.sig**2 > 1e10,
                               ], axis=0)
            if np.sum(bad_flux):
                msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
                spectrum.data['flux'][spectrum.select][bad_flux] = 0.
                spectrum.data['sig'][spectrum.select][bad_flux] = 0.
            # Append
            spectra_list.append(spectrum)
    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra

#def ech_flux(files,norder=5):


def ech_coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  wave_grid_min=None, wave_grid_max=None,v_pix=None,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=True, qafile=None, outfile=None,
                  do_cr=True, **kwargs):


    ech_kwargs = {'echelle':True,'wave_grid_min': wave_grid_min, 'wave_grid_max': wave_grid_max, 'v_pix': v_pix}
    kwargs.update(ech_kwargs)
    spec1d = coadd.coadd_spectra(spectra, wave_grid_method=wave_grid_method, niter=niter,
                        scale_method=scale_method, do_offset=do_offset, sigrej_final=sigrej_final,
                        do_var_corr=do_var_corr, qafile=qafile, outfile=outfile,
                        do_cr=do_cr, **kwargs)
    return spec1d



#scifiles = ['/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698_FLUX.fits',
#            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100949.328_FLUX.fits',
#            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T101642.428_FLUX.fits',
#            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T102337.058_FLUX.fits']
#objids = ['OBJ0001','OBJ0002','OBJ0002','OBJ0001']
norder =5
datapath = '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/'
cat = np.genfromtxt(datapath+'J0252_objinfo.txt',dtype=str)
filenames = cat[:,0]
scifiles = []
for i in range(len(filenames)):
    filename = datapath+filenames[i]
    scifiles += [filename.replace('.fits','_FLUX.fits')]
objids = cat[:,1]

#scifiles = ['/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T100254.698_FLUX.fits',
#            '/Users/feige/Work/Observations/NIRES/NIRES_Barth/J0252/reduce0930/Science/spec1d_J0252-0503_NIRES_2018Oct01T102337.058_FLUX.fits']
#objids = ['OBJ0001','OBJ0001']

plt.figure()
for i in range(len(scifiles)):
    sciframe = scifiles[i]
    spectra = ech_load_spec([sciframe],objid = [objids[i]],norder=norder,extract='OPT',flux=True)
    for iord in range(norder-1):
        plt.plot(spectra[iord+1].wavelength,spectra[iord+1].flux)

plt.ylim([-0.5,2.0])
plt.show()

kwargs={}
spectra = ech_load_spec(scifiles,objid=objids,norder=norder,extract='OPT',flux=True)
spec1d = ech_coadd_spectra(spectra, wave_grid_method='velocity', niter=5,
                  wave_grid_min=9400.0, wave_grid_max=None,v_pix=None,
                  scale_method='auto', do_offset=False, sigrej_final=3.,
                  do_var_corr=True, qafile='test', outfile=None,
                  do_cr=True,**kwargs)
plt.figure()
plt.plot(spec1d.wavelength,spec1d.flux)
plt.plot(spec1d.wavelength,spec1d.sig,'-',color='0.7')
plt.show()
from IPython import embed
embed()