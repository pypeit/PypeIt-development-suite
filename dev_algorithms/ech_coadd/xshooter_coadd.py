import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import coadd1d
from pypeit import inputfiles
from pypeit.par import pypeitpar
from pypeit.scripts import scriptbase
from pypeit.spectrographs.util import load_spectrograph
from pypeit import sensfunc
from pypeit.core import coadd, flux_calib




coadd1d_file_vis = './vlt_xshooter_vis.coadd1d'
coadd1d_file_nir = './vlt_xshooter_nir.coadd1d'


# Load the file
# config_lines, spec1dfiles, objids = read_coaddfile(args.coadd1d_file)
coadd1dFile_vis = inputfiles.Coadd1DFile.from_file(coadd1d_file_vis)
coadd1dFile_nir = inputfiles.Coadd1DFile.from_file(coadd1d_file_nir)


# Read in spectrograph from spec1dfile header
header_vis = fits.getheader(coadd1dFile_vis.filenames[0])
spectrograph_vis = load_spectrograph(header_vis['PYP_SPEC'])

header_nir = fits.getheader(coadd1dFile_nir.filenames[0])
spectrograph_nir = load_spectrograph(header_nir['PYP_SPEC'])

# Parameters
spectrograph_vis_def_par = spectrograph_vis.default_pypeit_par()
par_vis = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_vis_def_par.to_config(),
                                             merge_with=coadd1dFile_vis.cfg_lines)
par_nir = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_nir.default_pypeit_par().to_config(),
                                             merge_with=coadd1dFile_nir.cfg_lines)

# Write the par to disk
#print("Writing the parameters to {}".format(par_outfile))
#par_nir.to_config(args.par_outfile)
sensfile_vis = par_vis['coadd1d']['sensfuncfile']
sensfile_nir = par_nir['coadd1d']['sensfuncfile']
coaddfile = 'test_echelle.fits'

# Instantiate the coadd1d object
coadd1d_vis = coadd1d.CoAdd1D.get_instance(coadd1dFile_vis.filenames, coadd1dFile_vis.objids,
                                           spectrograph=spectrograph_vis, par=par_vis['coadd1d'],
                                           sensfile=sensfile_vis, debug=True, show=True)
coadd1d_nir = coadd1d.CoAdd1D.get_instance(coadd1dFile_nir.filenames, coadd1dFile_nir.objids,
                                           spectrograph=spectrograph_nir, par=par_nir['coadd1d'],
                                           sensfile=sensfile_nir, debug=True, show=True)

waves_vis, fluxes_vis, ivars_vis, gpms_vis, header_out_vis = coadd1d_vis.load_arrays()
waves_nir, fluxes_nir, ivars_nir, gpms_nir, header_out_nir = coadd1d_nir.load_arrays()

weights_sens_vis = sensfunc.SensFunc.sensfunc_weights(sensfile_vis, waves_vis)
weights_sens_nir = sensfunc.SensFunc.sensfunc_weights(sensfile_nir, waves_nir)


shape_vis = waves_vis.shape
shape_nir = waves_nir.shape
nspec_vis, norders_vis, nexp_vis = shape_vis
nspec_nir, norders_nir, nexp_nir = shape_nir
nspec_out = np.fmax(nspec_vis, nspec_nir)
norders_out = norders_vis + norders_nir
nexp = nexp_vis # TODO Generalize to differing numbers of exposures???

shape_out = (nspec_out, norders_out, nexp)
waves = np.zeros(shape_out)
fluxes = np.zeros(shape_out)
ivars = np.zeros(shape_out)
gpms = np.zeros(shape_out, dtype=bool)
weights_sens = np.zeros(shape_out)

waves[:nspec_vis, :norders_vis, :] = waves_vis
fluxes[:nspec_vis, :norders_vis, :] = fluxes_vis
ivars[:nspec_vis, :norders_vis, :] = ivars_vis
gpms[:nspec_vis, :norders_vis, :] = gpms_vis
weights_sens[:nspec_vis, :norders_vis, :] = weights_sens_vis

waves[:nspec_nir, norders_vis:, :] = waves_nir
fluxes[:nspec_nir, norders_vis:, :] = fluxes_nir
ivars[:nspec_nir, norders_vis:, :] = ivars_nir
gpms[:nspec_nir, norders_vis:, :] = gpms_nir
weights_sens[:nspec_nir, norders_vis:, :] = weights_sens_nir


par = par_nir['coadd1d']
wave_grid_mid, (wave_coadd, flux_coadd, ivar_coadd, gpm_coadd), order_stacks = coadd.ech_combspec(
    waves, fluxes, ivars, gpms, weights_sens, nbest=par['nbest'],
                         sn_smooth_npix=par['sn_smooth_npix'],
                         wave_method=par['wave_method'],
                         spec_samp_fact=par['spec_samp_fact'],
                         ref_percentile=par['ref_percentile'],
                         maxiter_scale=par['maxiter_scale'],
                         sigrej_scale=par['sigrej_scale'],
                         scale_method=par['scale_method'],
                         sn_min_medscale=par['sn_min_medscale'],
                         sn_min_polyscale=par['sn_min_polyscale'],
                         maxiter_reject=par['maxiter_reject'],
                         lower=par['lower'], upper=par['upper'],
                         maxrej=par['maxrej'], sn_clip=par['sn_clip'],
                         debug=False, show=True, debug_scale=False, show_exp=False)

# This is hacky do this right
coadd1d_vis.wave_grid_mid, coadd1d_vis.wave_coadd, coadd1d_vis.flux_coadd, coadd1d_vis.ivar_coadd, \
coadd1d_vis.gpm_coadd = wave_grid_mid, wave_coadd, flux_coadd, ivar_coadd, gpm_coadd
coadd1d_vis.header = header_out_vis
coadd1d_vis.save(coaddfile)


