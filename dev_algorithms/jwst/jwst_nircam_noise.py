import os
import numpy as np
import scipy
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.stats import sigma_clipped_stats

from IPython import embed

# set environment variables
os.environ['CRDS_PATH'] = '/Users/joe/crds_cache/jwst_pub'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds-pub.stsci.edu'
from matplotlib import pyplot as plt
from astropy.io import fits
from gwcs import wcstools

## JWST imports
# The calwebb_spec and spec3 pipelines
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline
# individual steps
from jwst.assign_wcs import AssignWcsStep
from jwst.background import BackgroundStep
from jwst.imprint import ImprintStep
from jwst.msaflagopen import MSAFlagOpenStep
from jwst.extract_2d import Extract2dStep
from jwst.srctype import SourceTypeStep
from jwst.wavecorr import WavecorrStep
from jwst.flatfield import FlatFieldStep
from jwst.pathloss import PathLossStep
from jwst.barshadow import BarShadowStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleSpecStep
from jwst.extract_1d import Extract1dStep
from jwst import datamodels

DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']

# PypeIt imports
from jwst_utils import compute_diff, get_cuts, jwst_nircam_proc, jwst_nircam_subimgs
from pypeit.display import display
from pypeit import specobjs
from pypeit import slittrace
from pypeit.utils import inverse, fast_running_median
from pypeit.core import findobj_skymask
from pypeit.core import skysub, coadd
from pypeit.core import procimg
from pypeit.core import flat

from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import pypeitimage
from pypeit import calibrations
from pypeit import find_objects
from pypeit import extraction
from pypeit import msgs
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit.images.pypeitimage import PypeItImage
from pypeit.scripts.show_2dspec import show_trace


DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']

# NIRCAM 2 flat field exposures
rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRCAM_WFSS/GTO/'
output_dir = '/Users/joe/jwst_redux/redux/NIRCAM_WFSS/GTO/calwebb/output/'
pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRCAM_WFSS/GTO/calwebb/pypeit/'

scifile1 = os.path.join(rawpath_level2, 'jw01243001001_02101_00001_nrcalong_rate.fits')
scifile2 = os.path.join(rawpath_level2, 'jw01243001001_02101_00002_nrcalong_rate.fits')
scifile3 = os.path.join(rawpath_level2, 'jw01243001001_02101_00003_nrcalong_rate.fits')
scifiles = [scifile1, scifile2, scifile3]


# TODO: This needs to be defined by the user
scipath = os.path.join(pypeit_output_dir, 'Science')
if not os.path.isdir(scipath):
    msgs.info('Creating directory for Science output: {0}'.format(scipath))
    os.makedirs(scipath)

# NIRCAM Mirage configuration file
condir = '/Users/joe/jwst_redux/Raw/NIRCAM_WFSS/GTO/sens/'
filt = 'F356W'
module='a'
grism = 'R'
configfile = os.path.join(condir, 'NIRCAM_' + filt + '_mod' + module.upper() + '_' + grism + '.conf')
detname ='a'
# print(configfile)


# TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
param_dict = {
    'extract_2d': {'skip': True},
    'bkg_subtract': {'skip': True},
    'flat_field': {'save_interpolated_flat': True},
    'wfss_contam': {'skip': True},
    'photom': {'skip': True},
    'extract_1d': {'skip': True},
}

basenames = []
for sci in scifiles:
    basenames.append(os.path.basename(sci).replace('_rate.fits', ''))

# Run the spec2 pipeline
runflag = False
if runflag:
    for sci in scifiles:
        step = AssignWcsStep()
        step.save_results = True
        step.output_dir = output_dir
        result = step.run(sci)

        step = FlatFieldStep()
        step.output_dir = output_dir
        # step.override_flat = new_flat_file
        step.save_results = True
        result = step.run(result)

#        spec2 = Spec2Pipeline(steps=param_dict)
#        spec2.save_results = True
#        spec2.output_dir = output_dir
#        result = spec2(sci)


# Some pypeit things
spectrograph = load_spectrograph('jwst_nircam')
par = spectrograph.default_pypeit_par()
det_container = spectrograph.get_detector_par(1)

pypeline = 'MultiSlit'
par['rdx']['redux_path'] = pypeit_output_dir
qa_dir = os.path.join(pypeit_output_dir, 'QA')
par['rdx']['qadir'] = 'QA'
png_dir = os.path.join(qa_dir, 'PNGs')
if not os.path.isdir(qa_dir):
    msgs.info('Creating directory for QA output: {0}'.format(qa_dir))
    os.makedirs(qa_dir)
if not os.path.isdir(png_dir):
    os.makedirs(png_dir)


# TODO Fix this, currently does not work if target names have - or _
filename_first = os.path.basename(scifiles[0])
filename_last = os.path.basename(scifiles[-1])
split_first = filename_first.split('_')
split_last = filename_last.split('_')
prefix_first = ''
for ii in range(len(split_first) - 3):
    prefix_first += split_first[ii] + "_"
out_filename = prefix_first + split_first[-3] + "-" + split_last[-3]

show = True
show_peaks=True
show_skysub_fit=True

flat_output_files = []
for base in basenames:
    flat_output_files.append(os.path.join(output_dir, base + '_flatfieldstep.fits'))


proc=False
if proc:
    RA, DEC = 15.05428, 28.040513

    rate_obj0, science0, sciivar0, gpm0, dq_gpm0, base_var, count_scale, finitemask, tilts, waveimg, spat_img = jwst_nircam_proc(
        flat_output_files[0], configfile, RA, DEC, kludge_err=1.0, noise_floor=0.01, saturation=65000)
    rate_obj1, science1, sciivar1, gpm1, dq_gpm1, base_var, count_scale, finitemask, tilts, waveimg, spat_img = jwst_nircam_proc(
        flat_output_files[1], configfile, RA, DEC, kludge_err=1.0, noise_floor=0.01, saturation=65000)

    # KLUDGETOWN
    science1 = science1[1:, :]
    sciivar1 = sciivar1[1:, :]
    gpm1 = gpm1[1:]

    gpm_diff = gpm0 & gpm1
    diff = science0 - science1
    sig0 = np.sqrt(inverse(sciivar0))
    sig1 = np.sqrt(inverse(sciivar1))
    var_diff = inverse(sciivar0) + inverse(sciivar1)
    sig_diff = np.sqrt(var_diff)
    ivar_diff = inverse(var_diff)
    chi = diff * np.sqrt(ivar_diff)

    xlo = 0
    xhi = 12
    ylo = 0
    yhi = science1.shape[0]

else:
    xlo = 1130
    xhi = 1800
    ylo = 800
    yhi = 865

    sci_obj0 = datamodels.open(flat_output_files[0])
    sci_obj1 = datamodels.open(flat_output_files[1])

    science0 = sci_obj0.data
    err0 = sci_obj0.err
    sci_var_poisson0 = sci_obj0.var_poisson
    sci_var_rnoise0 = sci_obj0.var_rnoise
    # Good pixel masks
    gpm_sci0 = np.logical_not(sci_obj0.dq & DO_NOT_USE)
    gpm_sci1 = np.logical_not(sci_obj1.dq & DO_NOT_USE)

    science1 = sci_obj1.data
    err1 = sci_obj1.err
    sci_var_poisson1 = sci_obj1.var_poisson
    sci_var_rnoise1 = sci_obj1.var_rnoise

    gpm_diff = gpm_sci0 & gpm_sci1
    diff = science0 - science1
    var_diff_calwebb = err0**2 + err1**2

    var_diff = sci_var_poisson0 + sci_var_rnoise0 + sci_var_poisson1 + sci_var_rnoise1

    sig_diff =np.sqrt(var_diff)
    ivar_diff = inverse(var_diff)
    chi = diff*np.sqrt(ivar_diff)


#nspat, nspec = sci_rate.shape
cuts = get_cuts(science0[ylo:yhi, xlo:xhi])
#cuts_var = get_cuts(sci_err[ylo:yhi, xlo:xhi]**2)

# yvals = ylo + np.arange(yhi - ylo)
#slit_left = np.full(nspec, ylo)
#slit_righ = np.full(nspec, yhi)
#spec_val = xlo + np.arange(xhi - xlo)
#viewer_sci, ch_sci = display.show_image(sci_rate.T, cuts=get_cuts(sci_rate), chname='raw', clear=True)
#display.show_slits(viewer_sci, ch_sci, slit_left, slit_righ, spec_vals=spec_val, pstep=1)

display.show_image(science0[ylo:yhi, xlo:xhi], cuts=cuts, chname='science', wcs_match=True)
#display.show_image(sig0[ylo:yhi, xlo:xhi]**2, cuts=cuts_var, chname='sci_var', wcs_match=True)
#display.show_image(sci_var_poisson[ylo:yhi, xlo:xhi], cuts=cuts_var, chname='sci_var_poi', wcs_match=True)
#display.show_image(sci_var_rnoise[ylo:yhi, xlo:xhi], cuts=cuts_var, chname='sci_var_rn', wcs_match=True)
display.show_image(science1[ylo:yhi, xlo:xhi], cuts=cuts, chname='bkg1', wcs_match=True)
#display.show_image(sig1[ylo:yhi, xlo:xhi], cuts=cuts, chname='bkg1', wcs_match=True)
display.show_image(diff[ylo:yhi, xlo:xhi], cuts=get_cuts(diff), chname='diff', wcs_match=True)
#display.show_image(sig_diff[ylo:yhi, xlo:xhi], cuts=get_cuts(sig_diff), chname='sigma', wcs_match=True)
display.show_image(chi[ylo:yhi, xlo:xhi], cuts=(-5.0,5.0), chname='chi', wcs_match=True)

# Grab only non-flagged pixels
chi_sub = chi[ylo:yhi, xlo:xhi]
gpm_diff_sub = gpm_diff[ylo:yhi, xlo:xhi]
chi_flat = chi_sub[gpm_diff_sub].flatten()
# Compute the mean and standard deviation of the chi distribution
mean, med, sigma = sigma_clipped_stats(chi_flat, sigma_lower=5.0,sigma_upper=5.0)
sig_range = 5.0*sigma

n_bins = 50
binsize = 2.0 * sig_range / n_bins
bins_histo = -sig_range + np.arange(n_bins) * binsize + binsize / 2.0
xvals = np.arange(-10.0, 10, 0.02)
gauss = scipy.stats.norm(loc=0.0, scale=1.0)

plt.figure(figsize=(12, 8))

plt.hist(chi_flat, bins_histo, density=True, histtype='step', align='mid', color='k', linewidth=3,
         label='Chi distribution, sigma={:5.3f}'.format(sigma))
plt.plot(xvals, gauss.pdf(xvals), 'c-', lw=3, label='sigma=1')
plt.ylabel('Residual distribution')
plt.xlabel('chi')
plt.xlim([-sig_range, sig_range])
plt.legend(fontsize=13, loc=2)
plt.show()

#rate_test = sci_rate[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]
#var_test = (sci_err[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]**2)
#var_poi_test = sci_var_poisson[ylo:yhi, xlo:xhi][gpm_sci[ylo:yhi, xlo:xhi]]

#plt.plot(rate_test, var_poi_test, color='black', marker='o', markersize=3, linestyle='None', label='Poisson')
#plt.plot(rate_test, 1/2946.0*rate_test, color='blue',linestyle='solid', label='1/4000*rate')
#plt.xlim([-0.02,0.15])
#plt.ylim([-1e-5,1.5e-4])
#plt.show()
