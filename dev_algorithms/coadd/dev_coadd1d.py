import os
import coadd1d

datapath = os.path.join(os.getenv('HOME'),'Dropbox/PypeIt_Redux/GMOS/R400_Flux/')
fnames = [datapath+'spec1d_flux_S20180903S0136-J0252-0503_GMOS-S_1864May27T160716.387.fits',\
          datapath+'spec1d_flux_S20180903S0137-J0252-0503_GMOS-S_1864May27T160719.968.fits',\
          datapath+'spec1d_flux_S20180903S0138-J0252-0503_GMOS-S_1864May27T160723.353.fits',\
          datapath+'spec1d_flux_S20180903S0141-J0252-0503_GMOS-S_1864May27T160727.033.fits',\
          datapath+'spec1d_flux_S20180903S0142-J0252-0503_GMOS-S_1864May27T160730.419.fits',\
          datapath+'spec1d_flux_S20181015S0140-J0252-0503_GMOS-S_1864May27T185252.770.fits']
gdobj = ['SPAT1073-SLIT0001-DET03','SPAT1167-SLIT0001-DET03','SPAT1071-SLIT0001-DET03','SPAT1072-SLIT0001-DET03',\
         'SPAT1166-SLIT0001-DET03','SPAT1073-SLIT0001-DET03']

# parameters for load_1dspec_to_array
ex_value = 'OPT'
flux_value = True

# Reading data
waves,fluxes,ivars,masks = coadd1d.load_1dspec_to_array(fnames,gdobj=gdobj,order=None,ex_value=ex_value,flux_value=flux_value)

# Coadding
wave_stack, flux_stack, ivar_stack, mask_stack, scale_array = \
    coadd1d.long_comb(waves, fluxes, ivars, masks,wave_method='pixel', scale_method='median', maxiter_reject = 5, \
                      qafile='J0252_gmos', outfile='J0252_gmos.fits', verbose=False, debug=True)
