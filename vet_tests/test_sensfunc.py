import os

from IPython import embed

import numpy as np
from pathlib import Path

from astropy.io import fits
from astropy import table

from pypeit import sensfunc
from pypeit import inputfiles
from pypeit import coadd1d
from pypeit.spectrographs.util import load_spectrograph

def sensfunc_io_ir_test():
    sensfunc_io_tester('IR')

def sensfunc_io_uvis_test():
    sensfunc_io_tester('UVIS')

# TODO -- THERE IS NO TEST HERE?!
#   THIS WAS WRITTEN BY MILAN
def sensfunc_io_tester(algorithm, redux_out):

    # Remove any existing file from previous runs that were interrupted
    test_file = 'test_sens.fits'
    if os.path.isfile(test_file):
        os.remove(test_file) 

    # Random spec1d file of a standard from cooked
    spec1dfile = os.path.join(redux_out,
                             'shane_kast_blue', '600_4310_d55',
                             'shane_kast_blue_A', 'Science',
                              'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits')

    # Determine the spectrograph
    header = fits.getheader(spec1dfile)
    spectrograph = load_spectrograph(header['PYP_SPEC'])
    par = spectrograph.default_pypeit_par()
    par['sensfunc']['algorithm'] = algorithm

    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(spec1dfile, test_file, par['sensfunc'])

    assert sensobj.std_dict['name'] == 'FEIGE66', 'incorrect standard star found'

    # Assign junk values just to test that I/O works
    sensobj.sens = sensobj.empty_sensfunc_table(*sensobj.wave_cnts.T.shape)
    sensobj.wave = sensobj.wave_cnts.copy()
    sensobj.zeropoint = sensobj.counts.copy()
    sensobj.throughput = np.ones(sensobj.counts.shape, dtype=float)

    assert sensobj.sens['SENS_WAVE'].shape == sensobj.wave_cnts.T.shape, 'shape is wrong'

    # Write it
    sensobj.to_file(test_file, overwrite=True)

    # Check the extensions
    hdu = fits.open(test_file)
    ext = [h.name for h in hdu]
    assert ext == ['PRIMARY', 'SENS', 'WAVE', 'ZEROPOINT', 'THROUGHPUT'], \
                'incorrect extensions written'
    del hdu

    # Read it back in and check that the data are the same
    _sensobj = sensfunc.SensFunc.from_file(test_file)

    assert np.array_equal(_sensobj.wave, sensobj.wave), 'I/O error'
    assert np.array_equal(_sensobj.zeropoint, sensobj.zeropoint), 'I/O error'
    assert isinstance(_sensobj.sens, table.Table), 'sens table has wrong type'
    assert _sensobj.algorithm == sensobj.algorithm, \
        'algorithm mismatch after writing to and reading from disk'

    os.remove(test_file)

def test_sensfunc_from_onespec(redux_out):
    # Paths
    this_rdx = Path(redux_out).resolve() / 'keck_lris_red_mark4' / 'long_600_10000_d680'
    sci_dir = this_rdx / 'Science'

    ###### FIRST RUN COADD1D #######
    # set coadded file name
    coadd1d_fname = 'GD153_coadd1d.fits'
    if Path(coadd1d_fname).exists():
        # Remove results from previous runs
        os.unlink(coadd1d_fname)
    # read spec1d files to coadd
    spec1d_files = [str(f) for f in sorted(sci_dir.glob('spec1d*GD153*.fits'))]

    ## create the coadd1d config file
    cfg_lines = ['[rdx]']
    cfg_lines += ['  spectrograph = keck_lris_red_mark4']
    cfg_lines += ['[coadd1d]']
    cfg_lines += [f'  coaddfile = {coadd1d_fname}']
    cfg_lines += ['  wave_method = linear']
    cfg_lines += ['  flux_value = False']

    all_specfiles, all_obj = [], []
    # NOTE: For simplicity, this assumes that the spec1d files have only one object detected
    for ii in range(len(spec1d_files)):
        txtinfofile = spec1d_files[ii].replace('.fits', '.txt')
        meta_tbl = table.Table.read(txtinfofile,format='ascii.fixed_width')
        all_specfiles.append(spec1d_files[ii])
        all_obj.append(meta_tbl['name'][ii])
    data = table.Table()
    data['filename'] = all_specfiles
    data['obj_id'] = all_obj

    coadd1dFile = inputfiles.Coadd1DFile(config=cfg_lines, file_paths=[str(sci_dir)], data_table=data)
    # get par and spectrograph
    spectrograph, par, _ = coadd1dFile.get_pypeitpar()
    # coadd
    coAdd1d = coadd1d.CoAdd1D.get_instance(coadd1dFile.filenames,
                                           coadd1dFile.objids,
                                           spectrograph=spectrograph,
                                           par=par['coadd1d'], chk_version=par['rdx']['chk_version'])
    # Run
    coAdd1d.run()
    # Save to file
    coAdd1d.save(coadd1d_fname)

    # checks
    assert coAdd1d.coaddfile == coadd1d_fname, 'coadd1d file name is wrong'
    assert Path(coadd1d_fname).exists(), 'coadd1d file does not exist'

    ###### THEN RUN SENSFUNC #######
    # set sensfunc file name
    sens_file = 'GD153_sensfunc.fits'
    # loop on the algorithms
    algorithms = ['IR', 'UVIS']
    for algorithm in algorithms:
        if Path(sens_file).exists():
            # Remove results from previous runs
            os.unlink(sens_file)

        # set algorithm
        par['sensfunc']['algorithm'] = algorithm

        # Instantiate the SensFunc class for the requested algorithm
        sensobj = sensfunc.SensFunc.get_instance(coadd1d_fname, sens_file, par['sensfunc'], write_qa=False)
        sensobj.run()
        sensobj.to_file(sens_file, overwrite=True)
        # Read it back in and make some checks
        _sensobj = sensfunc.SensFunc.from_file(sens_file)
        # check algorithm
        assert _sensobj.algorithm == algorithm, 'algorithm mismatch after writing to and reading from disk'
        # check used archival standard star
        assert Path(_sensobj.std_cal).name == 'fGD153.dat', 'incorrect archival standard star name'
        # various checks on the sensfunc table
        assert _sensobj.sens is not None, 'sensfunc table is None'
        assert len(_sensobj.sens) == 1, 'sensfunc table has wrong length'
        assert _sensobj.sens['SENS_WAVE'][_sensobj.sens['SENS_WAVE']>0].size > 0, 'SENS_WAVE has no values'
        assert _sensobj.sens['SENS_ZEROPOINT'][_sensobj.sens['SENS_ZEROPOINT']>0].size > 0, \
            'SENS_ZEROPOINT has no values'
        assert _sensobj.sens['SENS_FLUXED_STD_FLAM'][_sensobj.sens['SENS_FLUXED_STD_FLAM']>0].size > 0, \
            'SENS_FLUXED_STD_FLAM has no values'
        assert _sensobj.sens['SENS_LOG10_BLAZE_FUNCTION'][_sensobj.sens['SENS_LOG10_BLAZE_FUNCTION']>0].size == 0, \
            'SENS_LOG10_BLAZE_FUNCTION should be all zero'
        assert _sensobj.wave[_sensobj.wave > 0].size > 0, 'wave has no values'
        assert _sensobj.zeropoint[_sensobj.zeropoint > 0].size > 0, 'zeropoint has no values'
        assert _sensobj.throughput[_sensobj.throughput > 0].size > 0, 'throughput has no values'
        os.unlink(sens_file)

    # remove coadded file
    os.unlink(coadd1d_fname)