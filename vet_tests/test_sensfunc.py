import os

from IPython import embed

import numpy as np
from pathlib import Path

from astropy.io import fits
from astropy import table

from pypeit import sensfunc
from pypeit import inputfiles
from pypeit.scripts import coadd_1dspec
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

    coadd1d_fname = 'GD153_coadd1d.fits'

    rdx = Path(redux_out).resolve()
    this_rdx = Path(redux_out).resolve() / 'keck_lris_red_mark4' / 'long_600_10000_d680'
    sci_dir = this_rdx / 'Science'

    # spec1d files to coadd
    spec1d_files = [str(f) for f in sorted(sci_dir.glob('spec1d*GD153*.fits'))]

    ## coadd1d pypeit file
    cfg_file = 'GD153_coadd1d.coadd1d'
    if Path(cfg_file).exists():
        # Remove results from previous run
        os.unlink(cfg_file)
    # make it
    cfg_lines = ['[coadd1d]']
    cfg_lines += [f'  coaddfile = {coadd1d_fname}']
    cfg_lines += ['  wave_method = linear']
    cfg_lines += ['  flux_value = False']
    #TODO: TO REMOVE
    cfg_lines += ['[rdx]']
    cfg_lines += ['  chk_version = False']

    all_specfiles, all_obj = [], []
    for ii in range(len(spec1d_files)):
        txtinfofile = spec1d_files[ii].replace('.fits', '.txt')
        meta_tbl = table.Table.read(txtinfofile,format='ascii.fixed_width')
        all_specfiles.append(spec1d_files[ii])
        all_obj.append(meta_tbl['name'][ii])
    data = table.Table()
    data['filename'] = all_specfiles
    data['obj_id'] = all_obj

    coadd1dFile = inputfiles.Coadd1DFile(config=cfg_lines, file_paths=[str(sci_dir)], data_table=data)
    # Write
    coadd1dFile.write(cfg_file)

    parser = coadd_1dspec.CoAdd1DSpec.get_parser()
    args = parser.parse_args([cfg_file])
    coadd_1dspec.CoAdd1DSpec.main(args)

    # create sensfunc
    sens_file = 'GD153_sensfunc.fits'
    if Path(sens_file).exists():
        # Remove results from previous run
        os.unlink(sens_file)
    spectrograph = load_spectrograph('keck_lris_red_mark4')
    par = spectrograph.default_pypeit_par()

    # First with IR algorithm
    par['sensfunc']['algorithm'] = 'IR'

    # Instantiate the relevant class for the requested algorithm
    sensobj = sensfunc.SensFunc.get_instance(coadd1d_fname, sens_file, par['sensfunc'])
    sensobj.run()
    sensobj.to_file(sens_file, overwrite=True)

    # Then with UVIS algorithm
    os.unlink(sens_file)
    par['sensfunc']['algorithm'] = 'UVIS'
    sensobj = sensfunc.SensFunc.get_instance(coadd1d_fname, sens_file, par['sensfunc'])
    sensobj.run()
    sensobj.to_file(sens_file, overwrite=True)

    # MAKE some asserts

    # unlink everything
    os.unlink(coadd1d_fname)
    os.unlink(cfg_file)
    os.unlink('coadd1d.par')
    os.unlink(sens_file)