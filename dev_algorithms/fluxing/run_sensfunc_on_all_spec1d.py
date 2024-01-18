from pathlib import Path
import sys

from pypeit.sensfunc import SensFunc
from multiprocessing import Pool
import os
from pypeit import io
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit import inputfiles


def get_primary_hdr(spectrograph, hdul):
    
    # Construct a primary FITS header that includes the spectrograph's
    #   config keys for inclusion in the output sensfunc file
    primary_hdr = io.initialize_header()
    add_keys = (
        ['PYP_SPEC', 'DATE-OBS', 'TELESCOP', 'INSTRUME', 'DETECTOR']
        + spectrograph.configuration_keys() + spectrograph.raw_header_cards()
    )
    for key in add_keys:
        if key.upper() in hdul[0].header.keys():
            primary_hdr[key.upper()] = hdul[0].header[key.upper()]

    return primary_hdr


def run_sensfunc(options):
    spec1d_file, config_lines, output_file = options
    try:
        pid = os.getpid()
        with io.fits_open(str(spec1d_file)) as hdul:
            spectrograph = load_spectrograph(hdul[0].header['PYP_SPEC'])
            primary_hdr = get_primary_hdr(spectrograph, hdul)
            config_specific_par = spectrograph.config_specific_par(hdul)

        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=config_specific_par.to_config(),
                                                 merge_with=(config_lines,))

        # SensFunc.run() will try to flux using the new sensfunc
        # for QA purposes. Because we didn't want to extrapolate in the
        # genreated sensfuncs we need to add this to prevent failures in this fluxing
        par['fluxcalib']['extrap_sens'] = True
        par_output = output_file.parent / (output_file.stem + ".par")
        par['sensfunc'].to_config(str(par_output), section_name='sensfunc', include_descr=False)      

        sensobj = SensFunc.get_instance(str(spec1d_file), str(output_file), par['sensfunc'],par_fluxcalib=par['fluxcalib'])

        # Generate the sensfunc
        sensobj.run()
        # Write it out to a file, including the new primary FITS header
        sensobj.to_file(str(output_file), primary_hdr=primary_hdr, overwrite=True)
        

    except Exception as e:
        print(f"Failed to run sensfunc on file {spec1d_file} to {output_file}, exception {e}")
        with open(f"failures.{pid}.txt", "a") as f:
            print(f"Failed to run sensfunc on file {spec1d_file} to {output_file}, exception {e}",file=f)
        return "FAILED"

    return "OK"



if __name__ == '__main__':

    src_dir = Path(sys.argv[1])
    src_pattern = sys.argv[2]
    dest_dir = Path(sys.argv[3])
    sens_file = sys.argv[4]
    jobs = []

    config_lines = inputfiles.SensFile.from_file(sens_file).cfg_lines

    for spec1d_file in src_dir.rglob(src_pattern):
        dest_file = dest_dir / spec1d_file.name.replace("spec1d", "sens", 1)
        dest_file.parent.mkdir(parents=True, exist_ok=True)
        options = [spec1d_file, config_lines, dest_file]
        jobs.append(options)

    with Pool(processes=8) as p:
        p.map(run_sensfunc, jobs)