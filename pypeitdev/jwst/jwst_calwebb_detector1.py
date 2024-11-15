import os
# Use single cores (forcing it for numpy operations)
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import traceback
import configparser
from pypeit import msgs
from IPython import embed
import sys 
from multiprocessing import Pool
from tqdm import tqdm



def imap_unordered_bar(func, args, nproc):
    """
    Display progress bar.
    """
    p = Pool(processes=nproc)
    res_list = []
    with tqdm(total = len(args)) as pbar:
        for i, res in tqdm(enumerate(p.imap_unordered(func, args))):
            pbar.update()
            res_list.append(res)
    pbar.close()
    p.close()
    p.join()
    return res_list

def mk_stpipe_log_cfg(output_dir, log_name):
    """
    Create a configuration file with the name log_name, where
    the pipeline will write all output.
    
    Parameters
    ----------
    outpur_dir: str 
        path of the output directory
    log_name: str 
        name of the log to record screen output

    Returns
    -------
        nothing
    """
    config = configparser.ConfigParser()
    config.add_section("*")
    config.set("*", "handler", "file:" + log_name)
    config.set("*", "level", "INFO")
    pipe_log_config = os.path.join(output_dir, "pipeline-log.cfg")
    config.write(open(pipe_log_config, "w"))

def run_det1(uncal_file, output_dir):
    """
    Run the Detector1 pipeline on the given file.

    Parameters
    ----------
        uncal_file: str 
            name of uncalibrated file to run
        outpur_dir: str 
            path of the output directory
    Returns
    -------
        nothing
    """
    log_name = os.path.basename(uncal_file).replace('.fits', '')
    mk_stpipe_log_cfg(output_dir, log_name+'.log')
    from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
    pipe_success = False
    try:
        det1 = Detector1Pipeline()
        det1.call(uncal_file, output_dir=output_dir, logcfg="pipeline-log.cfg", save_results=True)
        pipe_success = True
        msgs.info('\n * Pipeline finished for file: ', uncal_file, ' \n')
    except Exception:
        msgs.warn('\n *** OH NO! The detector1 pipeline crashed! *** \n')
        pipe_crash_msg = traceback.print_exc()
    if not pipe_success:
        crashfile = open(log_name+'_pipecrash.txt', 'w')
        msgs.info('Printing file with full traceback')
        print(pipe_crash_msg, file=crashfile)

def run_calwebb_detector1(uncalfiles, output_dir, cores2use=1, overwrite=False):
    """
    Run the Detector1 pipeline on the provided files with the option to use multiprocessing.

    Parameters
    ----------
    uncalfiles: list
        list of uncalibrated files to run on. 
    output_dir: str
        path of the output directory
    cores2use: int or str 
        Number of available cores that will be used for multi-processing 
          - If int, the number of cores to use.
          - If str, one of 'quarter', 'half', 'all' to use a fraction of all available cores. 
          - The default value is 1, which results in no multi-processing.
        
        Note that these fractions refer to the total available cores and on most CPUs these include physical and 
        virtual cores. The clock time for the step is reduced almost linearly by the number of physical cores used 
        on all machines.  For example, on an Intel CPU with six real cores and six virtual cores, setting maximum_cores to 
        'half' results in a decrease of a factor of six in the clock time for the step to run.  Depending on the system, the clock 
        time can also decrease even more with maximum_cores set to 'all'. Setting the number of cores to an integer can be useful 
        when running on machines with a large number of cores where the user is limited in how many cores they can use.
    
    Returns
    -------
    nothing
    """

    # get the cores to use
    if isinstance(cores2use, int):
        system_cores = cores2use
        if system_cores > os.cpu_count():
            msgs.warn('The number of cores to use is larger than the available cores {}. '
                    'Using all available cores instead.'.format(os.cpu_count()))
            system_cores = os.cpu_count()
    elif isinstance(cores2use, str):
        if cores2use == 'quarter':
            system_cores = int(os.cpu_count()/4)
        elif cores2use == 'half':
            system_cores = int(os.cpu_count()/2)
        elif cores2use == 'all':
            system_cores = os.cpu_count()
        else:
            msgs.error('Invalid value for cores2use. Please use one of: quarter, half, all')
    
    files_to_run = []
    # Run the stage1 pipeline
    for uncal in uncalfiles:
        ratefile = os.path.join(output_dir,  os.path.basename(uncal).replace('_uncal', '_rate'))
        if os.path.isfile(ratefile) and not overwrite:
            msgs.info('Using existing rate file: {0}'.format(ratefile))
        else: 
            files_to_run.append(uncal)
    
    msgs.info('Runing the calwebb Detector1 pipeline on {} files'.format(len(files_to_run)))
    for file in files_to_run:
        print(' **** ', file)

    # the output list should be the same length as the files to run
    outptd = [output_dir for _ in range(len(files_to_run))]

    msgs.info('**** Using {} cores for multiprocessing.'.format(system_cores))
    # set the pool and run multiprocess
    args = [(file, out) for file in files_to_run for out in outptd]
    #with multiprocessing.Pool(system_cores) as pool:
    #    pool.starmap(run_det1, zip(files_to_run, outptd))
    
    output = imap_unordered_bar(run_det1, args, system_cores)

    print('\n * Finished multiprocessing! \n')
    

if __name__ == '__main__':
    sys.exit(run_calwebb_detector1([], ''))