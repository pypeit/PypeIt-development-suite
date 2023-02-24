from pathlib import Path
import sys

from pypeit.scripts.sensfunc import SensFunc
from multiprocessing import Pool
import os

def run_sensfunc(options):
    try:
        pid = os.getpid()
        parser = SensFunc.get_parser()
        args = parser.parse_args(options)
        SensFunc.main(args)
    except Exception as e:
        fits_file = options[0]
        print(f"Failed to run sensfunc on file {fits_file} to {dest_file}, exception {e}")
        with open(f"failures.{pid}.txt", "a") as f:
            print(f"Failed to run sensfunc on file {fits_file} to {dest_file}, exception {e}",file=f)
        return "FAILED"

    return "OK"



if __name__ == '__main__':

    file_list = sys.argv[1]
    src_dir = Path(sys.argv[2])
    dest_dir = Path(sys.argv[3])
    sens_file = sys.argv[4]
    jobs = []

    with open(file_list, "r") as f:
        files = [Path(x.strip()) for x in f.readlines()]

    for file in files:
        spec1d_file = src_dir / file.parent / ("spec1d_" + file.name)
        dest_file = dest_dir / file.parent / ("sens_" + file.name)
        dest_file.parent.mkdir(parents=True, exist_ok=True)
        options = [f"{spec1d_file}", "-s", sens_file, "--outfile", str(dest_file)]
        jobs.append(options)

    with Pool(processes=8) as p:
        p.map(run_sensfunc, jobs)