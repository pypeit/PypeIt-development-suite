""" Run the QL tests """
import os

from pypeit.scripts import ql_keck_deimos

def afternoon():
    cmds = [os.path.join(os.environ("PYPEIT_DEV"), 'DEIMOS_QL_TST', 'raw'),
             '--calibs_only',  '--root DE.' '--det 3']
    pargs = ql_keck_deimos.parse_args(cmds)
    # Run
    ql_keck_deimos.main(pargs)

def one_slit():
    cmds = [os.path.join(os.environ("PYPEIT_DEV"), 'DEIMOS_QL_TST', 'raw'),
            '--science DE.20130409.20629.fits',  '--slit_spat 3:175']
    pargs = ql_keck_deimos.parse_args(cmds)
    # Run
    ql_keck_deimos.main(pargs)

# Command line execution
if __name__ == '__main__':
    afternoon()
