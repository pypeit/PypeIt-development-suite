""" Run the QL tests """
import os

from pypeit.scripts import ql_keck_deimos

from IPython import embed

raw_path = os.path.join(os.getenv("PYPEIT_DEV"), 'DEIMOS_QL_TST', 'raw')
redux_path = os.path.join(os.getenv("PYPEIT_DEV"), 'DEIMOS_QL_TST')


def afternoon():
    QLKD = ql_keck_deimos.QLKeckDEIMOS()
    cmds = [raw_path, '--calibs_only',
            '--root=DE.', '-d=3', '--redux_path={}'.format(redux_path)]
    pargs = QLKD.parse_args(cmds)
    # Run
    QLKD.main(pargs)

def one_slit():
    QLKD = ql_keck_deimos.QLKeckDEIMOS()
    cmds = [raw_path, '--science=DE.20130409.20629.fits',  '--slit_spat=3:763',
            '--redux_path={}'.format(redux_path)]
    pargs = QLKD.parse_args(cmds)
    # Run
    QLKD.main(pargs)

# Command line execution
if __name__ == '__main__':
    #afternoon()
    one_slit()
