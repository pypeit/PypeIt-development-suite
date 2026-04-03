import numpy as np

from pypeit.core.arc import fit2darc
from pypeit import wavecalib
from pypeit import edgetrace
from IPython import embed


def load_data():
    dirc = 'PypeIt_Templates/'
    setups = [['564_l', '564_u'],
             ['580_l', '580_u']]
    alldata = [[None for dd in range(len(setups[0]))] for all in setups]
    allorders = [[None for dd in range(len(setups[0]))] for all in setups]
    for ii in range(len(setups)):
        # Load the data for each detector
        for dd in range(2):
            grat, chip = setups[ii][dd].split('_')
            # Load the edges
            edges = edgetrace.EdgeTraceSet.from_file(f'{grat}/Edges_{chip}.fits.gz')
            edge_fit = edges.edge_fit
            slitcen = np.mean(edge_fit.reshape((edge_fit.shape[0], edge_fit.shape[1] // 2, 2)), axis=2)
            nspec, norders = slitcen.shape
            wvcal_file = f"{dirc}/WaveCalib_{grat}{chip}.fits"
            wvcal = wavecalib.WaveCalib.from_file(wvcal_file)
            orders = wvcal.ech_orders
            assert orders.size==norders
            thisdata = np.zeros((nspec, norders, 3))
            for ord in range(norders):
                wave = wvcal.wv_fit2d[0].eval(np.linspace(0.0, 1.0, nspec), x2=np.full(nspec, orders[ord]))/orders[ord]
                spec = np.arange(nspec)
                spat = slitcen[:,ord]
                thisdata[:, ord, 0] = wave
                thisdata[:, ord, 1] = spec
                thisdata[:, ord, 2] = spat
            alldata[ii][dd] = thisdata.copy()
            allorders[ii][dd] = orders.copy()
    return alldata, allorders


if __name__ == '__main__':

    # Load the data
    data, orders = load_data()
