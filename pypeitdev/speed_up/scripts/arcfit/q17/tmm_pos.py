import os
import numpy as np
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet

def flat_files(instr):
    if instr == 'keck_deimos':
        return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '1200G_Cooper', f)
                for f in ['d0115_0023.fits.gz', 'd0115_0024.fits.gz', 'd0115_0025.fits.gz']]
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_mosfire', 'J_multi', f)
            for f in ['m191015_0002.fits', 'm191015_0003.fits', 'm191015_0004.fits']]

def _trace_at(edges, pos, side):
    """Index of the left (side=-1) or right (side=1) trace nearest ``pos``."""
    x = edges.edge_fit[edges.pca.reference_row, :]
    cand = np.arange(x.size)
    return cand[np.argmin(np.abs(x[cand] - pos))]

def test_add_missing_slits():
    removed = {'keck_deimos': [(458, -1), (1736, 1), (496, -1), (561, 1)],
               'keck_mosfire': [(349, -1), (2032, 1), (973, -1), (1100, 1)]}
    recovered = {'keck_deimos': [(459, -1), (1737, 1), (497, -1), (561, 1)],
                 'keck_mosfire': [(348, -1), (2032, 1), (973, -1), (1101, 1)]}
    for name in ['keck_deimos', 'keck_mosfire']:
        instrument = load_spectrograph(name)
        par = instrument.default_pypeit_par()
        traceImage = buildimage.buildimage_fromlist(instrument, 1, par['calibrations']['traceframe'],
                                                    flat_files(name))
        par = instrument.config_specific_par(traceImage.files[0])
        edges = EdgeTraceSet(traceImage, instrument, par['calibrations']['slitedges'], auto=False)
        edges.initial_trace(); edges.centroid_refine(); edges.fit_refine(); edges.merge_traces()
        edges.pca_refine(); edges.peak_refine(rebuild_pca=True)
        indx = np.zeros(edges.ntrace, dtype=bool)
        for pos, side in removed[name]:
            i = _trace_at(edges, pos, side)
            assert round(edges.edge_fit[edges.pca.reference_row, i]) == pos, f'{name} trace {pos}'
            indx[i] = True
        edges.remove_traces(edges.synced_selection(indx, mode='ignore'), rebuild_pca=True)
        edges.maskdesign_matching()
        edges.sync()
        for pos, side in recovered[name]:
            i = _trace_at(edges, pos, side)
            got = edges.edge_fit[edges.pca.reference_row, i]
            assert round(got) == pos, f'{name} trace {pos} not recovered ({got:.1f})'
        print(name, 'OK')

test_add_missing_slits()
