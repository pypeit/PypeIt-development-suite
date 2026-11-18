
import os

from IPython import embed

import numpy as np

from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet


# This test check that `maskdesign_matching` method is properly assigning DEIMOS slit-mask design IDs
# to traced slits

# Load flats files
def flat_files(instr='keck_deimos'):
    if instr == 'keck_deimos':
        return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '1200G_Cooper', ifile)
                for ifile in ['d0115_0023.fits.gz', 'd0115_0024.fits.gz', 'd0115_0025.fits.gz']]
    elif instr == 'keck_mosfire':
        return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_mosfire', 'J_multi', ifile)
                    for ifile in ['m191015_0002.fits', 'm191015_0003.fits', 'm191015_0004.fits']]


def test_maskdef_id():
    instr_names = ['keck_deimos', 'keck_mosfire']
    for name in instr_names:
        # Load instrument
        instrument = load_spectrograph(name)
        par = instrument.default_pypeit_par()

        # working only on detector 1
        det=1

        # Built trace image
        traceImage = buildimage.buildimage_fromlist(instrument, det, par['calibrations']['traceframe'],
                                                    flat_files(instr=name))

        # load specific config parameters
        par = instrument.config_specific_par(traceImage.files[0])
        trace_par = par['calibrations']['slitedges']

        # Run edge trace
        edges = EdgeTraceSet(traceImage, instrument, trace_par, auto=True)

        slits = edges.get_slits()
        # Check that the `maskdef_id` assigned to the first and last slits is correct
        if name == 'keck_deimos':
            assert slits.maskdef_id[0] == 1098247, 'DEIMOS maskdef_id not found or wrong'
            assert slits.maskdef_id[-1] == 1098226, 'DEIMOS maskdef_id not found or wrong'
            # `maskdef_id` is the slit_id ("dSlitId") from the DEIMOS slit-mask design. If the matching is done
            # correctly these are the slit_id values that the first and last slits should have. These values are
            # coming from a reduction run with DEEP2 IDL-based pipeline.
        elif name == 'keck_mosfire':
            assert slits.maskdef_id[0] == 11, 'MOSFIRE maskdef_id not found or wrong'
            assert slits.maskdef_id[-1] == 1, 'MOSFIRE maskdef_id not found or wrong'
            # `maskdef_id` is the slit_id ("Slit_Number") from the MOSFIRE slit-mask design.
            # The expected values are determined by visual inspection.


def _trace_at(edges, pos):
    """
    Return the index of the edge trace nearest ``pos`` at the PCA reference row.

    Parameters
    ----------
    edges : :class:`~pypeit.edgetrace.EdgeTraceSet`
        Edge traces with a PCA decomposition.
    pos : :obj:`float`
        Spatial position (pixels) at the PCA reference row.

    Returns
    -------
    :obj:`int`
        Index of the nearest trace.
    """
    return int(np.argmin(np.abs(edges.edge_fit[edges.pca.reference_row, :] - pos)))


def test_add_missing_slits():
    instr_names = ['keck_deimos', 'keck_mosfire']
    for name in instr_names:
        # Load instrument
        instrument = load_spectrograph(name)
        par = instrument.default_pypeit_par()

        # working only on detector 1
        det = 1

        # Built trace image
        traceImage = buildimage.buildimage_fromlist(instrument, det,
                                                    par['calibrations']['traceframe'],
                                                    flat_files(instr=name))

        # load specific config parameters
        par = instrument.config_specific_par(traceImage.files[0])
        trace_par = par['calibrations']['slitedges']

        # Running the EdgeTraceSet steps (one-by-one)

        # Initialize EdgeTraceSet
        edges = EdgeTraceSet(traceImage, instrument, trace_par, auto=False)
        # Perform the initial edge detection and trace identification
        edges.initial_trace()
        # Initial trace can result in no edges found
        if not edges.is_empty:
            # Refine the locations of the trace using centroids of the
            # features in the Sobel-filtered image.
            edges.centroid_refine()
        # Initial trace can result in no edges found, or centroid
        # refinement could have removed all traces (via `check_traces`)
        if not edges.is_empty:
            # Fit the trace locations with a polynomial
            edges.fit_refine()
            # Use the fits to determine if there are any discontinous
            # trace centroid measurements that are actually components
            # of the same slit edge
            edges.merge_traces()
        # Check if the PCA decomposition is possible; this should catch
        # long slits
        if edges.par['auto_pca'] and edges.can_pca():
            # Use a PCA decomposition to parameterize the trace
            # functional forms
            edges.pca_refine()
            # Use the results of the PCA decomposition to rectify and
            # detect peaks/troughs in the spectrally collapsed
            # Sobel-filtered image, then use those peaks to further
            # refine the edge traces
            edges.peak_refine(rebuild_pca=True)

        # Traces to remove, identified by their position at the PCA reference
        # row.  Two traces NOT from the same slit, then a left/right pair
        # from the same slit.  The traces are located by position rather than
        # by index because the number of intermediate traces found by the
        # edge detection (and hence the index of a given edge) can change
        # with small changes in the Sobel-image peak fitting.
        removed = {'keck_deimos': [458, 1736, 496, 561],
                   'keck_mosfire': [349, 2032, 973, 1100]}
        indx = np.zeros(edges.ntrace, dtype=bool)
        for pos in removed[name]:
            i = _trace_at(edges, pos)
            assert round(edges.edge_fit[edges.pca.reference_row, i]) == pos, \
                        f'wrong {name} trace position; expected {pos}'
            indx[i] = True
        edges.remove_traces(edges.synced_selection(indx, mode='ignore'), rebuild_pca=True)

        # Run slitmask matching algorithm
        edges.maskdesign_matching()

        # Sync left and right edges
        edges.sync()

        # Check the values of the traces that have been recovered.  These
        # values are obtained by running PypeIt with the "adding missing
        # traces" functionality.  These are within a few pixels (max 4
        # pixels) from their original values. The correctness of these values
        # was also tested by visual inspection of the EdgeTraceSet image.
        recovered = {'keck_deimos': [459, 1737, 497, 561],
                     'keck_mosfire': [348, 2032, 973, 1101]}
        for pos in recovered[name]:
            i = _trace_at(edges, pos)
            assert round(edges.edge_fit[edges.pca.reference_row, i]) == pos, \
                        f'{name} trace position {pos} not recovered'


def test_overlapped_slits():
    # Load flats frames that have overlapping alignment slits.
    deimos_flats = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                 '1200G_alignslit', ifile)
                        for ifile in ['DE.20110529.11732.fits', 'DE.20110529.11813.fits',
                                      'DE.20110529.11895.fits']]

    # Load instrument
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()

    # working only on detector 2
    det = 2

    # Built trace image
    traceImage = buildimage.buildimage_fromlist(keck_deimos, det,
                                                par['calibrations']['traceframe'], deimos_flats)

    # load specific config parameters
    par = keck_deimos.config_specific_par(traceImage.files[0])
    trace_par = par['calibrations']['slitedges']

    # Run edge trace
    edges = EdgeTraceSet(traceImage, keck_deimos, trace_par, auto=True)
    slits = edges.get_slits()
    # Check that the total number of expected slits and the number of alignment slits are correct.
    assert len(slits.maskdef_designtab['MASKDEF_ID'].data) == 22, \
                'wrong number of slits for this detector'
    assert np.sum(slits.maskdef_designtab['ALIGN'].data) == 3, 'wrong number of alignment slits'
    # These number have been verified by visual inspection.


