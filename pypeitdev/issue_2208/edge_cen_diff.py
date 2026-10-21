"""
Locate the EDGE_CEN differences between two Edges calibration files.

For each trace (column) with any difference, reports the trace ID and side,
the slit-mask ID, any insertion flags, the number of differing rows, the row
ranges, and the maximum ULP distance.  Also says whether the corresponding
EDGE_FIT column differs.

Usage::

    python edge_cen_diff.py <Edges_a.fits.gz> <Edges_b.fits.gz>

Written for PypeIt Issue #2208.
"""
import argparse

import numpy as np
from astropy.io import fits

from pypeit.edgetrace import EdgeTraceBitMask

from compare_calibs import compare_arrays


def row_ranges(rows):
    """Collapse sorted row indices into 'a-b' ranges."""
    if rows.size == 0:
        return ''
    breaks = np.where(np.diff(rows) > 1)[0]
    starts = np.r_[rows[0], rows[breaks + 1]]
    ends = np.r_[rows[breaks], rows[-1]]
    return ','.join(f'{s}' if s == e else f'{s}-{e}' for s, e in zip(starts, ends))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('file_a')
    parser.add_argument('file_b')
    args = parser.parse_args()

    bm = EdgeTraceBitMask()
    with fits.open(args.file_a) as ha, fits.open(args.file_b) as hb:
        # Arrays are (nspec, ntrace)
        cen_a = ha['EDGE_CEN'].data.astype(np.float32)
        cen_b = hb['EDGE_CEN'].data.astype(np.float32)
        fit_a = ha['EDGE_FIT'].data
        fit_b = hb['EDGE_FIT'].data
        msk = ha['EDGE_MSK'].data
        traceid = ha['TRACEID'].data
        maskdef = ha['MASKDEF_ID'].data if 'MASKDEF_ID' in ha else None
        for ext in ['SOBELSIG', 'TRACE_IMAGE', 'TRACEID', 'MASKDEF_ID', 'EDGE_MSK']:
            if ext in ha and ext in hb:
                same = compare_arrays(ha[ext].data, hb[ext].data)['identical']
                print(f'{ext:12s} {"identical" if same else "DIFFERS"}')

    if cen_a.shape != cen_b.shape:
        print(f'EDGE_CEN shapes differ: {cen_a.shape} vs {cen_b.shape}')
        return

    nspec, ntrace = cen_a.shape
    diff = cen_a.view(np.int32) != cen_b.view(np.int32)
    print(f'\nEDGE_CEN: {diff.sum()} of {diff.size} elements differ, '
          f'in {np.any(diff, axis=0).sum()} of {ntrace} traces')
    res = compare_arrays(cen_a, cen_b)
    print(f'max|d|={res["maxabs"]}  maxULP={res["maxulp"]}')

    print(f'\n{"col":>3s} {"tid":>4s} {"side":5s} {"mdef":>5s} {"ndiff":>6s} {"fitdiff":>7s} '
          f'{"maxULP":>6s}  flags / rows')
    for j in range(ntrace):
        rows = np.where(diff[:, j])[0]
        fitdiff = int(np.sum(fit_a[:, j] != fit_b[:, j]))
        if rows.size == 0 and fitdiff == 0:
            continue
        flags = set()
        for f in bm.insert_flags:
            if np.all(bm.flagged(msk[:, j], flag=f)):
                flags.add(f)
        ulp = 0
        if rows.size:
            ulp = int(np.max(np.abs(cen_a[rows, j].view(np.int32).astype(np.int64)
                                    - cen_b[rows, j].view(np.int32).astype(np.int64))))
        side = 'left' if traceid[j] < 0 else 'right'
        md = '' if maskdef is None else str(maskdef[j])
        print(f'{j:3d} {traceid[j]:4d} {side:5s} {md:>5s} {rows.size:6d} {fitdiff:7d} {ulp:6d}  '
              f'{",".join(sorted(flags)) or "-"}  rows={row_ranges(rows)}')


if __name__ == '__main__':
    main()
