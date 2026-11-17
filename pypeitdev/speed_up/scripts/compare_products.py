#!/usr/bin/env python
"""
Compare two directories of PypeIt science products element by element.

Written for the PR B (detector parallelism) validation, which requires that
``--ncpu 1`` and ``--ncpu N`` produce bit-identical ``spec1d``/``spec2d`` files,
and that two ``--ncpu N`` runs of the same data also agree bit for bit.

Reports, per file, how many arrays were compared and how many differ, and for
any that do, the maximum absolute and relative difference and the number of
differing elements -- so a failure says *how far off* it is, not merely that it
is off.

Usage::

    python compare_products.py DIR_A DIR_B [--label-a serial --label-b 'ncpu 4']
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def _diff_stats(a, b):
    """
    Summarize the difference between two arrays.

    Args:
        a, b (`numpy.ndarray`_): Arrays to compare.

    Returns:
        :obj:`tuple`: Number of differing elements, maximum absolute
        difference, and maximum relative difference.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    finite = np.isfinite(a) & np.isfinite(b)
    neq = a != b
    # NaNs in the same places are not a difference
    neq &= finite | (np.isnan(a) != np.isnan(b))
    ndiff = int(np.sum(neq))
    if ndiff == 0:
        return 0, 0.0, 0.0
    d = np.abs(a[finite] - b[finite])
    scale = np.maximum(np.abs(a[finite]), np.abs(b[finite]))
    maxabs = float(np.max(d)) if d.size else 0.0
    rel = d[scale > 0] / scale[scale > 0]
    maxrel = float(np.max(rel)) if rel.size else 0.0
    return ndiff, maxabs, maxrel


def compare_fits(fa, fb, verbose=True):
    """
    Compare every numeric image extension of two FITS files.

    Args:
        fa, fb (`Path`_): The two files.
        verbose (:obj:`bool`, optional): Print a line per differing extension.

    Returns:
        :obj:`tuple`: Number of extensions compared and number that differ.
    """
    ncomp = 0
    ndiff = 0
    with fits.open(fa) as ha, fits.open(fb) as hb:
        names_a = [h.name for h in ha]
        names_b = [h.name for h in hb]
        if names_a != names_b:
            print(f'  !! {fa.name}: different extensions')
            print(f'     only in A: {sorted(set(names_a) - set(names_b))}')
            print(f'     only in B: {sorted(set(names_b) - set(names_a))}')
            return 0, 1
        for ea, eb in zip(ha, hb):
            if ea.data is None:
                continue

            # Binary tables (spec1d holds one row group per extracted object,
            # with the spectra as columns) are compared column by column.  These
            # must not be skipped: for spec1d they are the entire content.
            if isinstance(ea, (fits.BinTableHDU, fits.TableHDU)):
                if len(ea.data) != len(eb.data):
                    print(f'  !! {fa.name}:{ea.name}: {len(ea.data)} rows vs {len(eb.data)}')
                    ndiff += 1
                    continue
                for col in ea.columns.names:
                    ca, cb = ea.data[col], eb.data[col]
                    if ca.dtype.kind not in 'fiub':
                        continue
                    ncomp += 1
                    if np.array_equal(ca, cb):
                        continue
                    ndiff += 1
                    n, mabs, mrel = _diff_stats(ca, cb)
                    if verbose:
                        print(f'  DIFFERS {ea.name}:{col:18s} {n:>10,}/{ca.size:,} elements  '
                              f'max|d|={mabs:.3e}  maxrel={mrel:.3e}')
                continue

            if not isinstance(ea.data, np.ndarray) or ea.data.dtype.kind not in 'fiu':
                continue
            ncomp += 1
            if np.array_equal(ea.data, eb.data):
                continue
            ndiff += 1
            n, mabs, mrel = _diff_stats(ea.data, eb.data)
            if verbose:
                print(f'  DIFFERS {ea.name:24s} {n:>10,}/{ea.data.size:,} elements  '
                      f'max|d|={mabs:.3e}  maxrel={mrel:.3e}')
    return ncomp, ndiff


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir_a', type=str)
    ap.add_argument('dir_b', type=str)
    ap.add_argument('--label-a', type=str, default='A')
    ap.add_argument('--label-b', type=str, default='B')
    cli = ap.parse_args()

    A = Path(cli.dir_a)
    B = Path(cli.dir_b)
    print(f'A = {cli.label_a}: {A}')
    print(f'B = {cli.label_b}: {B}')
    print()

    files = sorted([p for p in A.glob('spec*.fits')])
    if len(files) == 0:
        print('No spec1d/spec2d FITS files found in A.')
        return 2

    tot_comp = 0
    tot_diff = 0
    missing = 0
    for fa in files:
        fb = B / fa.name
        if not fb.exists():
            print(f'  !! {fa.name}: missing from B')
            missing += 1
            continue
        ncomp, ndiff = compare_fits(fa, fb)
        tot_comp += ncomp
        tot_diff += ndiff
        status = 'identical' if ndiff == 0 else f'{ndiff} DIFFER'
        print(f'{fa.name[:64]:66s} {ncomp:3d} arrays  {status}')

    print()
    print(f'TOTAL: {tot_comp} arrays compared across {len(files)} files, '
          f'{tot_diff} differ, {missing} files missing')
    verdict = 'BIT-IDENTICAL' if (tot_diff == 0 and missing == 0) else 'DIFFERENCES FOUND'
    print(f'VERDICT: {verdict}')
    return 0 if (tot_diff == 0 and missing == 0) else 1


if __name__ == '__main__':
    sys.exit(main())
