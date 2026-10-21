"""
Compare the FITS products of two PypeIt reductions, array by array.

Walks ``Calibrations/`` (and ``Science/`` when present) in two reduction
directories and, for every image extension and every binary-table column,
reports whether the data are identical, how many elements differ, the
maximum absolute difference, and the maximum ULP distance in the array's own
floating-point type.

Files are listed in calibration dependency order (Arc, Tiltimg, WaveCalib,
Edges, Slits, Flat, ..., then science products).

Usage::

    python compare_calibs.py <dir_a> <dir_b> [--all] [--files Edges Slits]

Written for PypeIt Issue #2208.
"""
import argparse
import re
from pathlib import Path

import numpy as np
from astropy.io import fits

# Calibration dependency order; anything not listed sorts after these
ORDER = ['Bias', 'Dark', 'Arc', 'Tiltimg', 'WaveCalib', 'Edges', 'Slits', 'Align',
         'Flat', 'Scattlight', 'Tilts', 'spec2d', 'spec1d']


def file_rank(path):
    name = path.name
    for i, prefix in enumerate(ORDER):
        if name.startswith(prefix):
            return (i, name)
    return (len(ORDER), name)


def ordered_ints(a):
    """Map floats to integers whose difference is the ULP distance."""
    itype = {4: np.int32, 8: np.int64}[a.dtype.itemsize]
    i = a.view(itype).astype(np.int64)
    if itype is np.int32:
        return np.where(i < 0, np.iinfo(np.int32).min - i, i)
    # int64: INT64_MIN - i cannot overflow for negative i
    return np.where(i < 0, np.iinfo(np.int64).min - i, i)


def compare_arrays(a, b):
    """
    Compare two arrays.

    Returns:
        dict: ``identical``, ``ndiff``, ``size``, ``maxabs`` and ``maxulp``
        (the last two are None when not applicable).
    """
    out = dict(identical=False, ndiff=None, size=None, maxabs=None, maxulp=None, note='')
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        out['note'] = f'shape {a.shape} vs {b.shape}'
        return out
    if a.dtype != b.dtype:
        out['note'] = f'dtype {a.dtype} vs {b.dtype}'
    out['size'] = a.size
    if a.dtype.kind == 'f' and b.dtype.kind == 'f':
        # Native byte order, same dtype
        dt = np.promote_types(a.dtype, b.dtype).newbyteorder('=')
        a = a.astype(dt)
        b = b.astype(dt)
        both_nan = np.isnan(a) & np.isnan(b)
        # Bitwise comparison, so that -0.0 vs 0.0 counts as a difference
        same = (a.view(f'i{dt.itemsize}') == b.view(f'i{dt.itemsize}')) | both_nan
        diff = ~same
        out['ndiff'] = int(np.sum(diff))
        out['identical'] = out['ndiff'] == 0
        if not out['identical']:
            fin = diff & np.isfinite(a) & np.isfinite(b)
            if np.any(fin):
                out['maxabs'] = float(np.max(np.abs(a[fin] - b[fin])))
                oa = ordered_ints(a[fin])
                ob = ordered_ints(b[fin])
                same_sign = (oa >= 0) == (ob >= 0)
                ulp = np.where(same_sign, np.abs(oa - ob), np.iinfo(np.int64).max)
                out['maxulp'] = int(np.max(ulp))
            if np.any(diff & ~fin):
                out['note'] += ' nonfinite-mismatch'
        return out
    if a.dtype.kind in 'fiub' and b.dtype.kind in 'fiub':
        diff = a != b
        out['ndiff'] = int(np.sum(diff))
        out['identical'] = out['ndiff'] == 0
        if not out['identical'] and a.dtype.kind != 'b':
            out['maxabs'] = float(np.max(np.abs(a[diff].astype(float) - b[diff].astype(float))))
        return out
    # Strings, objects, ...
    diff = np.asarray(a != b)
    out['ndiff'] = int(np.sum(diff))
    out['identical'] = out['ndiff'] == 0
    return out


def iter_arrays(hdul):
    """Yield (label, array) for every image extension and table column."""
    for i, hdu in enumerate(hdul):
        name = hdu.name or f'HDU{i}'
        if isinstance(hdu, fits.BinTableHDU):
            if hdu.data is None:
                continue
            for col in hdu.columns.names:
                yield f'{name}:{col}', np.asarray(hdu.data[col])
        elif hdu.data is not None:
            yield name, np.asarray(hdu.data)


def compare_files(fa, fb):
    """Compare two FITS files; return a list of (label, result)."""
    results = []
    with fits.open(fa) as ha, fits.open(fb) as hb:
        arrays_b = dict(iter_arrays(hb))
        seen = set()
        for label, a in iter_arrays(ha):
            seen.add(label)
            if label not in arrays_b:
                results.append((label, dict(identical=False, note='missing in B')))
                continue
            results.append((label, compare_arrays(a, arrays_b[label])))
        for label in arrays_b:
            if label not in seen:
                results.append((label, dict(identical=False, note='missing in A')))
    return results


def fmt(res):
    if res.get('ndiff') is None:
        return f"DIFFERS  {res.get('note', '')}"
    s = f"DIFFERS  {res['ndiff']:>9d}/{res['size']:<9d}"
    if res.get('maxabs') is not None:
        s += f"  max|d|={res['maxabs']:.4g}"
    if res.get('maxulp') is not None:
        s += f"  maxULP={res['maxulp']}"
    if res.get('note'):
        s += f"  ({res['note'].strip()})"
    return s


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dir_a', type=Path)
    parser.add_argument('dir_b', type=Path)
    parser.add_argument('--all', action='store_true', help='Also list identical arrays')
    parser.add_argument('--files', nargs='*', default=None,
                        help='Only compare files whose names start with these prefixes')
    args = parser.parse_args()

    for sub in ['Calibrations', 'Science']:
        da = args.dir_a / sub
        db = args.dir_b / sub
        if not da.is_dir() or not db.is_dir():
            continue
        files_a = {p.name: p for p in da.iterdir() if re.search(r'\.fits(\.gz)?$', p.name)}
        files_b = {p.name: p for p in db.iterdir() if re.search(r'\.fits(\.gz)?$', p.name)}
        names = sorted(set(files_a) | set(files_b), key=lambda n: file_rank(Path(n)))
        if args.files:
            names = [n for n in names if any(n.startswith(f) for f in args.files)]
        print(f'===== {sub}')
        for name in names:
            if name not in files_a or name not in files_b:
                print(f'{name:45s}  only in {"A" if name in files_a else "B"}')
                continue
            results = compare_files(files_a[name], files_b[name])
            ndiff = sum(not r['identical'] for _, r in results)
            status = 'identical' if ndiff == 0 else f'{ndiff} DIFFER'
            print(f'{name:45s}  {len(results):4d} arrays  {status}')
            for label, res in results:
                if not res['identical']:
                    print(f'    {label:40s} {fmt(res)}')
                elif args.all:
                    print(f'    {label:40s} identical')


if __name__ == '__main__':
    main()
