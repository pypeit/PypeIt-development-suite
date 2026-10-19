"""
Re-run PypeIt slit-edge tracing from saved inputs, hashing the trace state
after every EdgeTraceSet method call.

The TraceImage is read from an Edges calibration file written by
``run_pypeit``, and the spectrograph and parameters are rebuilt from the
``.pypeit`` file exactly as ``run_pypeit`` does
(``PypeItFile.get_pypeitpar``).  Every method defined by
:class:`~pypeit.edgetrace.EdgeTraceSet` is wrapped; after each call returns,
the *in-memory* (float64) trace arrays and PCA models are hashed (SHA-256),
along with the call depth and the parent call.

Each repetition writes:

    - ``<out>/<tag>_rep<N>.json``: the ordered list of calls, with hashes
    - ``<out>/<tag>_rep<N>_steps.npz``: the full state after each top-level
      step of ``auto_trace`` (for measuring the size of differences)
    - ``<out>/<tag>_rep<N>_final.npz``: the final state

Usage::

    python edge_harness.py run --pypeit <file.pypeit> --edges <Edges_*.fits.gz> \\
        --out <dir> --tag <tag> [--nrep 5]
    python edge_harness.py compare <json_a> <json_b> [...]
    python edge_harness.py check --final <..._final.npz> --edges <Edges_*.fits.gz> [...]

Written for PypeIt Issue #2208.
"""
import argparse
import functools
import hashlib
import inspect
import json
import os
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

# Arrays of the EdgeTraceSet state that are hashed after every call
STATE_ARRAYS = ['edge_img', 'edge_cen', 'edge_err', 'edge_msk', 'edge_fit', 'traceid',
                'maskdef_id', 'sobelsig']
PCA_ATTRS = ['pca', 'left_pca', 'right_pca']
# Methods that only read or display; not wrapped
SKIP = {'show', 'to_file', 'qa_plot', 'get_slits', 'from_hdu', 'from_file', 'to_hdu',
        '_bundle', '_parse', '_init_internals', '__init__'}


def _hash(a):
    if a is None:
        return None
    a = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(str(a.dtype).encode())
    h.update(str(a.shape).encode())
    h.update(a.tobytes())
    return h.hexdigest()[:16]


def pca_arrays(pca, prefix):
    """Flatten a TracePCA into a dict of named arrays."""
    out = {}
    if pca is None:
        return out
    for key in ['reference_row', 'trace_coo', 'pca_coeffs', 'pca_components', 'pca_mean']:
        val = getattr(pca, key, None)
        if val is not None:
            out[f'{prefix}:{key}'] = np.asarray(val)
    models = getattr(pca, 'pca_coeffs_model', None)
    if models is not None:
        for i, m in enumerate(models):
            for key in ['xval', 'yval', 'fitc', 'gpm']:
                val = getattr(m, key, None)
                if val is not None:
                    out[f'{prefix}:model{i}:{key}'] = np.asarray(val)
    return out


def state_arrays(ets):
    """All hashed arrays of an EdgeTraceSet."""
    out = {}
    for key in STATE_ARRAYS:
        val = getattr(ets, key, None)
        if val is not None:
            out[key] = np.asarray(val)
    for key in PCA_ATTRS:
        out.update(pca_arrays(getattr(ets, key, None), key))
    return out


def alignment(ets):
    """Byte offset modulo 64 of the main arrays' data pointers."""
    out = {}
    for key in ['edge_cen', 'edge_fit', 'edge_err', 'sobelsig']:
        val = getattr(ets, key, None)
        if isinstance(val, np.ndarray):
            out[key] = int(val.ctypes.data % 64)
    return out


class Recorder:
    """Wraps EdgeTraceSet methods and records the state after each call."""

    def __init__(self):
        self.calls = []
        self.stack = []
        self.steps = {}

    def reset(self):
        self.calls = []
        self.stack = []
        self.steps = {}

    def wrap(self, cls):
        for name, member in list(cls.__dict__.items()):
            if name in SKIP or name.startswith('__') or not inspect.isfunction(member):
                continue
            setattr(cls, name, self._wrapper(name, member))

    def _wrapper(self, name, func):
        rec = self

        @functools.wraps(func)
        def wrapped(obj, *args, **kwargs):
            parent = rec.stack[-1] if rec.stack else None
            rec.stack.append(name)
            try:
                return func(obj, *args, **kwargs)
            finally:
                rec.stack.pop()
                arrays = state_arrays(obj)
                entry = dict(i=len(rec.calls), name=name, depth=len(rec.stack), parent=parent,
                             hashes={k: _hash(v) for k, v in arrays.items()},
                             align=alignment(obj))
                rec.calls.append(entry)
                # Save the full state after each top-level step of auto_trace
                if parent == 'auto_trace':
                    rec.steps[f'{entry["i"]:04d}_{name}'] = arrays
        return wrapped


def run(args):
    for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
        if os.environ.get(v) != '1':
            raise SystemExit(f'{v} must be exported as 1 before starting Python')

    from pypeit import log
    from pypeit.inputfiles import PypeItFile
    from pypeit import edgetrace
    from pypeit.images.buildimage import TraceImage

    log.setLevel('WARNING')

    pf = PypeItFile.from_file(args.pypeit)
    spectrograph, par, _ = pf.get_pypeitpar()
    epar = par['calibrations']['slitedges']

    # The raw trace frame(s), as passed to buildimage_fromlist by Calibrations
    ftype = pf.data['frametype']
    trace_files = [f for f, t in zip(pf.filenames, ftype) if 'trace' in t.split(',')]

    def load_traceimg():
        # A fresh TraceImage for every repetition; copy.deepcopy cannot be used
        # on develop because of the BitMaskArray __getattr__ recursion (fixed
        # in #2209)
        with fits.open(args.edges) as hdu:
            traceimg = TraceImage.from_hdu(hdu, chk_version=False)
        traceimg.files = trace_files
        return traceimg

    rec = Recorder()
    rec.wrap(edgetrace.EdgeTraceSet)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for rep in range(args.nrep):
        rec.reset()
        traceimg = load_traceimg()
        t0 = time.perf_counter()
        ets = edgetrace.EdgeTraceSet(traceimg, spectrograph, epar, qa_path=None, auto=True)
        dt = time.perf_counter() - t0
        stem = out / f'{args.tag}_rep{rep}'
        with open(f'{stem}.json', 'w') as f:
            json.dump(dict(tag=args.tag, rep=rep, pid=os.getpid(), wall=dt,
                           success=bool(ets.success), calls=rec.calls), f)
        np.savez_compressed(f'{stem}_steps.npz',
                            **{f'{s}|{k}': v for s, d in rec.steps.items() for k, v in d.items()})
        np.savez_compressed(f'{stem}_final.npz', **state_arrays(ets))
        print(f'{args.tag} rep{rep} pid={os.getpid()} wall={dt:.1f}s success={ets.success} '
              f'ncalls={len(rec.calls)} final edge_cen={_hash(ets.edge_cen)}', flush=True)
        del ets


def compare(args):
    runs = []
    for fn in args.json:
        with open(fn) as f:
            runs.append(json.load(f))
    ref = runs[0]
    print(f'reference: {args.json[0]}')
    for fn, r in zip(args.json[1:], runs[1:]):
        print(f'--- vs {fn}')
        n = min(len(ref['calls']), len(r['calls']))
        if len(ref['calls']) != len(r['calls']):
            print(f'  call counts differ: {len(ref["calls"])} vs {len(r["calls"])}')
        first = None
        for i in range(n):
            a, b = ref['calls'][i], r['calls'][i]
            if a['name'] != b['name']:
                print(f'  call sequence diverges at #{i}: {a["name"]} vs {b["name"]}')
                first = i
                break
            keys = sorted(set(a['hashes']) | set(b['hashes']))
            bad = [k for k in keys if a['hashes'].get(k) != b['hashes'].get(k)]
            if bad:
                first = i
                print(f'  first differing call #{i}: {a["name"]} (depth {a["depth"]}, '
                      f'parent {a["parent"]}); differing: {", ".join(bad)}')
                break
        if first is None:
            print('  all calls identical')
            continue
        # Context: the top-level auto_trace steps and their status
        print('  top-level steps:')
        for i in range(n):
            a, b = ref['calls'][i], r['calls'][i]
            if a['parent'] != 'auto_trace':
                continue
            bad = [k for k in a['hashes'] if a['hashes'].get(k) != b['hashes'].get(k)]
            print(f'    #{i:4d} {a["name"]:22s} {"identical" if not bad else "DIFFERS: " + ", ".join(bad)}')


def check(args):
    """Compare a harness final state with the EDGE_CEN in Edges files."""
    final = np.load(args.final)
    cen = final['edge_cen'].astype(np.float32)
    for fn in args.edges:
        with fits.open(fn) as h:
            ref = h['EDGE_CEN'].data.astype(np.float32)
            sob = h['SOBELSIG'].data
        same_shape = ref.shape == cen.shape
        ndiff = int(np.sum(ref.view(np.int32) != cen.view(np.int32))) if same_shape else None
        sob_same = np.array_equal(sob.astype(np.float32).view(np.int32),
                                  final['sobelsig'].astype(np.float32).view(np.int32))
        print(f'{fn}: EDGE_CEN shape {"ok" if same_shape else "MISMATCH"}, '
              f'{ndiff} float32 elements differ; SOBELSIG {"identical" if sob_same else "DIFFERS"}')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('run')
    p.add_argument('--pypeit', required=True)
    p.add_argument('--edges', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--tag', required=True)
    p.add_argument('--nrep', type=int, default=1)
    p.set_defaults(func=run)
    p = sub.add_parser('compare')
    p.add_argument('json', nargs='+')
    p.set_defaults(func=compare)
    p = sub.add_parser('check')
    p.add_argument('--final', required=True)
    p.add_argument('--edges', nargs='+', required=True)
    p.set_defaults(func=check)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
