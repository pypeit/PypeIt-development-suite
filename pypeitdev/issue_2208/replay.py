"""
Replay captured :func:`pypeit.core.pca.pca_decomposition` calls in isolation.

Each input was saved by ``pca_patch`` (``edge_harness.py run --capture``).
For every captured call, the replay calls the real PypeIt function
``--nrep`` times, with the input vectors copied to each byte offset in
``--offsets`` within an over-allocated buffer.  It reports the number of
distinct outputs (hash of coefficients + components), the max |diff| from
the first output, and a hash of the first output that can be compared
across processes and environments.

Usage::

    python replay.py <pca_call_*.npz> [...] [--nrep 5] [--offsets 0 8 16 ...] \\
        [--fix full|seed]

Written for PypeIt Issue #2208.
"""
import argparse
import hashlib
import os
import time

import numpy as np


def at_offset(a, offset):
    """Copy ``a`` into a buffer so its data starts ``offset`` bytes past 64-byte alignment."""
    raw = np.empty(a.nbytes + 128, dtype=np.uint8)
    start = (-raw.ctypes.data) % 64 + offset
    out = raw[start:start + a.nbytes].view(a.dtype).reshape(a.shape)
    out[...] = a
    assert out.ctypes.data % 64 == offset % 64
    return out


def digest(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()[:16]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('files', nargs='+')
    parser.add_argument('--nrep', type=int, default=5)
    parser.add_argument('--offsets', type=int, nargs='+', default=[0])
    parser.add_argument('--fix', default=None, choices=['full', 'seed'])
    args = parser.parse_args()

    from pypeit import log
    log.setLevel('WARNING')
    import pca_patch
    pca_patch.install(fix=args.fix)
    from pypeit.core import pca as pcamod
    import threadpoolctl
    arch = sorted({i.get('architecture') for i in threadpoolctl.threadpool_info()
                   if i['internal_api'] == 'openblas'})

    print(f'pid={os.getpid()} fix={args.fix} openblas={arch} '
          f'threads={os.environ.get("OPENBLAS_NUM_THREADS")}')
    for fn in args.files:
        d = np.load(fn)
        vectors = d['vectors']
        npca = None if int(d['npca']) < 0 else int(d['npca'])
        pcav = None if np.isnan(d['pca_explained_var']) else float(d['pca_explained_var'])
        mean = None if d['mean'].size == 0 else d['mean']
        outputs = []
        t0 = time.perf_counter()
        for off in args.offsets:
            for _ in range(args.nrep):
                pca_patch.CALLS.clear()
                v = at_offset(vectors, off)
                m = None if mean is None else at_offset(mean, off)
                coeffs, comps, pmean, _ = pcamod.pca_decomposition(v, npca=npca,
                                                                   pca_explained_var=pcav, mean=m)
                outputs.append((off, digest(coeffs, comps, pmean), coeffs, comps,
                                pca_patch.CALLS[-1]['solvers']))
        dt = (time.perf_counter() - t0) / len(outputs)
        ref = outputs[0]
        maxd = max(max(np.max(np.abs(o[2] - ref[2])), np.max(np.abs(o[3] - ref[3])))
                   for o in outputs)
        distinct = len({o[1] for o in outputs})
        per_offset = {off: len({o[1] for o in outputs if o[0] == off}) for off in args.offsets}
        same_as_captured = ref[1] == digest(d['coeffs'], d['components'], d['pca_mean'])
        print(f'{os.path.basename(fn)} shape={vectors.shape} solvers={ref[4]} '
              f'{distinct}/{len(outputs)} distinct, distinct per offset={per_offset}, '
              f'max|d|={maxd:.2e}, first={ref[1]}, equals captured={same_as_captured}, '
              f'{1e3 * dt:.1f} ms/call')


if __name__ == '__main__':
    main()
