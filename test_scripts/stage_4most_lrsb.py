"""Stage the 4MOST LRS-B red commissioning fixture for ESO OB 5148227."""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from urllib.request import urlopen
from zipfile import ZipFile


INSTRUMENT = 'vista_4most_lrsb_red'
SETUP = 'OB5148227'
SCIENCE = 'QMOST.2026-11-27T00:47:46.453.fits.fz'
TRACE_ENV = 'PYPEIT_4MOST_TRACE'

ESO_ARCHIVE = 'https://dataportal.eso.org/dataPortal/file/'
RAW_SHA256 = {
    SCIENCE: 'c106a03764b198db3954d60a9c5adbc8bcb9bbc21a5871a14e2eb3e59621ac6b',
    'QMOST.2026-11-27T08:43:26.164.fits.fz':
        '68088a186dfbd4a0a66db1444aeae7d2aa1ae463b539d56aa6c603a23a6f6056',
    'QMOST.2026-11-27T08:45:19.991.fits.fz':
        '3b042e355f42c227f9eccde15e445e27810eb5449eb336c07aa61849a43d3a6a',
    'QMOST.2026-11-27T08:47:14.979.fits.fz':
        '7a25119497a74780b770c7e76b381f5642ee96f6a6f8444cf6ccf9cc8e2dc591',
    'QMOST.2026-11-27T09:02:29.940.fits.fz':
        '0f697295329bea64b921ca7bba50c61bd5a424ad7df61d8ac4c6ddcc3e98b84f',
    'QMOST.2026-11-27T09:05:01.781.fits.fz':
        'd739dcd5c62157763e47809bf4610afa79d6cd1105f249b55c16a14fb1f74f6a',
    'QMOST.2026-11-27T10:00:05.525.fits.fz':
        '2597111dae832ada19ae915eafe4df29ff39d1b46e43785ccdf343a57c154317',
}
DAY_FLAT = 'QMOST.2026-11-27T10:00:05.525.fits.fz'

SEED_URL = ('https://raw.githubusercontent.com/mcoughlin/PypeIt/'
            '488c7d677775ea7070ad3e691ab1767dcd118f17/reproduction.zip')
SEED_SHA256 = 'ebe1c469571c7704b690937b7d366c0967247760d86afdb7efc22a9e13b267bf'
SLITS = 'Slits_A_0_DET01.fits.gz'
CALIBRATIONS = (SLITS, 'WaveCalib_A_0_DET01.fits', 'Tilts_A_0_DET01.fits',
                'FiberFlat_A_0_DET01.fits')
TRACE = 'qmost_trace.npz'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def stage_file(source, destination, expected_hash):
    if destination.is_file() and sha256(destination) == expected_hash:
        return
    print(f'Staging {destination.name}', flush=True)
    with NamedTemporaryFile(dir=destination.parent, delete=False) as temporary:
        pending = Path(temporary.name)
        try:
            with (source.open('rb') if isinstance(source, Path)
                  else urlopen(source, timeout=120)) as stream:
                shutil.copyfileobj(stream, temporary)
            temporary.close()
            if sha256(pending) != expected_hash:
                raise ValueError(f'SHA-256 mismatch for {source}')
            pending.replace(destination)
        finally:
            pending.unlink(missing_ok=True)


def build_lamp_flat(raw_flat, calib_dir):
    import numpy as np
    from pypeit import flatfield, slittrace
    from pypeit.images import buildimage
    from pypeit.spectrographs.util import load_spectrograph

    spectrograph = load_spectrograph(INSTRUMENT)
    par = spectrograph.default_pypeit_par()['calibrations']['traceframe']
    par['process']['use_biasimage'] = False
    par['process']['mask_cr'] = False
    image = buildimage.buildimage_fromlist(spectrograph, 1, par, [str(raw_flat)])
    slits = slittrace.SlitTraceSet.from_file(calib_dir / SLITS)
    flat = flatfield.FlatImages(pixelflat_raw=image.image.astype(np.float32),
                                spat_id=slits.spat_id, PYP_SPEC=INSTRUMENT)
    flat.set_paths(calib_dir, 'A', ['0'], 'DET01')
    flat.to_file(overwrite=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outputdir', type=Path, default=Path('REDUX_OUT'),
                        help='Same output root as pypeit_test -o (default: REDUX_OUT).')
    parser.add_argument('--raw-source', type=Path, help='Use a local cache of the ESO FITS files.')
    parser.add_argument('--seed-archive', type=Path, help='Use a local copy of reproduction.zip.')
    args = parser.parse_args()
    if not os.environ.get('PYPEIT_DEV'):
        parser.error('Set PYPEIT_DEV to the Development Suite checkout first.')

    raw_dir = Path(os.environ['PYPEIT_DEV']) / 'RAW_DATA' / INSTRUMENT / SETUP
    raw_dir.mkdir(parents=True, exist_ok=True)
    for name, expected_hash in RAW_SHA256.items():
        source = (ESO_ARCHIVE + name.removesuffix('.fits.fz') if args.raw_source is None
                  else args.raw_source / name)
        stage_file(source, raw_dir / name, expected_hash)

    seed_dir = Path(os.environ['PYPEIT_DEV']) / 'CALIBS' / INSTRUMENT / SETUP
    seed_dir.mkdir(parents=True, exist_ok=True)
    archive = seed_dir / 'reproduction.zip'
    stage_file(args.seed_archive or SEED_URL, archive, SEED_SHA256)

    redux_dir = args.outputdir.resolve() / INSTRUMENT / SETUP
    calib_dir = redux_dir / 'Calibrations'
    calib_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive) as seeds:
        for name in (*CALIBRATIONS, TRACE):
            destination = (redux_dir if name == TRACE else calib_dir) / name
            with seeds.open(f'portable_seed/{name}') as stream, destination.open('wb') as output:
                shutil.copyfileobj(stream, output)
    build_lamp_flat(raw_dir / DAY_FLAT, calib_dir)
    print(f'{TRACE_ENV}={redux_dir / TRACE}')


if __name__ == '__main__':
    main()
