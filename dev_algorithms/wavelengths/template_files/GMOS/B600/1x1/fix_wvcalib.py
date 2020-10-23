import glob

from pypeit import wavecalib

from IPython import embed

# Grab files
wvcalib_files = glob.glob('*/*/wvcalib.fits')

# Loop
for ifile in wvcalib_files:
    wvcalib = wavecalib.WaveCalib.from_file(ifile)
    # Reset slit
    if len(wvcalib.spat_ids) > 1:
        raise ValueError("Uh oh")
    wvcalib.wv_fits[0].spat_id = wvcalib.spat_ids[0]

    # Write
    wvcalib.to_file(ifile, overwrite=True)
