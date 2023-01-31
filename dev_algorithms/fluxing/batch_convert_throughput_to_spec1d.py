from pathlib import Path
import sys

from astropy.table import Table, vstack
from astropy.io import fits
import numpy as np

from pypeit.spectrographs.keck_deimos import load_wmko_std_spectrum

# Copy of SpecObj med_s2n that masks bad counts
def masking_med_s2n(specobj):
    """Return median S/N of the spectrum
    Uses OPT_COUNTS if present and then BOX_COUNTS

    Returns:
        float
    """
    SN = 0.
    mask = specobj['OPT_MASK']
    if specobj['OPT_COUNTS'] is not None:
        SN = np.median(specobj['OPT_COUNTS'][mask] * np.sqrt(specobj['OPT_COUNTS_IVAR'][mask]))

    return SN


root_dir = Path(sys.argv[1])
dest_dir = Path(sys.argv[2])
pad = True if (len(sys.argv) >= 4 and sys.argv[3].lower() == "pad") else False
combined_table = None

for fits_file in [x.relative_to(root_dir) for x in root_dir.rglob("*.fits")]:
    dest_file = dest_dir / fits_file.parent / f'spec1d_{fits_file.name}'
    try:

        dest_file.parent.mkdir(parents=True, exist_ok=True)
        sobjs = load_wmko_std_spectrum(str(root_dir / fits_file), str(dest_file), pad)

        hdul = fits.open(str(root_dir / fits_file))
        t = Table(hdul[1].data)
        t2 = Table(hdul[2].data)        
        t.add_column([masking_med_s2n(sobjs[1])], index=0, name = "Det 7 s2n")
        t.add_column([masking_med_s2n(sobjs[0])], index=0, name = "Det 3 s2n")
        t.add_column([len(t2['COUNTS'])], index=0, name = "# counts")
        t.add_column([len(t2['COUNTS'][t2['COUNTS']<0])], index=0, name = "# counts < 0")
        t.add_column([fits_file], index=0, name = "Filename")
        if combined_table is None:
            combined_table = t
        else:
            combined_table = vstack([combined_table, t])

    except Exception as e:
        print(f"Failed to make spec1d file {dest_file} from {fits_file}, exception {e}")

combined_table.write("./throughput_index.csv", format="csv", overwrite=True)
