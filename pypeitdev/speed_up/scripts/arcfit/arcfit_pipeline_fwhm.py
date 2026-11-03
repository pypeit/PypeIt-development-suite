"""
The fwhm detect_lines was actually run with in the wavelength calibration of a
dev-suite reduction: the slit's measured fwhm (stored in WaveFit.fwhm) when
[[wavelengths]] fwhm_fromlines is True, else the fwhm parameter
(autoid.set_fwhm).  Read from the setup's .par file next to Calibrations/.
"""
import re
from pathlib import Path


def wavelengths_par(wavecalib_file):
    """Return (fwhm, fwhm_fromlines) from the .par file of the reduction."""
    setup_dir = Path(wavecalib_file).resolve().parent.parent
    pars = sorted(setup_dir.glob('*.par'))
    if not pars:
        return 4.0, True
    txt = pars[-1].read_text()
    m = re.search(r'\[\[wavelengths\]\](.*?)(?=\n\s*\[\[)', txt, re.S)
    block = m.group(1) if m else txt
    fwhm = re.search(r'^\s*fwhm\s*=\s*([\d.]+)', block, re.M)
    froml = re.search(r'^\s*fwhm_fromlines\s*=\s*(\w+)', block, re.M)
    return (float(fwhm.group(1)) if fwhm else 4.0), (froml.group(1) == 'True' if froml else True)


def pipeline_fwhm(wavecalib_file, wavefit):
    fwhm_par, fromlines = wavelengths_par(wavecalib_file)
    if fromlines and wavefit.fwhm is not None:
        return float(wavefit.fwhm)
    return fwhm_par
