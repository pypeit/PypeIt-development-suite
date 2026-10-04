.. include:: ../include/links.rst

.. _moircs:

*************
Subaru MOIRCS
*************

Overview
========

This file summarizes several instrument-specific settings for the
Subaru/MOIRCS spectrograph.

MOIRCS (Multi-Object InfraRed Camera and Spectrograph) is a near-infrared
imager and multi-object spectrograph covering 0.9--2.5 µm.  The beam is split
into two optically independent channels, each with its own camera, grism and
2048×2048 Hawaii-2RG detector (0.116 arcsec/pixel).  The field is split along
the dispersion direction, so each detector sees a different part of the slit
mask.

.. warning::

    MOIRCS support is still under development.  So far, only MOS mode with
    the ``HK500`` and ``VB_K`` grisms has been tested: ``HK500`` with the
    data in the PypeIt development suite and ``VB_K`` with one mask
    observed in 2026.  Please report problems and send data for other
    grisms to the PypeIt developers.

Detectors and files
===================

MOIRCS writes each detector of an exposure to its own FITS file: chip 1 to
the odd frame number and chip 2 to the next (even) number, with the same
``EXP-ID`` header card.  For example, ``MCSP00237323.fits`` (chip 1) and
``MCSP00237324.fits`` (chip 2).

PypeIt treats the pair as **one exposure with two detectors** (DET01 and
DET02):

- Only the chip-1 files are listed in the :ref:`pypeit_file`.
  :ref:`pypeit_setup` removes the chip-2 files from the metadata table.
- When DET02 is reduced, PypeIt reads the chip-2 file automatically.  It
  must be in the same directory as the chip-1 file, and its ``DET-ID`` and
  ``EXP-ID`` are checked.  A missing or mismatched chip-2 file raises an
  error.
- The reduced spec2d and spec1d files hold both detectors.
- To reduce only one detector, set ``detnum`` in the :ref:`pypeit_file`,
  e.g.:

  .. code-block:: ini

      [rdx]
          spectrograph = subaru_moircs
          detnum = 2

The two detectors are *not* combined into a mosaic.  Because the channels
are optically independent, slit edges, tilts and wavelength solutions are
determined separately for each detector.

Multi-extension files
---------------------

Recent MOIRCS files (e.g. ``MCSA00352031.fits``, taken in 2026) hold more
than one HDU: the primary HDU is the processed image (the final minus the
initial Fowler reads), and the extensions are the individual reads.  A file
has 3 HDUs for a single Fowler sample (``DET-NSMP = 1``) and 21 HDUs for
``DET-NSMP = 10``.  PypeIt reads only the primary HDU, for both detectors.

Read noise
----------

MOIRCS uses Fowler sampling, and the number of samples (``DET-NSMP``)
depends on the frame: e.g. 10 for science frames and 1 for flats, arcs and
short standard-star exposures.  The read noise is set from the header of
each file as :math:`17.5/\sqrt{N_{\rm SMP}}` e-, i.e. 5.53 e- for
``DET-NSMP = 10`` and 17.5 e- for ``DET-NSMP = 1``.  ``DET-NSMP = 10`` is
assumed if the card is missing.

.. note::

    The detector values (gain, read noise, dark current and saturation)
    are still to be confirmed by the instrument scientist.  The read noise
    measured in the reference pixels of the ``VB_K`` test data is lower
    (about 11--14 e- for one sample and 5--6 e- for ten), and it does not
    scale as :math:`1/\sqrt{N_{\rm SMP}}`.

PypeIt File
===========

Run :ref:`pypeit_setup` with the ``-b`` flag so that the ``calib``,
``comb_id`` and ``bkg_id`` columns are added for A-B differencing:

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/raw/ -b -c A

Here is the :ref:`data_block` for the development-suite ``HK500`` data
(some columns removed for clarity):

.. code-block:: console

             filename |                 frametype |       target | dispname |        decker | exptime | lampstat01 | dithpat | dithpos | dithoff | calib | comb_id | bkg_id
    MCSP00237323.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |        off |   LINE2 |       A |     1.5 |     0 |       1 |      2
    MCSP00237325.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |        off |   LINE2 |       B |    -1.5 |     0 |       2 |      1
    MCSP00237177.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237179.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237181.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237163.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237165.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237167.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1

For ``VB_K`` data with standard stars and ThAr arcs, see the
:ref:`moircs_howto`.  Rows that :ref:`pypeit_setup` cannot type (e.g. the
ThAr arcs; see :ref:`moircs_howto_vbk_thar`) are written commented out
(``# MCSA... | None | ...``) and are ignored by :ref:`run-pypeit`.

.. _moircs_config_report:

Configuration
=============

A configuration is defined by the grism (``dispname``, from ``DISPERSR``),
the slit mask (``decker``, from ``SLIT``) and the binning (``BIN-FCT1/2``).
The mask image is taken without a grism (``DISPERSR = HOLE``), so it forms
its own configuration.

.. _moircs_frames_report:

Frames
======

Frame types are set from the ``DATA-TYP`` and ``OBJECT`` header cards and
the exposure time:

=======================================================  ===========================================
Header values                                            PypeIt frame type
=======================================================  ===========================================
``DATA-TYP = OBJECT``, exptime > 20 s                    ``science``, ``arc``, ``tilt``
``DATA-TYP = OBJECT``, exptime ≤ 20 s                    ``standard``, ``arc``, ``tilt``
``DATA-TYP = STANDARD_STAR``, exptime > 20 s             ``science``
``DATA-TYP = STANDARD_STAR``, exptime ≤ 20 s             ``standard``
``DATA-TYP = DOMEFLAT``, ``OBJECT = DOMEFLAT``           ``pixelflat``, ``illumflat``, ``trace``
``DATA-TYP = DOMEFLAT``, ``OBJECT = DOMEFLAT_OFF``       ``lampoffflats``
``DATA-TYP = DOMEFLAT``, ``OBJECT = MASKIMAGE``          not typed (ignored)
``DATA-TYP = INSTFLAT``, ``OBJECT = TH-AR``              not typed (ignored); see :ref:`moircs_howto_vbk_thar`
=======================================================  ===========================================

Standard stars are separated from science frames by their exposure time
only, whether their ``DATA-TYP`` is ``OBJECT`` or ``STANDARD_STAR`` (as in
the ``VB_K`` data); check the frame types in the
:ref:`pypeit_file`.  ``STANDARD_STAR`` frames are not used as arc or tilt
frames.

Lamp-off flats
--------------

The lamp-off dome flats are subtracted from the lamp-on flats.  This
removes the thermal emission from the telescope and dome, which is
significant in the *K* band.

.. _moircs_wavecalib:

Wavelength Calibration
======================

Wavelengths and tilts are measured from the OH sky lines in the science
frames.  The default line list is ``OH_NIRES``.

For ``HK500``, the default method is ``reidentify``, using the archive
``subaru_moircs_HK500.fits``.  This holds the OH spectra and wavelength
solutions of 19 slits from the development-suite mask.  Each slit of a new
mask covers a different part of the spectrum, depending on its position in
the mask, and reidentifying against many archived slits handles this well.
On the development-suite data, every science slit on both detectors is
calibrated, with an RMS of 0.24--0.43 pixels.

The sky emission in ``HK500`` spectra covers roughly 1.3--2.3 µm; the
wavelength solution is an extrapolation outside that range.

``VB_K``
--------

``VB_K`` (nominally 1.8--2.5 µm, about 1.94 Å/pixel) disperses past the
detector edges, so each slit sees only part of the band, set by its
position in the mask.  The defaults are:

- ``method = reidentify`` with the archive ``subaru_moircs_VB_K.fits``,
  which holds the OH spectra and solutions of 28 slits (both detectors) of
  one mask and covers 1.90--2.54 µm;
- ``lamps = OH_MOSFIRE_K`` and ``fwhm = 6.0``: the lines are about 6 pixels
  wide with slits of about 0.7 arcsec (R ~ 1900);
- ``match_toler = 0.75`` (pixels).  The archive spectra overlap each slit
  only in part, and with the default tolerance (2 pixels) lines
  misidentified at a slit end can survive the fit.  In a leave-one-out
  test, this gave a solution with a normal RMS (0.36 pixels) that was off
  by up to 35 Å at its red end.

On the test mask, all 14 science slits on each detector are calibrated with
41--56 lines and an RMS of 0.06--0.24 pixels.  Two independent checks agree
with these solutions: the Ar lines in ThAr arcs taken through the same mask
(to about 0.1 pixel per slit, after a flexure shift; see
:ref:`moircs_flexure`), and the slit positions in the mask image.

Things to keep in mind:

- The archive is built from the same mask it was tested on.  Other masks
  have not been tested; check the solutions of the slits at the ends of the
  band (see :ref:`moircs_howto`).  A wrong solution can have a normal RMS
  and be off by tens of Å at one end of the slit.
- Slits that reach below 1.90 µm have no archive match.  The
  ``holy-grail`` method (``[calibrations][wavelengths] method =
  holy-grail``) also calibrates every slit of the test mask and can be
  used instead.
- The OH line lists end at 2.499 µm, so the reddest part of some slits
  (up to ~250 pixels) is extrapolated.  The light falls to half its peak
  at about 2.47 µm.

Other grisms use the ``holy-grail`` method until templates are built.

.. _moircs_flexure:

Flexure
=======

Spectral flexure correction is turned off (``spec_method = skip``) because
the wavelength solution comes from the science frames themselves.

Frames taken at other pointings, however, are shifted relative to the
science frames.  In the ``VB_K`` test data, the standard stars (taken 2.3 h
after the science frames) and the ThAr arcs (taken 1.2 h before) are both
displaced by up to ~1 pixel (about 2 Å), with a shift plus a rotation of
0.3--0.6 pixel per 1000 pixels across each detector.  This matters for
standard stars, which use the science frames' wavelength solution by
default (see :ref:`moircs_sensitivity`).

Bad-pixel mask
==============

A static bad-pixel mask is applied for each detector.  It is derived from
the NAOJ MOIRCS bad-pixel masks for the detectors installed in 2016
(``mcsbadpix_oct2016``).  The beam-splitter shadow included in those masks
applies to imaging only and has been removed; about 0.07% of the pixels of
each detector are masked.

The mask still matches the ``VB_K`` flats taken in 2026 where it masks
(89--96% of the masked pixels are outliers in the flats), but about
720--760 newer defects per detector (~0.02% of the pixels) are not in it.
Static defects mostly cancel in A-B frames.

Slit masks
==========

Alignment-star boxes are traced as short slits (about 4.4 arcsec).  They
are flagged as box slits through ``minimum_slit_length_sci`` and are not
reduced as science slits.  Slits that overlap other slits in the spatial
direction (e.g. an alignment box next to a science slit at a different
position along the dispersion axis) may be merged or lost.

The lamp-on dome flats can saturate in the alignment-star boxes (in the
``VB_K`` test data, most of the box pixels are above the 33000 ADU
saturation level).  This does not affect the edge tracing, and the boxes
are not reduced.

The mask-design files are not yet used, so objects are not matched to
their targets in the mask design.

Background Subtraction
======================

The science frames are taken with the ``K_DITPAT``/``K_DITCNT``/
``K_DITWID`` dither cards.  For the two-position ``LINE2`` pattern,
position 1 is A and position 2 is B, with offsets of ±``K_DITWID``/2.
Each A frame is paired with the closest-in-time B frame for background
subtraction, and vice versa.  See :doc:`../A-B_differencing`.

Science and standard frames are paired separately, so several A-B
standard-star sets (e.g. one per slit position) are each paired within
themselves.

Single-sample frames (e.g. 15 s standards with ``DET-NSMP = 1``) show
bands of correlated read noise along the raw rows, i.e. along the
dispersion direction.  They leave faint residuals along the spectra, which
the object finding may pick up as faint objects.
For ``VB_K``, only one object per slit is kept in standard-star frames
(``[reduce][findobj] maxnumber_std = 1``), so that the standard star is
the object used for the sensitivity function.

.. _moircs_sensitivity:

Fluxing
=======

The default sensitivity function uses the ``IR`` algorithm with the PCA
telluric model (``TellPCA_3000_26000_R10000.fits``), as for the other
near-infrared spectrographs.  This has been tested on ``VB_K`` standards
only.

Standard stars
--------------

MOIRCS standards are taken through a slit of the science mask, so they
share the science configuration and calibration group.  By default they
also use the wavelength solution of the science frames (they are not
``arc`` frames; see :ref:`moircs_frames_report`).  Because of the flexure
between pointings (see :ref:`moircs_flexure`), the standards' own OH lines
can be offset by ~1 pixel.  In the ``VB_K`` test data, the telluric model
of the ``IR`` algorithm fits a shift of about this size, so the effect on
the sensitivity function is small.  A separate calibration
group for the standards (an optional manual step) is described in the
:ref:`moircs_howto`.

For ``VB_K``, each slit sees only part of the band, so the observer may
place the standard in two (or more) slits at different positions in the
mask.  Coadd the A and B frames of each set with
:ref:`pypeit_coadd_1dspec`, then pass the coadds to :ref:`pypeit_sensfunc`,
which splices sensitivity functions that cover different wavelength
ranges.

- The coadds must have the same length (same wavelength step and number of
  pixels): :ref:`pypeit_sensfunc` cannot yet splice inputs of different
  lengths.  Use ``wave_method = linear`` with the same ``dwave`` and the
  same ``wave_grid_max - wave_grid_min`` for every set.
- Set ``star_mag`` to the **V** magnitude of the star, which is how
  :ref:`pypeit_sensfunc` uses it, even though the data are in the *K*
  band.
- The pixel flat is normalized along each slit, so the flat-fielded counts
  keep each slit's throughput.  Standards in different slits can therefore
  differ by up to ~10% near the detector edges.  In the ``VB_K`` test
  data, this leaves a 0.14 mag step in the spliced sensitivity function,
  where it switches from one standard to the other.

Detector 2
----------

:ref:`pypeit_flux_calib` applies one sensitivity function to every object
in a spec1d file, whatever its detector.  A standard measured on one
detector is therefore also applied to the other one.  The two MOIRCS
channels have different throughputs: for ``VB_K``, the dome flats and the
night sky are 21--25% brighter on detector 2 than on detector 1, so
detector-2 spectra fluxed with a detector-1 standard are 21--25% too
bright.

.. warning::

    Transferring a sensitivity function between detectors is
    **experimental**, and the default is still under discussion with the
    instrument scientist.  The tools below are not part of PypeIt yet.

Two transfer modes have been tested on ``VB_K`` data, with a development
script in the PypeIt development suite
(``pypeitdev/subaru_moircs_vbk/transfer_sensfunc.py``):

- ``direct``: the zeropoint of the other detector equals that of the
  standard's detector, optionally scaled by a constant (e.g. 1.21 for
  ``VB_K``, from the dome flats);
- ``flatratio``: the zeropoint is scaled by the ratio of the lamp-on minus
  lamp-off dome-flat spectra of the two detectors, F2/F1(λ), where both
  detectors have slits at the same wavelength (it is left undefined
  elsewhere).

For ``VB_K``, F2/F1 is 1.21 to within ±3% over 1.95--2.50 µm, so the two
modes agree to a few percent.  The night sky gives a 3% higher ratio
(1.25), probably because each channel sees a different part of the dome
screen.  To apply a different sensitivity function to each detector, the
script splits each spec1d file into one file per detector, which
:ref:`pypeit_flux_calib` and :ref:`pypeit_coadd_1dspec` accept.

References
==========

- Ichikawa, T. et al. 2006, Proc. SPIE, 6269, 626916
- Suzuki, R. et al. 2008, PASJ, 60, 1347
