.. include:: ../include/links.rst

.. _moircs_howto:

===================
Subaru-MOIRCS HOWTO
===================

Overview
========

This tutorial walks through the reduction of Subaru/MOIRCS multi-object
spectroscopy with PypeIt.  The first part uses the ``HK500`` data in the
PypeIt development suite (``RAW_DATA/subaru_moircs/HK500``).  The second
part, :ref:`moircs_howto_vbk`, covers the ``VB_K`` grism, standard stars,
fluxing and ThAr arcs.  See :ref:`moircs` for the instrument-specific
settings.

.. warning::

    MOIRCS support is still under development; only the ``HK500`` and
    ``VB_K`` grisms in MOS mode have been tested so far.

If this is your first time using PypeIt, we recommend starting with the
:doc:`Shane Kast tutorial<kast_howto>`.

The data
========

The example data are one MOS mask (``MO17A_COSMOS2``) observed with the
``HK500`` grism:

- 3 lamp-on dome flats and 3 lamp-off dome flats;
- one mask image (not used);
- 2 science exposures of 180 s, taken as an A-B pair with a 3 arcsec
  dither (``K_DITPAT = LINE2``).

Each exposure is written as two files, one per detector: e.g.
``MCSP00237323.fits`` (chip 1) and ``MCSP00237324.fits`` (chip 2).  Keep the
two files of each exposure in the same directory.  PypeIt lists only the
chip-1 file and reads the chip-2 file itself when it reduces DET02.

Setup
=====

Run :ref:`pypeit_setup` with the ``-b`` flag, which adds the ``calib``,
``comb_id`` and ``bkg_id`` columns used for A-B differencing:

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/RAW_DATA/subaru_moircs/HK500 -b -c A

You will see warnings that the chip-2 files (and the mask-image frame) are
removed from the table; this is expected.  The resulting
``subaru_moircs_A/subaru_moircs_A.pypeit`` file has this :ref:`data_block`
(some columns removed for clarity):

.. code-block:: console

             filename |                 frametype |       target | dispname |        decker | exptime | dithpat | dithpos | dithoff | calib | comb_id | bkg_id
    MCSP00237323.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |   LINE2 |       A |    -1.5 |     0 |       1 |      2
    MCSP00237325.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |   LINE2 |       B |     1.5 |     0 |       2 |      1
    MCSP00237177.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237179.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237181.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237163.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237165.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237167.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1

Check that:

- the science frames are typed ``arc,science,tilt``; the OH lines in the
  science frames provide the wavelength calibration;
- the lamp-off flats are typed ``lampoffflats``;
- each A frame has the B frame as its background (``bkg_id``) and vice
  versa.  For longer sequences, each frame is paired with the closest-in-time
  frame at the other dither position; edit ``comb_id`` and ``bkg_id`` if you
  want a different pairing (see :doc:`../A-B_differencing`).

Main Run
========

.. code-block:: console

    cd subaru_moircs_A
    run_pypeit subaru_moircs_A.pypeit -o

Both detectors are reduced in one run.  To reduce only one of them, add
``detnum`` to the ``[rdx]`` block, e.g. ``detnum = 2``.

Inspecting the outputs
======================

Slit edges
----------

.. code-block:: console

    pypeit_chk_edges Calibrations/Edges_A_0_DET01.fits.gz

For the example mask, edge tracing finds 20 slits on DET01 (17 science slits
and 3 alignment-star boxes) and 19 on DET02 (15 science slits and 4 boxes).
The boxes (~4.4 arcsec long) are flagged as box slits and are not reduced
as science slits.  Many slits cover only part of the detector in the
spectral direction, depending on their position in the mask; this is
expected.

Flat field
----------

.. code-block:: console

    pypeit_chk_flats Calibrations/Flat_A_0_DET01.fits

The lamp-off flats are subtracted from the lamp-on flats before they are
used.

Wavelengths
-----------

.. code-block:: console

    pypeit_chk_wavecalib Calibrations/WaveCalib_A_0_DET01.fits

For ``HK500``, the OH lines are reidentified against an archive of solutions
from the example mask.  In the example data, every science slit on both
detectors is calibrated, with 50--72 lines and an RMS of 0.24--0.43 pixels
(the dispersion is ~7.6 Å/pixel).  The alignment-star boxes are not
calibrated.  The sky emission covers about 1.3--2.3 µm.

Science frames
--------------

.. code-block:: console

    pypeit_show_2dspec Science/spec2d_MCSP00237323-COSMOS2_MOIRCS_20170410T064929.408.fits --det 1

The A-B difference image shows each object as a positive trace with a
negative trace 26 pixels (3 arcsec) away.  The spec1d file of each exposure
holds the extracted spectra from both detectors:

.. code-block:: console

    pypeit_show_1dspec Science/spec1d_MCSP00237323-COSMOS2_MOIRCS_20170410T064929.408.fits

The object-finding threshold is lowered to ``snr_thresh = 5`` for MOIRCS;
faint targets may still need to be extracted manually (see
:doc:`../manual`).

Fluxing and coadding
====================

The ``HK500`` example has no standard star.  Fluxing (``IR`` algorithm with
the PCA telluric model) and coadding follow the standard PypeIt procedures;
see :doc:`../fluxing`, :doc:`../telluric` and :doc:`../coadd1d`.  Fluxing
has been tested with ``VB_K`` standards only; see
:ref:`moircs_howto_vbk_flux` for the MOIRCS-specific steps.

.. _moircs_howto_vbk:

VB_K example
============

The data
--------

The ``VB_K`` example is one MOS mask (``MO_CC0958PA200_1``) observed in
2026:

- 7 lamp-on and 7 lamp-off dome flats (7 s);
- 2 science exposures of 180 s, an A-B pair with a 3.1 arcsec dither;
- 2 standard-star sets (HIP59174, 15 s), each an A-B pair, with the star in
  a different slit of the science mask (see :ref:`moircs_howto_vbk_flux`);
- 4 ThAr arcs (3 s), taken through the mask in the evening;
- one mask image (not used).

The files (``MCSA*.fits``) hold several HDUs; PypeIt reads the primary HDU,
which is the processed image (see :ref:`moircs`).  As for ``HK500``, keep
the chip-1 and chip-2 files of each exposure in the same directory.

Setup
-----

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/VB_K/raw/ -b -c A

``-r`` accepts several directories if the frames are spread over more than
one.  The mask image (``DISPERSR = HOLE``) forms its own setup and is not
written with ``-c A``.  The :ref:`data_block` is (some columns removed for
clarity):

.. code-block:: console

               filename |                 frametype |         target | dispname |           decker | exptime | dithpat | dithpos | dithoff | calib | comb_id | bkg_id
    # MCSA00351935.fits |                      None |          TH-AR |     VB_K | MO_CC0958PA200_1 |     3.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    # MCSA00351937.fits |                      None |          TH-AR |     VB_K | MO_CC0958PA200_1 |     3.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    # MCSA00351939.fits |                      None |          TH-AR |     VB_K | MO_CC0958PA200_1 |     3.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    # MCSA00351941.fits |                      None |          TH-AR |     VB_K | MO_CC0958PA200_1 |     3.0 |    none |    none |     0.0 |     0 |      -1 |     -1
      MCSA00352031.fits |          arc,science,tilt | CC0958_PA200_1 |     VB_K | MO_CC0958PA200_1 |   180.0 |   LINE2 |       A |   -1.55 |     0 |       1 |      2
      MCSA00352033.fits |          arc,science,tilt | CC0958_PA200_1 |     VB_K | MO_CC0958PA200_1 |   180.0 |   LINE2 |       B |    1.55 |     0 |       2 |      1
      MCSA00351909.fits |              lampoffflats |   DOMEFLAT_OFF |     VB_K | MO_CC0958PA200_1 |     7.0 |    none |    none |     0.0 |     0 |      -1 |     -1
      ...
      MCSA00351921.fits |              lampoffflats |   DOMEFLAT_OFF |     VB_K | MO_CC0958PA200_1 |     7.0 |    none |    none |     0.0 |     0 |      -1 |     -1
      MCSA00351895.fits | pixelflat,illumflat,trace |       DOMEFLAT |     VB_K | MO_CC0958PA200_1 |     7.0 |    none |    none |     0.0 |     0 |      -1 |     -1
      ...
      MCSA00351907.fits | pixelflat,illumflat,trace |       DOMEFLAT |     VB_K | MO_CC0958PA200_1 |     7.0 |    none |    none |     0.0 |     0 |      -1 |     -1
      MCSA00352125.fits |                  standard |       HIP59174 |     VB_K | MO_CC0958PA200_1 |    15.0 |   LINE2 |       A |    -2.4 |     0 |       3 |      4
      MCSA00352127.fits |                  standard |       HIP59174 |     VB_K | MO_CC0958PA200_1 |    15.0 |   LINE2 |       B |     2.4 |     0 |       4 |      3
      MCSA00352141.fits |                  standard |       HIP59174 |     VB_K | MO_CC0958PA200_1 |    15.0 |   LINE2 |       A |    -2.5 |     0 |       5 |      6
      MCSA00352143.fits |                  standard |       HIP59174 |     VB_K | MO_CC0958PA200_1 |    15.0 |   LINE2 |       B |     2.5 |     0 |       6 |      5

Check that:

- the standards (``DATA-TYP = STANDARD_STAR``, 15 s) are typed
  ``standard``.  Standards and science frames are separated by exposure
  time only (standard ≤ 20 s);
- each standard set is paired within itself (``bkg_id``);
- the ThAr arcs are untyped and commented out (see
  :ref:`moircs_howto_vbk_thar`).

All frames are in calibration group 0, so the standards use the flats,
slit edges and OH wavelength solution of the science frames.

Main run
--------

.. code-block:: console

    cd subaru_moircs_A
    run_pypeit subaru_moircs_A.pypeit -o

Inspecting the outputs
----------------------

- **Slit edges.**  Edge tracing finds 19 slits on each detector: 14
  science slits and 5 alignment-star boxes, which are flagged as box slits.
  The lamp-on flats saturate in the boxes; this does not affect the
  tracing.
- **Wavelengths.**  ``VB_K`` uses ``reidentify`` with the
  ``subaru_moircs_VB_K.fits`` archive and the ``OH_MOSFIRE_K`` lines.  All
  14 science slits on each detector are calibrated, with 41--56 lines and
  an RMS of 0.06--0.24 pixels (~1.94 Å/pixel).  Each slit covers about
  0.4 µm, between 1.90 and 2.54 µm depending on its position in the mask.
  For a new mask, check with ``pypeit_chk_wavecalib`` and the sky spectra
  that the slits at the ends of the band are right: a wrong solution can
  have a normal RMS and still be off at one end.  If the archive fails,
  set ``method = holy-grail`` in the ``[calibrations][wavelengths]``
  block.
- **Objects.**  The science pair has one faint target (DET01).  In the
  standard frames, the star is in DET01 slit 337 (first set, 2.14--2.54
  µm) and slit 469 (second set, 1.90--2.29 µm), with S/N ~200; the DET02
  frames contain sky only.  Only one object per slit is kept in standard
  frames (``maxnumber_std = 1``).

.. _moircs_howto_vbk_flux:

Standard stars and fluxing
--------------------------

``VB_K`` disperses past the detector edges, so one slit covers only part
of the band.  The standard was observed in two slits so that, together,
the two sets cover 1.90--2.54 µm.

1. Coadd the A and B spectra of each set without fluxing, on linear
   wavelength grids of the same step and length (needed for the splice):

   .. code-block:: ini

       [coadd1d]
         coaddfile = coadd1d_std1_HIP59174.fits
         flux_value = False
         wave_method = linear
         dwave = 1.94
         wave_grid_min = 21430.0
         wave_grid_max = 25368.2

       coadd1d read
         path Science
         filename | obj_id
         spec1d_MCSA00352125-HIP59174_MOIRCS_20260529T082520.088.fits | SPAT0360-SLIT0337-DET01
         spec1d_MCSA00352127-HIP59174_MOIRCS_20260529T082555.088.fits | SPAT0318-SLIT0337-DET01
       coadd1d end

   and the same for the second set (``SLIT0469``, ``wave_grid_min =
   18995.0``, ``wave_grid_max = 22933.2``).  Run
   :ref:`pypeit_coadd_1dspec` on each file.

2. Compute the sensitivity function from both coadds.  The sens file
   gives the star's spectral type and **V** magnitude:

   .. code-block:: ini

       [sensfunc]
         star_type = A0
         star_mag = 7.45

   .. code-block:: console

       pypeit_sensfunc coadd1d_std2_HIP59174.fits coadd1d_std1_HIP59174.fits -s HIP59174.sens -o sens_HIP59174.fits

   :ref:`pypeit_sensfunc` splices the two sensitivity functions.  In this
   example the spliced zeropoint has a 0.14 mag step at 2.143 µm, where it
   switches from one standard to the other.  It comes from the different
   throughputs of the two slits, which the normalized pixel flat leaves in
   the counts.  The zeropoint is 18.3--18.55 mag over 1.95--2.45 µm.

3. Flux the science spectra with :ref:`pypeit_flux_calib` and coadd them
   with :ref:`pypeit_coadd_1dspec`, as usual.  :ref:`pypeit_flux_calib`
   overwrites the spec1d files, so keep a copy if you want to flux them
   again.

The standards use the OH wavelength solution of the science frames, taken
2.3 h earlier.  Their own OH lines are offset by about 1 pixel (about
-1 pixel in the two star slits), because of flexure between pointings.
The telluric model fits a shift of about this size.

.. warning::

    :ref:`pypeit_flux_calib` applies the sensitivity function to the
    objects of both detectors.  For ``VB_K``, a detector-1 standard makes
    detector-2 fluxes 21--25% too bright.  Transferring the sensitivity
    function to detector 2 is experimental and not yet part of PypeIt;
    see :ref:`moircs_sensitivity`.

Optional: a calibration group for the standards
-----------------------------------------------

To give the standards their own OH wavelength solution and tilts (instead
of the science frames' solution), edit the :ref:`data_block` by hand: type
the standards ``arc,standard,tilt`` in calibration group 1, and put the
flats in both groups (see :ref:`calibration-groups`):

.. code-block:: console

               filename |                 frametype | ... | calib | comb_id | bkg_id
      MCSA00352031.fits |          arc,science,tilt | ... |     0 |       1 |      2
      MCSA00352033.fits |          arc,science,tilt | ... |     0 |       2 |      1
      MCSA00351909.fits |              lampoffflats | ... |   0,1 |      -1 |     -1
      ...
      MCSA00351895.fits | pixelflat,illumflat,trace | ... |   0,1 |      -1 |     -1
      ...
      MCSA00352125.fits |         arc,standard,tilt | ... |     1 |       3 |      4
      MCSA00352127.fits |         arc,standard,tilt | ... |     1 |       4 |      3
      MCSA00352141.fits |         arc,standard,tilt | ... |     1 |       5 |      6
      MCSA00352143.fits |         arc,standard,tilt | ... |     1 |       6 |      5

The 15 s standard frames have enough OH lines (about 40 per slit).  This
step has not yet been tested with PypeIt.

.. _moircs_howto_vbk_thar:

Optional: ThAr arcs
-------------------

The ThAr arcs (``DATA-TYP = INSTFLAT``, ``OBJECT = TH-AR``) are not typed
by default, for these reasons:

- they share the science configuration, and PypeIt combines every ``arc``
  frame of a calibration group, so they would be combined with the science
  frames, which are the default (OH) arcs;
- PypeIt has no *K*-band Th line list, and the Ar lines alone are sparse
  (4--11 per slit in this example);
- the OH lines are taken at the same time as the science, so they include
  the flexure, whereas the arcs are offset by up to ~1 pixel (see
  :ref:`moircs_flexure`).

:ref:`pypeit_setup` writes their rows commented out (``# MCSA... | None |
...``), and :ref:`run-pypeit` ignores them.  In the example, the Ar lines
agree with the OH solutions to about 0.1 pixel per slit, once that offset
is removed.

To use the arcs anyway (not tested with PypeIt):

1. Uncomment their rows and set their frame type to ``arc``.
2. Remove ``arc`` from the science frames, but keep ``tilt``: the OH lines
   trace the tilts better than the sparse ThAr lines.
3. Set the line lists (e.g. ``lamps = Ar_IR_MOSFIRE``) and a method that
   does not use the OH archive (e.g. ``method = holy-grail``) in the
   ``[calibrations][wavelengths]`` block.
4. Combine the arcs without sigma clipping:

   .. code-block:: ini

       [calibrations]
           [[arcframe]]
               [[[process]]]
                   clip = False

   The lamp lights a different part of the mask in each ThAr frame, so the
   default clipped combine rejects the bright frames on the lines.  In the
   example, this distorts the line profiles and moves line centroids by up
   to 2.7 pixels.
