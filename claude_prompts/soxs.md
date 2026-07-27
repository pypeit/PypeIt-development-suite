# Develop the NTT/SOXS spectrograph

## Goals

We wish to develop the NTT/SOXS spectrograph for PypeIt.  

## Claude

### Skills

Consider using the skills in PypeIt-development-suite/.claude/skills/

## Context

Here are the context documents:
- Webpage for the instrument: https://www.eso.org/public/teles-instr/lasilla/ntt/soxs/
- SOXS pipeine: /home/xavier/Projects/PypeIt/Others/soxspipe

The PypeIt pipeline

## Coding

Here are guidelines for coding: 

- Use Python where possible
- Add inline comments to explain the effort
- Reuse existing code when possible
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods
- Add "Created by JXP and Claude" to any new methods.

## Testing

If you need to test the code: 

- use the files in PypeIt-development-suite/pypeitdev/ntt_soxs/
- run on the "pypeit14b" conda environment.

## Prompts

### Spectrograph Class

1. Read the context documents to understand the instrument and the code.  Your first task is to create a new spectrograph class in PypeIt/pypeit/spectrographs/ntt_soxs.py.  Model it after the other spectrograph classes in that folder.  Use the files in PypeIt-development-suite/pypeitdev/ntt_soxs/ as a guide.  Before doing so, ask me questions about the instrument and the code in Q&A below.  Use Fable if you can.  Log your work

## Q&A

Questions from Claude (2026-07-27), before writing `pypeit/spectrographs/ntt_soxs.py`.
Each has a proposed default — confirm or correct.

1. **Class structure.** I plan to model on `vlt_xshooter.py`: a `NTTSOXSSpectrograph`
   base class with two subclasses, `ntt_soxs_vis` and `ntt_soxs_nir`, both
   `pypeline = 'Echelle'` with `ech_fixed_format = True`. SOXS has no UVB arm and the
   acquisition camera is out of scope. OK?

2. **VIS arm as echelle.** The VIS arm is not a classical cross-dispersed echelle: it
   has 4 quasi-orders from stacked ion-etched gratings (per the soxspipe spectral
   format table: order 1 = 538–675 nm, 2 = 347–432, 3 = 431–541, 4 = 671–850), so
   wavelength is not monotonic with order number. I still plan to treat it as
   fixed-format Echelle with 4 orders. Agree, or would you rather treat VIS as
   MultiSlit with 4 "slits"?

3. **Git branch.** The PypeIt repo is currently on a branch named `soxs`, but it looks
   like it contains merges from the `state` (Binospec) work. Should I commit the new
   class to this `soxs` branch, or create a fresh branch off `develop`?

4. **Detector parameters.** Header values and soxspipe disagree slightly. soxspipe
   (from commissioning): VIS gain 1.0 e/ADU, RON 3.8 e; NIR gain 2.3 e/ADU, RON 7.3 e
   (header NIR cards are placeholders: GAIN=1.0, RON=0.0; header VIS says GAIN=0.909,
   CONAD=1.1, RON=3.3). Default: adopt the soxspipe commissioning values. OK?

5. **VIS readout formats.** All sample VIS frames are windowed 858×4096 (50-px
   prescan + 808 science columns of the 2048×4096 CCD44-82), 1x1 binning, 'Slow High
   Gain' readout. soxspipe also ships bad-pixel maps for 1x2, 2x1, 2x2 binning.
   Should the class support arbitrary binning/readout from the headers, or is
   858×4096 1x1 'Slow High Gain' the fixed operational mode I can assume?

6. **NIR readout & calibrations.** NIR frames are float32 post-processed up-the-ramp
   images in ADU (DIT×NDIT, e.g. the flux standard has DIT=150 s, NDIT=2). I'll take
   exptime = DIT and treat the frame as a single sampled-up-the-ramp image, X-shooter
   NIR style (no bias/overscan; darks used instead). The sample data contain **no
   bias or dark frames** — will those be added to `pypeitdev/ntt_soxs/`, or should
   the defaults be set up to reduce without them (e.g. NIR flats/arcs are not
   lamp-off subtracted here)?

7. **Wavelength solutions.** Lamps: VIS = ThAr, NIR = Ar/Hg/Ne/Xe (header
   `INS COMB LAMPS` is empty, so lamps are inferred from `DPR TYPE`/templates, as for
   X-shooter). No PypeIt reid templates exist for SOXS yet. For this first prompt I
   plan to wire plausible `reid_arxiv` placeholders + `method = 'holy-grail'`
   fallback, and build proper templates with `pypeit_identify` in a later prompt
   (per the `wavelength-calibration` skill). OK?

8. **soxspipe path.** The context path `/home/xavier/Projects/PypeIt/Others/soxspipe`
   does not exist on this machine. I cloned soxspipe (thespacedoctor/soxspipe) into my
   session scratchpad for reference. Want me to clone it permanently to the path named
   in this document (creating `Others/`), and should I update this document either way?

9. **Order tracing seed.** For the fixed-format methods (`order_spat_pos`, `orders`,
   `spec_min_max`) I need measured order positions. I plan to derive them by running
   edge tracing on the sample flats and/or from the soxspipe static tables
   (`SOXS_SPECTRAL_FORMAT_TAB_{VIS,NIR}.fits`, order centre/edge polynomials). Any
   preference or existing measurements I should use?

10. **Housekeeping.** `pypeitdev/ntt_soxs/{NIR,VIS}/` contain macOS AppleDouble
    droppings (`._SOXS.*.fits`) that are not valid FITS and may trip up
    `pypeit_setup` globs. OK to delete them?

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-07-27 (Context gathering for the SOXS spectrograph class; questions posted in Q&A)

Executed prompt 1 of "Spectrograph Class" up to its ask-first checkpoint: read the
context documents, surveyed the sample data, and posted 10 questions in the Q&A
section above. No code written yet, pending answers. What was learned:

- **Instrument** (ESO page + soxspipe): SOXS is a two-arm, R~4500 échelle
  spectrograph at the NTT Nasmyth focus, 350–2000 nm total coverage, built for
  transient follow-up. Both arms expose simultaneously (science frames are
  `OBJECT,ASYNC`); slit length 12 arcsec.
- **VIS arm**: e2v CCD44-82, raw frames windowed to 858×4096 (50-px prescan +
  808 science columns), 1x1, 'Slow High Gain' readout. 4 quasi-orders from stacked
  ion-etched gratings covering 347–850 nm; wavelength is NOT monotonic with order
  number. Orders run along the 4096-px axis (verified on a sample flat). Arc lamp
  is ThAr. soxspipe adopts gain 1.0 e/ADU, RON 3.8 e, saturation 65535.
- **NIR arm**: Teledyne Hawaii-2RG, 2048×2048, up-the-ramp readout delivered as
  float32 images in ADU (DIT/NDIT in `ESO DET SEQ1 DIT` / `ESO DET NDIT`); 4-px
  reference-pixel border. 15 orders (m=10–24) covering 800–2033 nm. Arc lamps
  Ar/Hg/Ne/Xe. soxspipe adopts gain 2.3 e/ADU, RON 7.3 e (commissioning values;
  the header GAIN/RON cards are placeholders).
- **Headers/frame typing**: standard ESO hierarchy. Arm in `ESO SEQ ARM`; frame
  types from `ESO DPR CATG/TYPE/TECH` (`LAMP,FLAT`, `LAMP,WAVE`, `STD,FLUX` with
  NODDING, `OBJECT` with STARE); slit in `ESO INS VISE ID` (VIS) / `ESO INS NISE
  NAME` (NIR) — arcs use SLIT0.5/UNDEFINED/BLANK. `ESO INS COMB LAMPS` is empty in
  all sample frames, so lamp identity cannot be read from the header. Sample data
  (7 frames/arm) contain no bias or dark frames.
- **Code**: no SOXS work exists in PypeIt yet. `vlt_xshooter.py` is the natural
  template (fixed-format echelle, ESO headers, NIR arm conventions); the
  `new-spectrograph` dev-suite skill documents the full procedure. The NTT
  telescope is already defined in `pypeit/telescopes.py` (used by `ntt_efosc2`).
- **soxspipe**: the path given in Context does not exist on this machine; cloned
  from GitHub (thespacedoctor/soxspipe) into the session scratchpad. Key resources:
  `soxspipe/resources/soxs_detector_parameters.yaml` (detector params, dispersion
  axes, line atlases) and
  `soxspipe/resources/static_calibrations/soxs/SOXS_SPECTRAL_FORMAT_TAB_{VIS,NIR}.fits`
  (order wavelength ranges, listed in the Q&A above).
- **Data gotcha**: `pypeitdev/ntt_soxs/{NIR,VIS}/` contain macOS AppleDouble files
  (`._SOXS.*.fits`) that are not valid FITS.