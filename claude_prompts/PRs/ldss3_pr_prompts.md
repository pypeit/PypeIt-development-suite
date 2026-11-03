# LDSS3 PR

## Goals

Perform a PR on the LDSS3 code for PypeIt

## Context

There is an open PR on the PypeIt repository for the LDSS3 code. The PR is located at `https://github.com/pypeit/PypeIt/pull/2201`

## Prompts

1. Read this file. Examine the PR and place a draft of your comments in the Draft section below.  Use Fable if you can. Log your work.

## Draft
*Draft review of PR #2201 (head `ac9a085`, branch `joaquinhernandezg:ldss3_multigrating`). Prepared with a Fable subagent; line numbers refer to the PR-head versions of the files. Not posted to GitHub.*

### Summary

Thanks, Joaquín — this is a big improvement over what was there. The amplifier-file handling (one file per exposure in `pypeit_setup`, companion found automatically), the telescope/MJD/binning/`dataext` fixes, the readout-mode gain/RN tables with the measured Fast values, and the much more careful frame typing and docs are all welcome. I checked the amp joining/orientation, the regex/`find_raw_files` logic (including `.fits.gz`), the MJD and binning parsing, the base-class signatures (`default_pypeit_par` classmethod, `config_specific_par(inp, inp_par)`, `find_raw_files` classmethod, `log` rather than `msgs`) and the line-list format (matches `OH_MODS_lines.dat`; `NIST=1` is what `wv_fitting` needs) — those all look right. There are a few things that need to change before merge, mostly around the saturation value and a couple of places where the new pieces contradict each other.

### Required changes

1. **`pypeit/spectrographs/magellan_ldss3.py:294` — saturation is double-scaled by gain.** `DetectorContainer.saturation` is in ADU (cf. `keck_deimos.py`, `magellan_mage.py`, both `65535.`). `PypeItImage.build_mask(saturation='default')` already multiplies `saturation * nonlinear * gain` per amplifier (via `map_detector_value`). With `65535 * min(gain)` the effective threshold becomes ~150–155k e⁻, above the physical 16-bit ceiling of 65535×gain ≈ 94–108k e⁻, so saturated pixels are never flagged. Please set `saturation = 65535.` and drop the comment block at L285–293. Related question: the Slow-mode gains of 0.16/0.19 e⁻/ADU (L98) look implausible — that would put the ADC ceiling at ~10k e⁻. Can you double-check the LCO table (units?)?

2. **`magellan_ldss3.py:588` vs `:187–190` — the "only c2 files" case is inconsistent.** `find_raw_files` deliberately keeps a lone `c2` file per exposure, and `get_rawimage` accepts it, but `check_frame_type` then requires `amp in ('1','None','')`, so a directory holding only `c2` files yields zero typed frames (verified with a synthetic table). Either drop the `primary` guard (`find_raw_files` already de-duplicates, and `get_rawimage` reads the companion regardless of which file is listed), or compute "primary" per exposure stem within `fitstbl`. Please also make the comment at L586–587 and the log message at L206–207 ("Only c1 files are read…") match.

3. **`magellan_ldss3.py:528` — `valid_configuration_values` can delete bias/dark frames.** `PypeItSetup.run` calls `clean_configurations()` *before* `get_frame_types()` (`pypeitsetup.py:367, 373`) and applies the `dispname` whitelist to every row, with no exception for config-independent frames. Any bias/dark whose `GRISM` is not exactly one of the three VPH names (e.g. `Open`) is removed before typing — contradicting the `config_independent_frames` docstring at L535–537 ("grism wheel left wherever it happened to be"). What `GRISM` values do your biases carry? If not always a VPH grism, drop `dispname` from `valid_configuration_values` (the `dispersed` mask in `check_frame_type` already keeps Open frames out of spectroscopic types).

4. **`magellan_ldss3.py:635–636` — frame typing regresses on the existing dev-suite MOS data.** The `magellan_ldss3_vph-red_mos.pypeit` flats have `OBJECT = 'ES1_34_QlQh'`; the old code matched `FlatQH`, but the new word list (`qh`, `ff`, `flatfield`) and substrings (`flat`, `quartz`, `dome`) do not match `qlqh`, so `pypeit_setup` now types those 4 s frames as `standard` (verified). Adding `'qh'` to `flat_substrings` (or `'qlqh'`/`'ql'` to `flat_words`) fixes it.

### Suggestions / questions

1. `magellan_ldss3.py:236` — when `hdu` is given but `amp_headers` is None (e.g. `get_detector_par(1)` from `datacube.py`, `ql.py`, or any caller passing an HDUList), the Fast-mode nominal values are used regardless of `SPEED`. Consider deriving from `hdu[0].header` in that branch. The two near-identical loops at L240–271 could be one loop or a small helper.

2. Single-amplifier fallback (`:704–706`, `:743`): a lone `c2` file is placed unflipped at the low-spatial edge and labelled amplifier 1 in `rawdatasec_img`, so the image is mirrored relative to the two-amp layout and the `bpm` column ranges (`:786–788`) land on the wrong columns. Is single-amp support needed? If yes, flip/label by `OPAMP`; if not, raising a `PypeItError` is simpler and safer than a warning.

3. `magellan_ldss3.py:155–163` / `:707–709` — if both `ccd0042c1.fits` and `ccd0042c1.fits.gz` exist, `amp_files` returns three entries and `get_rawimage` errors with "Expected 1 or 2 LDSS3 amplifier files" (verified). De-duplicate by amplifier index, preferring the extension of `raw_file`. (Nit: a stem containing `[`, `?` or `*` breaks `path.parent.glob`; `glob.escape` avoids it.)

4. `magellan_ldss3.py:814–825` — `x1`/`b1` are parsed but unused; the slices assume `DATASEC` starts at pixel 1 and `BIASSEC` follows immediately. Using the parsed ranges directly (`array[x1-1:x2]`, `array[b1-1:b2]`) removes that assumption.

5. Add a `raw_header_cards()` override returning e.g. `['GRISM', 'APERTURE', 'BINNING', 'FILTER', 'SPEED', 'OPAMP']` (as `keck_deimos.py` does), so spec1d/spec2d headers carry the config-identifying cards — `SPEED` especially, since gain/RN now depend on it.

6. Templates: `magellan_ldss3_VPH-ALL_7100.fits` is still in `pypeit/data/arc_lines/reid_arxiv/` but no longer referenced — please remove it in this PR. How were the four new templates built (`pypeit_identify` → `templates.py`?) — a short entry in `reid_arxiv/README` would help; please confirm they are the current format `full_template` expects. The sky template + `OH_LDSS3_vac` are documentation-only (not wired into any par), which is fine, but nothing exercises that they load.

7. Frame-typing edge cases worth a docs note: `'lamp'` is an arc word, so an `OBJECT` like `lamp flat` becomes `arc` (arcs win over flats at L614). And `science`/`standard` do not require `EXPTYPE == 'Object'`, so a `Bias`-typed frame with an odd name and `EXPTIME > 100` would become science.

8. Backward compatibility: `PypeItMetaData.merge` overwrites header-derived columns with the user table, so existing `.pypeit` files carrying `mjd` in JD keep those values until regenerated (JD "mjd" gives odd obstimes in output names). Also `filter1` is now a configuration key — dev-suite sets with calibrations and science through different filters will split into separate setups.

9. `config_specific_par` (`:465–479`): the three grism branches differ only in the filename — a dict lookup would be shorter; the `elif grating is not None` branch re-sets lamps/method already set in `default_pypeit_par`, so a warning alone suffices.

### Nits

- `doc/spectrographs/magellan_ldss3.rst:104` — section underline shorter than the title (Sphinx warning). `:16` "Currenly ony" → "Currently only". `:49–60` says matching is on delimited words, but `flat`, `quartz`, `dome` are substring matches and `arcs`, `ff`, `flatfield` are missing from the table — mirror `arc_words`/`flat_words`/`flat_substrings` (`.py:634–637`).
- `magellan_ldss3.py:6, 8` — module-docstring `:func:` refs should be fully qualified to resolve.
- `magellan_ldss3.py:721` — `det if det is not None else 1` is passed but `det` is unused in `get_detector_par`.
- `magellan_ldss3.py:29–41` — "LDS3-C" → "LDSS3-C".
- `doc/releases/2.1.0dev.rst:84–93` — fine under *Instrument-specific Updates*; the telescope/MJD/binning/`dataext` fixes could also get a one-liner under *Bug Fixes*.

### Before merge

- [ ] Fix saturation (Req. 1) and resolve Req. 2–4.
- [ ] Regenerate the five dev-suite pypeit files (`magellan_ldss3_{longslit,vph-red_longslit,vph-red_mos,vph-all_sci,vph-all_std}.pypeit`) with the new code; the two VPH-All files currently point at a user's Dropbox path. Only `VPH-ALL_Std`/`VPH-ALL_Sci` are registered in `test_scripts/setups.py`; add VPH-Blue and VPH-Red setups with raw data and register them.
- [ ] Run the dev suite for `magellan_ldss3` and check QA for edge tracing near the bpm columns and the amp-boundary step with the new gains.
- [ ] Add a lightweight unit test: `parse_amp_file`/`find_raw_files` on a temp dir (incl. `.fits.gz`), `check_frame_type` on a synthetic table, and that `OH_LDSS3_vac` and the four `reid_arxiv` files load.
- [ ] Docs build clean of new warnings.
- [ ] Confirm a gzipped `c1`/`c2` pair reduces end-to-end.

## Logs

### 2026-09-28 (Prompt 1: draft review of PR #2201)

- Fetched PR #2201 with `gh` (6 commits, +1098/−217; `magellan_ldss3.py` rewritten, docs, release notes, OH line list, four new `reid_arxiv` templates). No existing review comments on the PR.
- Saved the diff and PR-head files to the session scratchpad and had a Fable subagent review them against the PypeIt base classes, other spectrographs and the dev suite; it ran synthetic checks in `pypeit14` (scratch script, not committed). I re-verified the saturation (Req. 1) and `primary`-guard (Req. 2) findings directly.
- Learned: this checkout of PypeIt uses `log` (not `msgs`); `PypeItSetup.run` applies `valid_configuration_values` before frame typing, so a `dispname` whitelist can silently drop bias/dark frames; `DetectorContainer.saturation` is in ADU and gain-scaled downstream. Dev-suite LDSS3 coverage is only VPH-All (Sci/Std) in `test_scripts/setups.py`, and those pypeit files reference a Dropbox path.
- This prompt doc had no `## Logging` section, so I used the `### YYYY-MM-DD (summary)` format from `CLAUDE.md`.
- Nothing posted to GitHub; no git state changed.

### 2026-09-29 (Posted the PR #2201 review)

- At the user's request, posted the Draft section as a review on PR #2201 as `profxj`, with review type "Comment" (not "Request changes"): https://github.com/pypeit/PypeIt/pull/2201#pullrequestreview-5358481186
- Before posting, confirmed the PR head was still `ac9a085`, so the line numbers in the review are still correct. Replaced the internal preamble with a line saying which commit the line numbers refer to.

### 2026-09-29 (Inline comments on PR #2201)

- Posted a second "Comment" review as `profxj` with 20 inline comments, pinned to `ac9a085`: https://github.com/pypeit/PypeIt/pull/2201#pullrequestreview-5358522778. Each comment is tagged with its item in the summary review (Required 1–4, Suggestions 1–5, 7, 9, and the nits); Required 1 has a one-click `suggestion` block.
- Suggestion 6 (the leftover `VPH-ALL_7100` template, which is not in the diff) and Suggestion 8 (backward compatibility) have no single line to attach to, so they stay in the summary review only.
- Learned: GitHub only accepts inline comments on RIGHT-side lines inside diff hunks. The payload builder (scratchpad `build_inline.py`) checks each anchor against the diff before posting.
