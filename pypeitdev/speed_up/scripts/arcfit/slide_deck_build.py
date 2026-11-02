"""
Build pypeitdev/speed_up/Slides/Speed_Up_Arcfit_Slide_Deck.pptx (C.md Task 8).

Every number comes from Reports/Speed_Up_Arcfit_Report.md, the "After
workstream C" section of Reports/speed_up_results.md, or C.md (Q&A, Logs).
Figures come from slide_fig_*.py (drawn at placed size, fonts >= 20 pt).

Needs python-pptx (not in pypeit14b); e.g.
  pip install --target <scratch>/pylib python-pptx
  PYTHONPATH=<scratch>/pylib python slide_deck_build.py
The script ends with the font-size check: every run must be >= 20 pt.
"""
import os
import re
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

DEV = os.environ.get('PYPEIT_DEV', '/mnt/tank/Astronomy/PypeIt/PypeIt-development-suite').rstrip('/')
SLIDES = f'{DEV}/pypeitdev/speed_up/Slides'
FIG = f'{SLIDES}/figures'
OUT = f'{SLIDES}/Speed_Up_Arcfit_Slide_Deck.pptx'

INK = RGBColor(0x1B, 0x24, 0x30)
MUTED = RGBColor(0x5F, 0x63, 0x6B)
BLUE = RGBColor(0x2A, 0x78, 0xD6)
DBLUE = RGBColor(0x1B, 0x4F, 0x91)
ORANGE = RGBColor(0xEB, 0x68, 0x34)
BAND = RGBColor(0xE6, 0xEF, 0xFA)     # summary-line highlight
CARD = RGBColor(0xF3, 0xF4, 0xF6)     # card tint
WARM = RGBColor(0xFD, 0xEE, 0xE7)     # option-B / caution tint
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
DARK = RGBColor(0x1B, 0x24, 0x30)
FONT = 'Calibri'
BODY = 20           # minimum and body size
L, W = 0.5, 12.33   # content left edge and width (in)
TOP = 2.2           # content top (below title + summary)
BOT = 6.85          # content bottom (slide number below)

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
nslide = [0]


def _run(p, text, size=BODY, bold=False, color=INK, italic=False, font=FONT):
    r = p.add_run()
    # keep digit groups ("12 394") and number-unit pairs together on one line
    text = re.sub(r'(?<=\d) (?=\d{3}\b)', '\u00a0', text)
    r.text = text
    f = r.font
    f.size, f.bold, f.italic, f.name = Pt(size), bold, italic, font
    f.color.rgb = color
    return r


def text(slide, x, y, w, h, paras, size=BODY, color=INK, align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.TOP, space=6, name=None):
    """paras: list of str or list of (text, dict) runs, or dicts {runs, bullet, size}."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for m in ('margin_left', 'margin_right', 'margin_top', 'margin_bottom'):
        setattr(tf, m, Inches(0.05))
    first = True
    for para in paras:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = align
        p.space_after = Pt(space)
        if isinstance(para, str):
            para = [(para, {})]
        bullet = False
        if isinstance(para, dict):
            bullet = para.get('bullet', False)
            para = para['runs']
        if bullet:
            _run(p, '•  ' if bullet is True else f'{bullet}  ', size=size, bold=True,
                 color=ORANGE if bullet is True else BLUE)
        for t, o in para:
            _run(p, t, size=o.get('size', size), bold=o.get('bold', False),
                 color=o.get('color', color), italic=o.get('italic', False),
                 font=o.get('font', FONT))
    return tb


def box(slide, x, y, w, h, fill, rounded=True, name=None):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    shp.line.fill.background()
    shp.shadow.inherit = False
    if rounded:
        shp.adjustments[0] = 0.06
    if name:
        shp.name = name
    return shp


def new_slide(title, summary, dark=False):
    s = prs.slides.add_slide(BLANK)
    nslide[0] += 1
    bg = s.background.fill
    bg.solid()
    bg.fore_color.rgb = DARK if dark else WHITE
    text(s, L, 0.32, W, 0.75, [[(title, {'bold': True})]], size=34,
         color=WHITE if dark else INK, anchor=MSO_ANCHOR.MIDDLE, name='Title')
    if summary:
        box(s, L, 1.15, W, 0.88, RGBColor(0x2B, 0x38, 0x4A) if dark else BAND, name='Summary band')
        text(s, L + 0.2, 1.15, W - 0.4, 0.88, [[(summary, {'bold': False})]],
             color=WHITE if dark else INK, anchor=MSO_ANCHOR.MIDDLE, space=0, name='Summary')
    text(s, 12.0, 6.98, 0.83, 0.4, [str(nslide[0])], color=RGBColor(0xA0, 0xA4, 0xAB) if dark else MUTED,
         align=PP_ALIGN.RIGHT, space=0, name='Slide number')
    return s


def card(s, x, y, w, h, head, body, fill=CARD, head_color=INK, size=BODY, space=6):
    box(s, x, y, w, h, fill)
    paras = [[(head, {'bold': True, 'size': 22, 'color': head_color})]] + body
    text(s, x + 0.15, y + 0.12, w - 0.3, h - 0.24, paras, size=size, space=space)


def stat(s, x, y, w, h, big, label, fill=CARD, big_color=BLUE, big_size=40):
    box(s, x, y, w, h, fill)
    text(s, x + 0.15, y + 0.1, w - 0.3, h - 0.2,
         [[(big, {'bold': True, 'size': big_size, 'color': big_color})], label], space=2)


def table(s, x, y, w, rows, colw, rowh=0.46, header_fill=INK, bold_rows=(), first_bold=False,
          align_num=True):
    nr, nc = len(rows), len(rows[0])
    gt = s.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(rowh * nr))
    t = gt.table
    for j, cw in enumerate(colw):
        t.columns[j].width = Inches(cw)
    for i in range(nr):
        t.rows[i].height = Inches(rowh)
        for j in range(nc):
            c = t.cell(i, j)
            c.margin_left = c.margin_right = Inches(0.08)
            c.margin_top = c.margin_bottom = Inches(0.03)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid()
            c.fill.fore_color.rgb = header_fill if i == 0 else (CARD if i % 2 == 0 else WHITE)
            tf = c.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if (j == 0 or not align_num) else PP_ALIGN.RIGHT
            val = rows[i][j]
            bold = i == 0 or i in bold_rows or (first_bold and j == 0)
            color = WHITE if i == 0 else INK
            if isinstance(val, tuple):
                val, o = val
                bold = o.get('bold', bold)
                color = o.get('color', color)
            _run(p, val, size=BODY, bold=bold, color=color)
    # remove the default table style banding look
    tblPr = gt._element.graphic.graphicData.tbl.tblPr
    tblPr.set('bandRow', '0')
    tblPr.set('firstRow', '0')
    return gt


def picture(s, path, x, y, w):
    return s.shapes.add_picture(path, Inches(x), Inches(y), width=Inches(w))


def notes(s, txt):
    s.notes_slide.notes_text_frame.text = txt


# ---------------------------------------------------------------- 1 title
s = prs.slides.add_slide(BLANK)
nslide[0] += 1
s.background.fill.solid()
s.background.fill.fore_color.rgb = DARK
text(s, 0.8, 1.2, 11.7, 0.5, [[('PYPEIT SPEED-UP  ·  WORKSTREAM C', {'bold': True, 'color': ORANGE})]],
     size=22)
text(s, 0.8, 1.8, 11.7, 1.9, [[('Vectorized arc-line Gaussian fitting', {'bold': True})]],
     size=48, color=WHITE)
text(s, 0.8, 3.15, 11.7, 0.7, [[('+ the archive-sky FWHM cache in spectral flexure', {})]],
     size=28, color=RGBColor(0xC9, 0xD6, 0xE8))
text(s, 0.8, 4.6, 11.7, 1.6,
     [[('PypeIt PR C  ·  branch ', {}), ('speed_up_arcfit', {'font': 'Courier New'}),
       (' (on ', {}), ('speed_up_B', {'font': 'Courier New'}), (', PR #2209)', {})],
      [('Dev-suite branch ', {}), ('dev_speed_up', {'font': 'Courier New'})],
      [('2026-10-28', {})]],
     size=22, color=WHITE, space=8)
notes(s, 'Companion deck to the PR C description. Source: Reports/Speed_Up_Arcfit_Report.md.')

# ---------------------------------------------------------------- 2 replaced
s = new_slide('What was replaced, and why',
              'Fitting each arc line with its own curve_fit call cost ~840 s of a 12 394 s serial '
              'DEIMOS 600ZD reduction, almost all of it per-call overhead.')
cw, gap = (W - 2 * 0.3) / 3, 0.3
for k, (big, lab) in enumerate([('1.11 M', 'curve_fit calls (fit_gauss)'),
                                ('840 s', 'cumulative in fit_gauss'),
                                ('93 M', 'gauss_3deg calls (332 s self)')]):
    stat(s, L + k * (cw + gap), TOP, cw, 1.45, big, lab, big_size=44)
text(s, L, TOP + 1.75, W, 2.8, [
    {'bullet': True, 'runs': [('June 2026 serial profile of Keck/DEIMOS 600ZD_M_6500: 12 393.8 s wall.', {})]},
    {'bullet': True, 'runs': [('fit_arcspec', {'bold': True}), (' sits behind detect_lines, which feeds wavelength '
                              'calibration, tilts, object finding, slit edges and flexure.', {})]},
    {'bullet': True, 'runs': [('Each fit is a few pixels and 3 parameters, so the cost is per-call overhead: '
                              'a natural target for vectorization.', {})]},
    {'bullet': True, 'runs': [('Re-profiled on the current code with the legacy path: 1 114 679 calls, 846 s. '
                              'The baseline reproduces.', {})]},
], space=10)

# ---------------------------------------------------------------- 3 method
s = new_slide('The new fit: one batched cascade for all lines',
              'fit_arcspec now fits every detected line at once in four steps, keeping its signature '
              'and the −999 sentinel for failed fits.')
steps = [('1  Weighted log-parabola',
          ['ln y = a₀ + a₁u + a₂u² per line', 'Weights w² = y² match '
           'curve_fit\'s objective to first order', 'One stacked 3×3 solve']),
         ('2  Gauss–Newton polish',
          ['Levenberg steps on curve_fit\'s own objective, all pixels', 'Cap ARCFIT_GN_ITER = 10',
           'Per-line converged flag']),
         ('3  Moment-seeded rescue',
          ['For lines the parabola cannot seed', 'guess_gauss-style seed', 'Up to 20 iterations']),
         ('4  curve_fit fallback',
          ['Legacy loop kept verbatim', 'Last resort and reference', 'Hardened: no NaN leaks'])]
cw, gap = (W - 3 * 0.3) / 4, 0.3
for k, (head, body) in enumerate(steps):
    x = L + k * (cw + gap)
    card(s, x, TOP, cw, 4.45, head, [{'bullet': True, 'runs': [(b, {})]} for b in body],
         fill=WARM if k == 3 else CARD, head_color=ORANGE if k == 3 else DBLUE, space=10)
notes(s, 'Report section 2. Rescue added in Task 3 (Q4b); without it 6.7% of DEIMOS and 13.2% of sky '
         'lines fell back to curve_fit.')

# ---------------------------------------------------------------- 4 A vs B
s = new_slide('Shipped default A, and option B',
              'Option A hands every line the fast fit is unsure of back to curve_fit, so it accepts '
              'exactly the legacy lines; option B never does.')
hw = (W - 0.3) / 2
card(s, L, TOP, hw, 3.7, 'A  (shipped):  ARCFIT_MATCH_CURVEFIT = True', [
    'Hand a line to curve_fit when:',
    {'bullet': '1', 'runs': [('the polish has not converged in GN_ITER', {})]},
    {'bullet': '2', 'runs': [('the rescue has not converged in RESCUE_ITER', {})]},
    {'bullet': '3', 'runs': [('σ < 0.35 × window half-width (spikes)', {})]},
    [('Same accepted lines as legacy on all 370 dev-suite arc spectra.', {'bold': True})],
], head_color=DBLUE, space=8)
card(s, L + hw + 0.3, TOP, hw, 3.7, 'B  (option):  no hand-offs', [
    {'bullet': True, 'runs': [('fast_line_fit = True', {'bold': True}),
                              (' in [calibrations][wavelengths]: wavelength calibration only', {})]},
    {'bullet': True, 'runs': [('ARCFIT_MATCH_CURVEFIT = False', {'bold': True}), (': all callers', {})]},
    {'bullet': True, 'runs': [('Fit 13.3× faster than legacy on 92 DEIMOS arcs', {})]},
    {'bullet': True, 'runs': [('Changes which lines are accepted (next slides)', {})]},
], fill=WARM, head_color=ORANGE, space=8)

# ---------------------------------------------------------------- 5 fidelity
s = new_slide('Fidelity to legacy curve_fit',
              'With the polish capped at 10 iterations the vectorized fit tracks curve_fit to 4.4e-4 px '
              '(p99) on DEIMOS arcs, and under A the accepted lines are identical.')
text(s, L, TOP - 0.05, W, 0.45, [[('DEIMOS 600ZD arcs: 92 spectra, 10 170 lines (vectorized fit alone, as in B)',
                                    {'bold': True})]])
table(s, L, TOP + 0.45, W, [
    ['GN_ITER', 'fit_arcspec', 'speed-up', 'lines lost', 'median / p99 |Δc|'],
    ['0 (seed only, no rescue)', '—', '—', '—', '0.009 / 0.25 px'],
    ['3', '0.29 s', '16.1×', '122', '1.4e-6 / 0.14 px'],
    [('10 (shipped cap)', {'bold': True}), ('0.35 s', {'bold': True}), ('13.3×', {'bold': True}),
     ('51', {'bold': True}), ('1.4e-7 / 4.4e-4 px', {'bold': True})],
], colw=[3.9, 2.0, 1.7, 1.8, 2.93], rowh=0.5)
cw = (W - 2 * 0.3) / 3
y2 = TOP + 2.65
stat(s, L, y2, cw, 2.0, '0 / 0', 'lines gained / lost under A: 370 spectra, 33 354 legacy lines', big_size=36)
stat(s, L + cw + 0.3, y2, cw, 2.0, '< 1e-3 px', 'agreement of shared centroids under A', big_size=36)
stat(s, L + 2 * (cw + 0.3), y2, cw, 2.0, '0.87 → 0.33 px', 'Paranal sky p99 |Δc|, GN 3 → 10',
     big_size=36)
notes(s, 'Report section 3. Seed-only on the Paranal sky: median 0.040 px, p99 2.0 px, 13.2% fall back.')

# ---------------------------------------------------------------- 6 line-count changes
s = new_slide('Line-count changes under option B',
              'Over 370 spectra option B gains 83 and loses 56 of the legacy lines, for three reasons; '
              'centroids of shared lines are unaffected.')
cw = (W - 2 * 0.3) / 3
card(s, L, TOP, cw, 3.3, 'Negative σ (Q3)', [
    'Paranal sky: 56 legacy fits end at σ < 0 and are rejected; 49 pass every other cut.',
    [('B keeps 38: ', {'bold': True}), ('37 are 0.25–0.66 px spikes.', {})],
], space=10)
card(s, L + cw + 0.3, TOP, cw, 3.3, 'Saturated flat tops (Q6)', [
    'DEIMOS arcs: B rejects 51 legacy lines (129 at GN 0).',
    'All fail the width cut; 98–99% have raw peaks > 30 000 counts.',
], fill=WARM, head_color=ORANGE, space=10)
card(s, L + 2 * (cw + 0.3), TOP, cw, 3.3, 'Sub-resolution spikes', [
    [('B gains 83 lines ', {'bold': True}), ('with σ < 0.5× the expected width (median 0.33).', {})],
    'Legacy curve_fit raised RuntimeError on them.',
], space=10)

# ---------------------------------------------------------------- 7 saturated = plateau edges
s = new_slide('The "saturated lines" are plateau edges',
              'Of the 46 DEIMOS 600ZD features B rejects, 39 sit on the vertical edge of a ~39 900-count '
              'plateau and 6 on its flat top, so the legacy σ is an artefact.')
picture(s, f'{FIG}/slide_sat_examples.png', L, TOP - 0.05, W)
text(s, L, TOP + 4.3, 11.4, 0.45, [[('None of the 46 enters a final wavelength fit; under A all three '
                                     'return the legacy (blue) fit.', {'bold': True})]])
notes(s, 'Report 5.1. Pipeline inputs: 86 slits, fwhm 7.0-7.7 px, 11-px window, width cut sigma < 9.3 px. '
         'Nearest pixel_fit line >= 7.2 px away (median 26 px). Plateaus: ~60 px at ~39 900 counts around '
         'cores hard-saturated at 79 503 counts, consistent with charge bleeding (not verified).')

# ---------------------------------------------------------------- 8 negative sigma
s = new_slide('Negative σ: a sign, not a shape',
              'gauss_3deg uses only σ², so a negative legacy σ is an arbitrary sign; '
              'B keeps such lines and A rejects them exactly as legacy does.')
picture(s, f'{FIG}/slide_negsigma_example.png', L, TOP - 0.05, W)
text(s, L, TOP + 4.3, 11.4, 0.45, [[('Most are sub-resolution spikes: expected σ ≈ 1.7 px for '
                                     'fwhm = 4.', {'bold': True})]])
notes(s, 'Report 5.2. Left: the one line-like case (centres differ by 0.012 px, widths by 12%). '
         'Right: the typical case, same optimum, opposite sign.')

# ---------------------------------------------------------------- 9 why A ships
s = new_slide('Why option A ships (Q8)',
              'B passed every vet test, but losing one or two saturated features flipped wavelength '
              'identifications on several DEIMOS slits.')
stat(s, L, TOP, 4.2, 3.05, '~30 px', 'blue-end misidentification on 600ZD MSC01 slit 2 under B; '
     'slit RMS moves only 0.588 → 0.617', fill=WARM, big_color=ORANGE, big_size=54)
text(s, L + 4.55, TOP, W - 4.55, 3.2, [
    {'bullet': True, 'runs': [('10 of 38 WaveCalib files move > 0.05 px (multi-slit DEIMOS, LRIS multi_300)', {})]},
    {'bullet': True, 'runs': [('Flips also on MSC01 slit 11 and MSC03 slits 0, 14, 15', {})]},
    {'bullet': True, 'runs': [('830G_L_8100: RMS 0.051 → 0.27 px', {})]},
    {'bullet': True, 'runs': [('One improvement: 1200G_M_7750 max RMS 0.42 → 0.13', {})]},
], space=10)
box(s, L, TOP + 3.4, W, 1.2, BAND)
text(s, L + 0.2, TOP + 3.4, W - 0.4, 1.2, [
    [('Decision: ship A, keep B as an option. ', {'bold': True}),
     ('A wrong identification is worse than a smaller speed-up; after the cache B saves only ~13 s.', {})]],
     anchor=MSO_ANCHOR.MIDDLE, space=0)

# ---------------------------------------------------------------- 10 discussion
s = new_slide('Discussion: saturated features buy identifications',
              'Saturated features never enter a final fit, but they anchor the full_template pattern '
              'match, so keeping them keeps the identifications right.')
cw, ch = (W - 0.3) / 2, 2.15
items = [('Anchors, not fit lines',
          'Synthetic arcs from all peaks are cross-correlated with the template; bright, reproducible '
          'plateau edges dominate.'),
         ('Dropping them flipped IDs',
          'On B\'s misidentified 600ZD slits the only input difference was 1–2 missing saturated '
          'features.'),
         ('Stable edge markers',
          'The legacy centroid is not a line wavelength, but it moves with the edge: good positional '
          'information.'),
         ('Future, outside PR C',
          'Detect flat-topped windows, model the edge, flag as "saturated": an ID anchor excluded from '
          'the fit.')]
for k, (h, b) in enumerate(items):
    x = L + (k % 2) * (cw + 0.3)
    y = TOP + (k // 2) * (ch + 0.3)
    card(s, x, y, cw, ch, f'{k + 1}  {h}', [b], fill=WARM if k == 3 else CARD,
         head_color=ORANGE if k == 3 else DBLUE)
notes(s, 'Report section 6. The plateau model would also remove the slowest curve_fit calls.')

# ---------------------------------------------------------------- 11 fallback + cache
s = new_slide('Fallback rate, and where the time went',
              'Few lines reach curve_fit, but the largest saving came from caching the archive-sky FWHM '
              'that flexure re-measured 170 times per DEIMOS run.')
sw = 4.3
stat(s, L, TOP, sw, 1.45, '0.05%', 'of DEIMOS lines reach curve_fit under B', big_size=30)
stat(s, L, TOP + 1.6, sw, 1.45, '~5%', 'of DEIMOS lines handed off under A', big_size=30)
stat(s, L, TOP + 3.2, sw, 1.45, '432 of 517', 'A hand-offs fail in curve_fit, as in legacy',
     big_size=30, fill=WARM, big_color=ORANGE)
rx = L + sw + 0.35
rw = W - sw - 0.35
text(s, rx, TOP - 0.05, rw, 0.45, [[('Serial DEIMOS, option A: before / after the cache', {'bold': True})]])
table(s, rx, TOP + 0.45, rw, [
    ['', 'before', 'after'],
    ['measure_fwhm calls / time', '170 / 571 s', ('1 / 3 s', {'bold': True})],
    ['fit_gauss calls / time', '319 k / 494 s', ('14 k / 42 s', {'bold': True})],
    ['fit_arcspec time', '526 s', ('57 s', {'bold': True})],
], colw=[3.6, 2.2, rw - 5.8], rowh=0.5)
text(s, rx, TOP + 2.65, rw, 2.0, [
    {'bullet': True, 'runs': [('lru_cache on ', {}), ('_archive_sky_fwhm(sky_file)', {'bold': True}),
                              (': bit-identical, one cache per --ncpu worker', {})]},
    {'bullet': True, 'runs': [('maxfev cannot be capped: successful fits reach 799 evaluations', {})]},
], space=8)
notes(s, 'Report 7.1-7.2 (Q10). Sky: 36 of 4 259 lines reach curve_fit under B (0.8%). 22 plausible lines '
         'need > 400 evaluations. measure_fwhm was ~3.4 s per call. test_archive_fwhm_cache.')

# ---------------------------------------------------------------- 12 vetted RMS
s = new_slide('Dev-suite vetting: the vetted RMS',
              'All variants reduced 35/35 setups; vet outcomes are identical for legacy, A and gn10, '
              'and only gn0 fails (Kast red 600_7500_d57).')
F = {'bold': True, 'color': ORANGE}
table(s, L, TOP, W, [
    ['setup (vetted quantity)', 'tol', 'legacy', 'A', 'gn10 (B)', 'gn0'],
    ['Kast red 300_7500_Ne', '0.05', '0.0308', '0.0308', '0.0308', '0.0413'],
    ['Kast red 600_7500_d57', '0.05', '0.0369', '0.0369', '0.0369', ('0.0633 fail', F)],
    ['Kast red 1200_5000_d57', '0.055', '0.0512', '0.0512', '0.0512', '0.0496'],
    ['DEIMOS 1200B_LVM_5200', '0.15', '0.1234', '= legacy', '0.1234', '0.1242'],
    ['DEIMOS 600ZD_M_6500', '0.35', '0.3356', '= legacy', '0.3350', '0.3100'],
    ['DEIMOS 900ZD_LVM_5500', '0.1', '0.0687', '= legacy', '0.0687', '0.0685'],
    ['LRIS-b multi_300_5000_d680', '0.32', '0.2970', '0.2970', '0.2969', '0.2877'],
    ['LRIS-b long_400_3400_d560', '0.21', '0.1947', '0.1947', '0.1947', '0.1853'],
    ['LRIS-b long_600_4000_d560', '0.08', '0.0621', '0.0621', '0.0621', '0.0772'],
], colw=[4.33, 1.2, 1.6, 1.6, 1.8, 1.8], rowh=0.45)
notes(s, 'Report section 8; values from the Task 5 log and REDUX_OUT_C/arcfit_metrics*.md. gn10 = B; '
         'gn0 = analytic seed + rescue, no hand-offs. Kast red 300_7500_Ne slit 1, DEIMOS rows a single slit, '
         'LRIS rows the max over slits.')

# ---------------------------------------------------------------- 13 WaveCalib by family
s = new_slide('Dev-suite vetting: WaveCalib by setup family',
              'A reproduces legacy WaveCalib (and every Tilts file) except one 600ZD slit; gn10 and gn0 '
              'move DEIMOS solutions by tens of pixels.')
table(s, L, TOP, W, [
    ['family (files)', 'A', 'gn10 (B)', 'gn0'],
    ['Kast blue/red (8)', 'identical', 'identical', ('600_7500_d57 fails', F)],
    ['LRIS blue (13)', 'identical', 'Δλ 0.096 px', 'Δλ ≤ 4 px'],
    ['other DEIMOS (10)', 'identical', ('830G_L_8100 0.27', F), 'Δλ ≤ 20 px'],
    ['DEIMOS 1200B/G_M (3)', 'identical', 'Δλ ≤ 28.5 px', 'Δλ ≤ 104 px'],
    ['600ZD MSC01–04 (4)', ('1 slit: 30.7 px', {'bold': True}), '31–52 px', '17–49 px'],
    ['lines in final fits', '21 710', '21 681', '21 700'],
], colw=[3.63, 2.9, 2.9, 2.9], rowh=0.5, align_num=False)
text(s, L, TOP + 3.75, W, 0.9, [
    [('Legacy: 21 702 lines in final fits. ', {'bold': True}),
     ('A\'s one change (MSC03 slit 11) is the known limitation, slide 16.', {})]])
notes(s, 'Report section 8 family table. gn10: 830G_L_8100 RMS 0.051 -> 0.27; 1200G_M_7750 max RMS '
         '0.42 -> 0.13. gn0: Kast 600_7500_d57 RMS 0.037 -> 0.063. Tilts: A median |dtilt| <= 1.3e-7 px, '
         'max <= 0.005 px; gn10 up to 0.88 px; gn0 up to 62 px (1200G_M_5500).')

# ---------------------------------------------------------------- 14 re-profiling
s = new_slide('Re-profiling serial DEIMOS 600ZD',
              'On the same code, serial DEIMOS 600ZD drops from 10 420 s (legacy) to 9 462 s with '
              'A + cache (−9.2%), and line fitting now costs 57 s.')
picture(s, f'{FIG}/slide_profile_deimos.png', L, TOP - 0.05, W)
text(s, L, TOP + 4.3, 11.4, 0.45, [[('Five serial cProfile runs on the same HEAD; the cache saves more '
                                     'than the fit change.', {'bold': True})]])
notes(s, 'Report section 9 / speed_up_results.md "After workstream C". fit_gauss calls: legacy 1.11 M, '
         'A 319 k, A + cache 14 k.')

# ---------------------------------------------------------------- 15 totals
s = new_slide('A + B + C totals',
              'Across the three workstreams, serial DEIMOS 600ZD went from 12 394 s to 9 462 s and '
              'serial Kast blue from 121 s to 88 s.')
cw = (W - 3 * 0.3) / 4
for k, (big, lab, hl) in enumerate([
        ('−23.7%', 'serial DEIMOS: 12 394 → 9 462 s (A+B+C)', True),
        ('−9.2%', 'C alone, same code: 10 420 → 9 462 s', False),
        ('−4.3%', 'DEIMOS --ncpu 4: 3 400 → 3 255 s', False),
        ('−27%', 'serial Kast blue: 121 → 88 s (mostly A)', False)]):
    stat(s, L + k * (cw + 0.3), TOP, cw, 2.1, big, lab, big_size=44,
         big_color=DBLUE if hl else BLUE)
box(s, L, TOP + 2.4, W, 2.25, CARD)
text(s, L + 0.2, TOP + 2.5, W - 0.4, 2.1, [
    [('Caveats', {'bold': True, 'size': 22})],
    {'bullet': True, 'runs': [('Walls carry a few % machine noise; the per-function numbers are robust.', {})]},
    {'bullet': True, 'runs': [('June baseline: different code, day and env; A+B+C totals are indicative.', {})]},
    {'bullet': True, 'runs': [('At --ncpu 4 the fit change alone is within the ~3% scatter.', {})]},
], space=4)
notes(s, 'Report section 9. --ncpu 4: legacy 3 400.2 s, A + cache 3 254.8 s on the same code; vs the B '
         'result 3 391.2 s it is -4.0%. Kast blue 600_4310_d55 serial: 121.2 s (June) -> 88.4 s. The A + cache '
         'profile overlapped the unit tests and doc build for ~15 min.')

# ---------------------------------------------------------------- 16 known limitation
s = new_slide('Known limitation: 600ZD MSC03 slit 11',
              'A moves one DEIMOS slit from a misidentified solution to the correct one, through the '
              'fragile full_template shift-and-stretch step.')
picture(s, f'{FIG}/slide_msc03_slit11.png', L, TOP - 0.05, W)
text(s, L, TOP + 4.3, 11.4, 0.45, [[('A 1.6e-8 change in the cross-correlation flips the optimizer '
                                     '(Q11: documented only).', {'bold': True})]])
notes(s, 'Report section 10. Same lines in all 15 detect_lines calls (centres within 2e-6 to 8.5e-4 px). '
         'Shift/stretch legacy 75.86/0.98896 vs A 98.59/0.97535; solutions differ by 38 px at pixel 0. '
         'Blue lines within 0.5 px of the line list: A 10 of 12 (median 0.14 px), legacy 2 of 12 (median '
         '7.2 px). A: RMS 0.251 with 81 lines vs legacy 0.276 with 73. Slit 9 may be misidentified in both '
         'runs (not verified).')

# ---------------------------------------------------------------- 17 switches and risks
s = new_slide('Switches, options and remaining risks',
              'One flag restores the legacy fit exactly; the main risks are untested detect_lines callers '
              'and the fragile identification step.')
lw = 5.2
card(s, L, TOP, lw, 3.9, 'Switches', [
    [('USE_VECTORIZED_ARCFIT = False', {'bold': True}), (': kill switch, bit-identical legacy', {})],
    [('ARCFIT_MATCH_CURVEFIT = False', {'bold': True}), (': option B, all callers', {})],
    [('fast_line_fit = True', {'bold': True}), (': option B, wavelength calibration only', {})],
], head_color=DBLUE, space=12)
card(s, L + lw + 0.3, TOP, W - lw - 0.3, 3.9, 'Risks', [
    {'bullet': '1', 'runs': [('Other detect_lines callers: only 4 instruments reduced; watch echelle and IFU', {})]},
    {'bullet': '2', 'runs': [('Identification fragility; multi-slit DEIMOS not reproducible (#2208)', {})]},
    {'bullet': '3', 'runs': [('curve_fit cost on sky-line-heavy spectra under A', {})]},
    {'bullet': '4', 'runs': [('NaN fits now return −999 (a fix)', {})]},
], fill=WARM, head_color=ORANGE, space=8)
notes(s, 'Report section 11. Kill switch checked by test_fit_arcspec_kill_switch. centerr is unused '
         'numerically; revalidate if a caller starts weighting by it.')

# ---------------------------------------------------------------- 18 closing
s = new_slide('What ships, and what remains before merge',
              'PR C ships the legacy-faithful vectorized fit and the FWHM cache; it needs PR B merged, '
              'a retarget to develop and one full dev-suite pass.', dark=True)
hw = (W - 0.3) / 2
DC = RGBColor(0x2B, 0x38, 0x4A)
box(s, L, TOP, hw, 3.4, DC)
text(s, L + 0.2, TOP + 0.15, hw - 0.4, 4.3, [
    [('Ships in PR C', {'bold': True, 'size': 24, 'color': ORANGE})],
    {'bullet': True, 'runs': [('Vectorized fit_arcspec, default A', {})]},
    {'bullet': True, 'runs': [('Option B via fast_line_fit', {})]},
    {'bullet': True, 'runs': [('Archive-sky FWHM cache', {})]},
    {'bullet': True, 'runs': [('Kill switch and unit tests', {})]},
], color=WHITE, space=10)
box(s, L + hw + 0.3, TOP, hw, 3.4, DC)
text(s, L + hw + 0.5, TOP + 0.15, hw - 0.4, 4.3, [
    [('Before merge', {'bold': True, 'size': 24, 'color': ORANGE})],
    {'bullet': '1', 'runs': [('PR B (#2209): re-review and merge', {})]},
    {'bullet': '2', 'runs': [('Retarget PR C to develop and sync', {})]},
    {'bullet': '3', 'runs': [('One full ./pypeit_test all at the tip of speed_up_arcfit', {})]},
], color=WHITE, space=10)

prs.save(OUT)
print('wrote', OUT, len(prs.slides), 'slides')

# ---------------------------------------------------------------- font-size check
chk = Presentation(OUT)
sizes, unset = [], 0


def runs_of(shape):
    if shape.has_text_frame:
        for p in shape.text_frame.paragraphs:
            for r in p.runs:
                yield r
    if getattr(shape, 'has_table', False) and shape.has_table:
        for row in shape.table.rows:
            for c in row.cells:
                for p in c.text_frame.paragraphs:
                    for r in p.runs:
                        yield r


for sl in chk.slides:
    for sh in sl.shapes:
        for r in runs_of(sh):
            if not r.text.strip():
                continue
            if r.font.size is None:
                unset += 1
            else:
                sizes.append(r.font.size.pt)
assert unset == 0, f'{unset} runs without an explicit size'
assert min(sizes) >= 20, f'font smaller than 20 pt: {min(sizes)}'
print(f'font check: {len(sizes)} runs, minimum {min(sizes):.0f} pt, {unset} runs with no explicit size')
