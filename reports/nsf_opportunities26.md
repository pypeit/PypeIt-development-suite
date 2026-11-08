# NSF Funding Opportunities for PypeIt — FY2026

*Prepared: 2026-06-20. Source: searches of nsf.gov / new.nsf.gov and the live
NSF solicitation pages (cited inline). This is a standalone companion to
`claude_prompts/nsf_proposals.md`; it consolidates all Task-1 findings and
incorporates the user's Task-1 Q&A decisions.*

## Guiding decisions (from Q&A)

- **AI/ML is central.** The user directed that AI be a core thrust of any
  proposal — not a peripheral feature. Programs and framings below are therefore
  prioritized by how well they support a genuine AI research/development agenda.
- **Target the FY2026 deadline cycle** (not a later year).
- **Lead program and proposing team are still TBD** — both deliberately deferred
  by the user. This report keeps the top candidates open and flags eligibility
  constraints (e.g. PI status) that will matter once the team is chosen.

## What PypeIt is (for proposal framing)

PypeIt is an open-source, pure-Python pipeline that reduces raw spectroscopic
data from many astronomical spectrographs (long-slit, multi-slit, cross-dispersed
echelle, and slicer IFUs) into calibrated science-ready spectra. It is widely
used research software / community cyberinfrastructure for ground-based optical/
IR astronomy, with an active developer and user community and an extensive
real-data regression suite (the PypeIt development suite). This identity —
sustained open research software for a broad scientific community, with clear
opportunities to inject modern ML/AI — is what makes it competitive for the
programs below.

---

## Ranked opportunities

### Tier 1 — Strongest fit

#### 1. CSSI — Cyberinfrastructure for Sustained Scientific Innovation (NSF 22-632)

The single best-aligned program. CSSI exists precisely to build and **sustain**
open, community research software/cyberinfrastructure, and it explicitly welcomes
AI/ML-driven methods as a cross-cutting theme. With AI now central to the plan,
a CSSI proposal that embeds machine learning into PypeIt's reduction chain is
squarely in scope and plays to the program's priorities.

- **Award classes:**
  - *Elements* — small focused teams that create/deploy a needed service. The
    natural vehicle for a targeted "AI in PypeIt" development effort.
  - *Framework Implementations* — larger, interdisciplinary teams solving common
    research problems; appropriate if the AI agenda is broad and multi-institution.
  - *Transition to Sustainability* — one-time award to put an existing,
    widely-used tool on durable footing; strong given PypeIt's established user base.
- **Program-wide funding:** up to $10M total (Elements), $20M (Frameworks), $4M
  (Transition to Sustainability); project support up to 10 years.
- **FY2026 deadline:** full-proposal window **Oct 1 – Nov 16, 2026**.
- **Solicitation:** https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation/nsf22-632/solicitation
- **Program page:** https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation
- **Contact:** CSSIQueries@nsf.gov
- **Why it fits an AI-central plan:** funds *both* the AI R&D and the engineering
  needed to deliver it as sustained community infrastructure — the rare program
  that rewards "research-grade ML" *and* "production-grade open software" together.

#### 2. NSF Trailblazer Engineering Impact Award (NSF 26-502)

Elevated in this report because AI is a **named priority area** and the user wants
AI to be central. Supports a single PI pursuing bold, novel research with real-
world impact.

- **Funding:** up to **$3M over 3 years** (≈$15M total FY26 budget).
- **Best framing:** PypeIt as the platform for a *novel AI research thrust*
  (e.g. learned wavelength solutions, ML calibration/extraction), not maintenance.
- **Caveat:** single-PI, research-forward; less suited to broad team/community
  infrastructure work than CSSI. Check eligibility/career-stage constraints in the
  solicitation once the lead PI is chosen.
- **Solicitation:** https://www.nsf.gov/funding/opportunities/trailblazer-nsf-trailblazer-engineering-impact-award/nsf26-502/solicitation

#### 3. PESOSE — Pathways to Enable Secure Open-Source Ecosystems (NSF 26-506)

The current successor to POSE (replaces NSF 24-606). Funds turning an open-source
**product** into a sustainable open-source **ecosystem** — governance, community,
contributor base, and security. Directly relevant to PypeIt's open-source nature
and growing community. Note: this funds *ecosystem-building*, **not** new science
or AI features per se — so it is a **companion** to (not a substitute for) an
AI-central CSSI/Trailblazer proposal.

- **Tracks:**
  - Track 1 *Scoping & Planning* — ≤ $300K, ≤ 1 yr, ~30 awards.
  - Track 2 *Establishing & Expanding* — ≤ $1.5M, ≤ 2 yr, ~10 awards.
  - Track 3 *Improving Safety, Security & Privacy* — ≤ $1.5M, ≤ 2 yr, ~10 awards.
- **FY2026 deadlines:** **Sep 1, 2026** and **Mar 2, 2027** (annual, first Tuesdays).
- **Eligibility note:** PIs at universities must hold tenured/tenure-track (or
  equivalent research) positions — relevant once the team is set.
- **Solicitation:** https://www.nsf.gov/funding/opportunities/pesose-pathways-enable-secure-open-source-ecosystems/nsf26-506/solicitation
- **POSE FAQ (still informative):** https://www.nsf.gov/pubs/2024/nsf24119/nsf24119.jsp

### Tier 2 — Good fit (astronomy-specific)

#### 4. AAG — Astronomy and Astrophysics Research Grants (NSF 22-624)

The MPS/AST workhorse grant. Explicitly "considers proposals for projects and
tools that *enable or enhance* astronomical research," so PypeIt development tied
to a clear science driver is eligible. With AI central, an AAG proposal would
pair an ML methods development effort with a concrete science program it unlocks.

- **Program page:** https://www.nsf.gov/funding/opportunities/aag-astronomy-astrophysics-research-grants
- **Best when:** the software/AI work is motivated by and bundled with a science
  program rather than pure infrastructure.

#### 5. ATI — Advanced Technologies and Instrumentation for the Astronomical Sciences (NSF 22-627)

Funds "hardware and/or software development and/or analysis to enable new types of
astronomical observations." A strong fit if PypeIt work is framed as *enabling
observations otherwise difficult/impossible* (new instrument support, advanced
IFU/echelle reduction, ML-enabled extraction).

- **FY2026 deadline:** full proposals **Oct 1 – Nov 15/16, 2026**.
- **Program page:** https://www.nsf.gov/funding/opportunities/ati-advanced-technologies-instrumentation-astronomical-sciences
- **Contact:** Matthew A. Bershady (mbershad@nsf.gov)

### Tier 3 — AI-emphasis programs (resource or out-of-scope)

#### 6. NAIRR — National AI Research Resource / NAIRR-OC (NSF 25-546)

Provides **compute/data access** for AI research rather than software-development
funding. Most useful as a **resource** to support an AI component of PypeIt (e.g.
training/serving models), not a direct funding target. Worth citing in an
AI-central proposal as the compute pathway.

- https://www.nsf.gov/funding/opportunities/nairr-oc-foundations-operating-national-artificial-intelligence/nsf25-546/solicitation

#### 7. CAIG — Collaborations in AI and Geosciences (NSF 25-530)

AI for Earth-system science; **astronomy is out of scope**. Listed only to record
that, despite its attractive AI framing, it is **not** a fit.

- https://www.nsf.gov/funding/opportunities/caig-collaborations-artificial-intelligence-geosciences/nsf25-530/solicitation

---

## AI/ML thrusts to anchor an AI-central proposal

Concrete machine-learning directions inside PypeIt that can form the research core
of a CSSI / Trailblazer / AAG proposal. Each maps to an existing, well-defined
stage of the reduction chain, so the AI work is grounded in real software and
real regression data:

- **Learned cosmic-ray / artifact rejection** and bad-pixel identification.
- **Automated wavelength-solution identification** — replacing/augmenting the
  current `pypeit_identify` + arc-template workflow with learned line matching and
  solution selection. (PypeIt's existing arc-template archive and the `arclines`
  data are a ready-made training corpus.)
- **Learned sky-subtraction and scattered-light modeling.**
- **Automated frame typing / setup classification** from raw FITS metadata and
  pixel data.
- **Anomaly detection / data-quality flagging** across the reduction chain, with
  uncertainty quantification.
- **An "agentic" reduction-diagnosis assistant** for non-expert users (building on
  the existing `diagnose-reduction` skill/workflow concept).

These also create a natural **sustainability + community** story (CSSI/PESOSE):
labeled training data from the dev suite, reproducible model artifacts, and
community contribution pathways for new instruments.

---

## Recommendation (given AI-central + FY2026)

1. **Primary: CSSI (Nov 16, 2026 window).** Best home for an AI-central effort
   that must also ship as sustained community software. Choose *Elements* for a
   focused AI thrust, or *Transition to Sustainability* if the emphasis is
   durability of the existing, widely-used pipeline.
2. **Strong alternative / parallel: Trailblazer (NSF 26-502)** if the team wants a
   single-PI, research-forward AI vehicle with a $3M ceiling — AI is a named
   priority there.
3. **Companion: PESOSE Track 1 (Sep 1, 2026 / Mar 2, 2027).** Low-cost planning
   award to formalize governance/sustainability; strengthens the CSSI
   sustainability narrative.
4. **Astronomy-science fallback: AAG / ATI (~Nov 2026)** when the work is tied to
   a specific science or instrumentation driver.

### Open items to resolve before drafting

- **Pick the lead program** (CSSI vs. Trailblazer as the AI-central anchor) —
  deferred by the user.
- **Confirm the proposing team / lead institution and PI eligibility** — deferred;
  note PESOSE's tenured/tenure-track PI requirement and any Trailblazer
  career-stage limits.

## Sources

- CSSI program page — https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation
- CSSI solicitation (NSF 22-632) — https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation/nsf22-632/solicitation
- PESOSE solicitation (NSF 26-506) — https://www.nsf.gov/funding/opportunities/pesose-pathways-enable-secure-open-source-ecosystems/nsf26-506/solicitation
- POSE FAQ (NSF 24-119) — https://www.nsf.gov/pubs/2024/nsf24119/nsf24119.jsp
- AAG program page — https://www.nsf.gov/funding/opportunities/aag-astronomy-astrophysics-research-grants
- ATI program page — https://www.nsf.gov/funding/opportunities/ati-advanced-technologies-instrumentation-astronomical-sciences
- Trailblazer Engineering Impact Award (NSF 26-502) — https://www.nsf.gov/funding/opportunities/trailblazer-nsf-trailblazer-engineering-impact-award/nsf26-502/solicitation
- NAIRR-OC (NSF 25-546) — https://www.nsf.gov/funding/opportunities/nairr-oc-foundations-operating-national-artificial-intelligence/nsf25-546/solicitation
- CAIG (NSF 25-530) — https://www.nsf.gov/funding/opportunities/caig-collaborations-artificial-intelligence-geosciences/nsf25-530/solicitation
