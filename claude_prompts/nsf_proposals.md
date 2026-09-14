# NSF Proposals for PypeIt

## Goals

Seek out NSF funding opportunities for PypeIt.  And then write them.

## Claude

### Skills

Consider using the skills in PypeIt-development-suite/.claude/skills/

Examine the following to gain context:

- The PypeIt code base located in the PypeIt repository: `PypeIt/`
- The PypeIt-development-suite code base located in the PypeIt-development-suite repository: `PypeIt-development-suite/`
- The `arclines` repository located in the `arclines/` repository

## Prompts

1. Read this doc.  Execute the first task under Tasks.   
2. Read this doc.  Execute the 2nd task under Tasks.   

## Searching

### Tasks

1.  Search the NSF website for funding opportunities for PypeIt.  

    - Emphasize AI opportunities
    - Emphasize open source opportunities

Put your findings in the Report section.
If you have any questions, please ask them in the Q&A section. Log your work

2.  See my answers in the Q&A section.  Generate a separate report named reports/nsf_opportunities26.md with all of you findings.

### Report

*Last updated: 2026-06-20. PypeIt is an open-source, pure-Python pipeline that
reduces raw spectroscopic data from many astronomical instruments — i.e. it is
research software / cyberinfrastructure for the astronomy community. The
opportunities below are ranked by fit. Budget caps and deadlines are taken from
the NSF solicitation pages cited; confirm against the live solicitation before
writing, since several solicitations are at/near their re-issue cycle.*

#### Tier 1 — Strongest fit

**1. CSSI — Cyberinfrastructure for Sustained Scientific Innovation (NSF 22-632)**
The single best-aligned program. CSSI exists precisely to fund and *sustain*
open, community research software/cyberinfrastructure. PypeIt is a textbook
candidate.
- **Award classes:** *Elements* (small focused teams building/deploying a needed
  service — the natural fit for a feature/AI push in PypeIt), *Framework
  Implementations* (larger interdisciplinary teams), and *Transition to
  Sustainability* (one-time award to put an existing widely-used tool on a
  durable footing — also highly relevant given PypeIt's existing user base).
- **Program-wide funding:** up to $10M total for Elements, $20M for Frameworks,
  $4M for Transition to Sustainability; support up to 10 years.
- **Deadline:** full-proposal window **Oct 1 – Nov 16, 2026** (per current cycle).
- **AI angle:** NSF lists AI as a cross-cutting focus area; a CSSI proposal that
  embeds ML/AI methods into PypeIt (see "AI hooks" below) is squarely in scope.
- Page: https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation
- Solicitation: https://www.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation/nsf22-632/solicitation
- Contact: CSSIQueries@nsf.gov

**2. PESOSE — Pathways to Enable Secure Open-Source Ecosystems (NSF 26-506)**
The current successor to the POSE program (replaces NSF 24-606). Explicitly about
turning an open-source *product* into a sustainable open-source *ecosystem*
(governance, community, contributor base, security) — directly relevant to
PypeIt's open-source nature and growing developer/user community.
- **Tracks:** Track 1 *Scoping & Planning* (≤ $300K, ≤ 1 yr, ~30 awards);
  Track 2 *Establishing & Expanding* (≤ $1.5M, ≤ 2 yr, ~10 awards);
  Track 3 *Improving Safety, Security & Privacy* (≤ $1.5M, ≤ 2 yr, ~10 awards).
- **Deadlines:** **Sep 1, 2026** and **Mar 2, 2027** (annual, first Tuesdays).
- **Note:** funds *ecosystem-building* (governance, sustainability, community
  growth), **not** new science features per se. Track 1 is a low-cost way to
  formalize PypeIt governance/sustainability planning; Track 2 to scale the
  contributor ecosystem.
- Solicitation: https://www.nsf.gov/funding/opportunities/pesose-pathways-enable-secure-open-source-ecosystems/nsf26-506/solicitation
- FAQ (POSE): https://www.nsf.gov/pubs/2024/nsf24119/nsf24119.jsp

#### Tier 2 — Good fit (astronomy-specific)

**3. AAG — Astronomy and Astrophysics Research Grants (NSF 22-624)**
The MPS/AST workhorse grant. Explicitly "considers proposals for projects and
tools that *enable or enhance* astronomical research," so a PypeIt-development
proposal tied to a clear science driver is eligible. Best when the software work
is motivated by/bundled with a science program rather than pure infrastructure.
- Page: https://www.nsf.gov/funding/opportunities/aag-astronomy-astrophysics-research-grants

**4. ATI — Advanced Technologies and Instrumentation for the Astronomical Sciences (NSF 22-627)**
Funds "hardware and/or software development and/or analysis to enable new types
of astronomical observations." A strong fit if PypeIt work is framed as enabling
observations otherwise difficult/impossible (e.g. new instrument support, IFU/
echelle data-reduction techniques).
- **Deadline:** full proposals **Oct 1 – Nov 15/16** annually.
- Page: https://www.nsf.gov/funding/opportunities/ati-advanced-technologies-instrumentation-astronomical-sciences
- Contact: Matthew A. Bershady (mbershad@nsf.gov)

#### Tier 3 — AI-emphasis programs (fit depends on framing)

**5. NSF Trailblazer Engineering Impact Award (NSF 26-502)** — single-PI, bold
projects with AI as a named priority area; up to **$3M over 3 years** (~$15M FY26
budget). Would require framing PypeIt around a novel AI research thrust, not
maintenance. https://www.nsf.gov/funding/opportunities/trailblazer-nsf-trailblazer-engineering-impact-award/nsf26-502/solicitation

**6. NAIRR — National AI Research Resource / NAIRR-OC (NSF 25-546)** — provides
*compute/data access* for AI research rather than software-dev funding; useful as
a **resource** to support an AI component of PypeIt (e.g. training models), not a
direct funding target. https://www.nsf.gov/funding/opportunities/nairr-oc-foundations-operating-national-artificial-intelligence/nsf25-546/solicitation

**7. CAIG — Collaborations in AI and Geosciences (NSF 25-530)** — AI for Earth-
system science; **astronomy is out of scope**. Listed only to flag it is *not* a
fit despite the attractive AI framing.

#### "AI hooks" that would strengthen an AI-emphasis proposal
Concrete ML/AI thrusts that could be woven into a CSSI/Trailblazer/AAG proposal:
- ML-based cosmic-ray / artifact rejection and bad-pixel identification.
- Automated wavelength-solution identification (vs. the current
  `pypeit_identify` / arc-template workflow) using learned line-matching.
- Learned sky-subtraction and scattered-light modeling.
- Automated frame typing / setup classification from raw metadata + pixel data.
- Anomaly detection / data-quality flagging across the reduction chain.
- An "agentic" reduction-diagnosis assistant (cf. the existing
  `diagnose-reduction` workflow) for non-expert users.

#### Recommendation
Lead with **CSSI Elements or Transition-to-Sustainability (Nov 2026 window)** as
the primary target — it matches PypeIt's identity (sustained open research SW)
and accommodates an AI emphasis. Pair with a **PESOSE Track 1** planning proposal
(low cost, Sep 2026 / Mar 2027) to formalize ecosystem/governance/sustainability,
which also strengthens the CSSI sustainability narrative. Use **AAG/ATI** as the
fallback/complement when the work can be tied to a specific science or
instrumentation driver.

### Q&A

**Q1 (Claude → user):** Which opportunity should I focus on drafting first? My
recommendation is **CSSI** (primary) + **PESOSE Track 1** (companion), but ATI/
AAG are viable if you prefer an astronomy-science framing. Please pick a target.

We will let you know later.

**Q2 (Claude → user):** How central should AI/ML be? Options: (a) AI as a core
research thrust (enables Trailblazer + strengthens CSSI), (b) AI as one feature
among several engineering goals, or (c) keep it minimal and emphasize
sustainability/infrastructure. This shapes which program leads.

AI should be central

**Q3 (Claude → user):** Who is the proposing team / lead institution (PI status
matters — e.g. PESOSE requires tenured/tenure-track or equivalent PIs)? Is UCSC
the lead, and are there co-PI institutions to include?

We will decide later.

**Q4 (Claude → user):** Should I treat the FY26 deadlines above as the targets
(CSSI ~Nov 16 2026; PESOSE Sep 1 2026 / Mar 2 2027; ATI ~Nov 15 2026), or are we
aiming at a later cycle?

Stick with the FY26 deadlines.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-06-20 (Searched NSF for PypeIt-relevant funding opportunities)

Executed Task 1: searched the NSF website (nsf.gov / new.nsf.gov) for funding
opportunities suited to PypeIt, emphasizing AI and open-source angles. Ran four
targeted web searches (CSSI, POSE/PESOSE, AAG, NSF AI programs) and fetched the
live solicitation/program pages for CSSI (NSF 22-632), PESOSE (NSF 26-506), and
ATI (NSF 22-627) to confirm award classes, budget caps, and deadlines.

Findings written to the Report section, ranked by fit:
- **Tier 1:** CSSI (best fit — sustained open research software; AI is in scope;
  Nov 2026 window) and PESOSE (open-source ecosystem sustainability; successor to
  POSE; Sep 2026 / Mar 2027).
- **Tier 2 (astronomy-specific):** AAG (tools that enable astronomical research)
  and ATI (software to enable new observations; ~Nov 2026).
- **Tier 3 (AI-emphasis):** Trailblazer (AI priority, up to $3M); NAIRR (a compute
  *resource*, not a funding target); CAIG (geosciences — explicitly NOT a fit).

Key learning: NSF's open-source program rebranded from POSE → **PESOSE** (NSF
26-506) and is now a 3-track structure (not Phase I/II); CSSI remains on the
NSF 22-632 solicitation with Elements / Framework / Transition-to-Sustainability
classes. Also drafted a list of candidate AI/ML "hooks" inside PypeIt (ML
cosmic-ray rejection, learned wavelength ID, sky subtraction, frame typing,
anomaly detection, agentic reduction-diagnosis) to support an AI-emphasis
proposal. Posted four questions in the new Q&A section to direct the next step
(target program, AI centrality, proposing team/PI eligibility, target deadline).

### 2026-06-20 (Generated standalone FY2026 opportunities report)

Executed Task 2: read the user's Q&A answers and generated a separate,
self-contained report at `reports/nsf_opportunities26.md` consolidating all Task-1
findings. Incorporated the user's decisions: **AI/ML is central** (so I elevated
CSSI and the Trailblazer AI-priority award, and expanded the AI/ML thrusts
section into the proposal's research core), and **target the FY2026 deadline
cycle**. Lead program and proposing team remain TBD per the user, so the report
keeps the top candidates open and flags PI-eligibility constraints (PESOSE
tenured/tenure-track requirement; Trailblazer career-stage limits to verify).

Report structure: guiding decisions from Q&A; PypeIt framing; ranked
opportunities (Tier 1: CSSI, Trailblazer, PESOSE; Tier 2: AAG, ATI; Tier 3:
NAIRR as a compute resource, CAIG as out-of-scope); AI/ML thrusts anchored to
real reduction stages (learned wavelength ID using the arc-template/`arclines`
corpus, ML cosmic-ray rejection, sky subtraction, frame typing, anomaly
detection, agentic diagnosis); a recommendation (lead with CSSI Nov 16 2026,
Trailblazer as the single-PI AI alternative, PESOSE Track 1 as companion); open
items; and a full sources list. Created the `reports/` directory under
`claude_prompts/` (did not previously exist).