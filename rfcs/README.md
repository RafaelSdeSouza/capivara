# CAPIVARA platform design — review index

**Status:** draft architecture documents for scientific review (20 November 2026). These documents define the proposed platform; they do **not** certify current implementation behavior or authorize package changes.

| Document | Decision it specifies |
| --- | --- |
| [RFC 0001 — Scientific platform architecture](0001-platform-architecture.md) | One scientific specification, native R and Python, module boundaries, representation semantics and research isolation |
| [RFC 0002 — Scientific API and mathematical contracts](0002-scientific-api-specification.md) | Complete 41-function R export inventory, intuition and assumptions, data/result schemas, three distinct segmentation objectives, Python API proposal, research operator gates |
| [RFC 0003 — Validation, documentation and release](0003-validation-documentation-release.md) | Scientific testing and parity, website/gallery information architecture, visualization standards, citation and staged release gates |
| [RFC 0004 — Scientific editorial and geometric statistics](0004-scientific-editorial-and-geometric-statistics.md) | Research-paper quality prose, figure-first visual directions, group-action and measurement-operator compatibility, mathematical proof obligations and experimental separation |

## What is decided

- Keep the published CAPIVARA spectral segmentation science and R public API reproducible.
- Treat full-flux/shape representations, missingness/support, objective and coordinate conventions as first-class scientific contracts.
- Develop a **native** Python package against shared reference fixtures; do not claim numerical parity until verified.
- Build a documentation site that teaches the underlying physical/statistical question, mathematics, example, assumptions and failure modes for every public function.
- Maintain optional experimental geometric/spectro-spatial operators **outside** the stable segmentation workflows until they earn separate mathematical and observational validation.
- No scientific inference is justified solely by a polished figure or a passing software test.

## Open review questions (not blocking this architectural baseline)

1. Should the *new, validated observed-entry* representation workflow be recommended to new users without changing existing R function defaults? The RFC recommends yes, with explicit validation-domain qualifiers.
2. Which minimum external cube adapters are essential for the first Python preview (MaNGA first versus MaNGA + MUSE)? Proposed MaNGA first, then MUSE; no universal FITS reader claim.
3. Should the conceptual documentation have one common gateway with distinct R pkgdown and Python Sphinx reference routes? Proposed yes, with no breaking redirects.
4. Which single physically meaningful estimand will justify a future experimental spectro-spatial operator? Not yet chosen; choose after the theoretical novelty and observational benchmark tests.

## Source material (read before implementing)

- [Original MNRAS CAPIVARA method](https://doi.org/10.1093/mnras/staf688)
- [Current public API audit](../docs/public_api_audit.md)
- [Earlier broad scientific coherence audit](../docs/CAPIVARA_2_COHERENCE_AUDIT.md) — a 2026-09 proposal, not an assertion that every defect remains.
- [Capability and implementation map](../docs/CAPIVARA_2_CAPABILITY_MAP.md) — originally described an older R revision; revalidate against current code.
- [Validated wavelength/LSF baseline](../docs/SCIENCE_BASELINE_FREEZE_V2.md)
- [Current semantic representation tutorial](../vignettes/semantic-representations.Rmd)
- [The separate experimental Python branch](https://github.com/RafaelSdeSouza/capivara/tree/feat/platform-foundations-r-python-math) — feasibility spike, not approved release.

**Scope fence:** This design branch contains only `rfcs/` Markdown documents. Nothing in the R namespace, Python prototype, website, tests, release configuration or package runtime is changed by accepting or reviewing these documents.