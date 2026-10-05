# Heat escaping through the base of internally heated convection at infinite Prandtl number

**Evidence Press, 5 October 2026 (revised after referee audit).**

Version 2.1.0-candidate. DOI: https://doi.org/10.5281/zenodo.23172948.

This bundle contains the manuscript, supporting code and stored computational evidence. Written continuum proofs are not formally verified. The paper proves a minimum-principle
bound F_B ≳ (R ln R)^{-1/3} and reports finite-range numerical tests with steady rolls. It addresses item 3, heat transport in
mantle-type convection, of *Open problems in geology suited to frontier AI*.

## Contents

| Path | What it is |
|---|---|
| `paper/main.pdf` | The paper. |
| `paper/paper.md` | Markdown version of the same text, with APA 7 references. Figures load from `../figs/`. |
| `paper/*.tex`, `paper/tab_*.tex`, `paper/refs_apa.tex` | LaTeX sources. `main.tex` inputs the section files; the tables and reference list are generated. |
| `paper/references.bib`, `references.json`, `apa.csl` | Bibliography (BibTeX and CSL-JSON) and the APA 7 citation style used to format it. |
| `paper/crossref.json` | Crossref metadata used to check the DOIs. |
| `paper/facts_*.json` | Numbers quoted in the text, as written by `constants.py` and `make_tables.py`. |
| `paper/novelty_report.md` | Novelty-gate report: search log (including the rerun for this revision), near-collisions, and the claim wording each permits. |
| `figs/` | Figures 1–4 as PDF (used by the paper) and PNG (used by the Markdown). |
| `code/` | All code (see `DEPENDENCIES.md`), with small intermediate data files: wall-function profiles, slab-correction grids, Carlson majorants, convergence checks, pickled symbolic Green's-function corrections. `code/cert/` holds the interval-arithmetic evaluators, used only as a cross-check; `code/superseded/` keeps the version-1 constant computation and other superseded or exploratory scripts, with a note on each. |
| `data/dns/` | Simulations. The continuation runs (`*_cont.*`) are the primary data: summaries (`.json`), time histories at every step (`_hist.npz`), mean and sub-block profiles (`_prof.npz`), final states (all eight repaired snapshots are included here). Also the version-1 spin-up runs, and the derived checks `dns_checks.json`, `chen_hypothesis_checks.json` and `resolution.json`. |
| `data/steady/` | Steady rolls: Rayleigh-number continuations (`rb_ns.json`, `ih_ns.json`, solution arrays `*_X_*.npy`), the table of largest maxima `rb_ns_best_found.json`, continuations in the period (`gscan_*.json`), refined maxima at Ra = 10^5, resolution checks, spectral tails, the benchmark test and logs. |
| `logs/` | Output of every check, test and build, including `verify_numbers.log`. |
| `requirements.txt`, `DEPENDENCIES.md`, `reproduce.sh` | Pinned versions; map from each result to its script, data and log; one-command rerun. |
| `MANIFEST.sha256` | SHA-256 checksum of every file in the bundle. |

## Main results

1. **Theorem (unconditional).** For uniform internal heating at infinite Prandtl number between isothermal plates,
   every classical solution satisfies F̄_B ≥ c (R𝒯 ln R𝒯)^{-1/3} (no-slip) and F̄_B ≥ c (R𝒯)^{-1/2}
   (stress-free), where 𝒯 is the long-time maximum temperature. In addition, 𝒯 ≤ 𝒯_max = 1 + (e/2π²) max{2,
   ln(2^{d+1}|Ω̃|)} for every incompressible flow. At fixed aspect ratio this improves the exponents of the best
   unconditional bounds from 2/3 to 1/3 (up to a logarithm) and from 40/29 to 1/2. The exact large-R asymptotics
   of the bound functional are determined.
2. **Status of the constants.** The constants are explicit functions of two slab-correction norms, C_ns and C_sf,
   which are proved finite. Their numerical values use sensitivity checks and a second route through Carlson majorants. Neither route supplies uniform error enclosures; the decimal values are estimates. With the computed values,
   c_ns ≈ 0.205 and c_sf ≈ 0.270 (d = 2); with the majorants, 0.204 and 0.254.
3. **Tests on computed flows.** Simulations and steady rolls satisfy every inequality of the proof with wide margins,
   and stay far above the bound: over the computed ranges, F_B ≈ 1.3 R^{-0.09} for both the published
   time-dependent simulations and F_B-minimising steady rolls. One finite-window profile violates the Chen triangle off-grid; this does not refute the stationary hypothesis.
4. **Rayleigh–Bénard at infinite Prandtl number.** Within the mirror-symmetric, layer-spanning roll family and the
   periods scanned, the wide-roll maximum of Nu over the period is not detected above the sampled interval between Ra = 10^5 and 1.8 × 10^5. Above
   that, the largest maxima found are narrow rolls with Nu ≈ 0.153 Ra^{0.298} up to Ra = 10^8. A rigorous
   construction with Nu ≥ c Ra^{1/3} remains open.

## Reproducing the results

Requirements: Python 3.11 or later with the packages in `requirements.txt` (the results were produced with
Python 3.13 and the pinned versions), pdflatex (TeX Live: natbib, amsthm, mathtools, hyperref, microtype, lmodern) and
pandoc 3.x. All paths are relative to the bundle.

```
./reproduce.sh          # about 20 minutes: tests, all derived checks, constants, tables, figures, paper, assertion suite
./reproduce.sh --full   # also recomputes the slab-correction grids and their convergence checks (about 1.5 hours more)
```

The script ends by running `code/verify_numbers.py`, an assertion suite that checks every number quoted in the paper
against the stored data (recomputing it where that is cheap) and exits with a nonzero status on any failure. The
simulations (`code/queue_dns.sh`, a few minutes to a quarter of an hour each) and the steady-roll continuations
(`code/steady_scan.py`, `code/gamma_scan.py`; 10–25 minutes per Rayleigh number above Ra = 3 × 10^7) are not rerun
by `reproduce.sh`; their stored outputs are what the checks read. `DEPENDENCIES.md` lists, for every result, the
script, data and log behind it.

## Notes

- Third-party papers consulted during the work are not redistributed. They are cited in the paper, whose reference
  list gives a DOI or URL for each; `paper/novelty_report.md` records how each was found and used.
- Scholarly creator: Anonymous. Publisher: Evidence Press. The work was produced by directed AI agents.

## Revision and replay boundaries

Read [REVIEW_RESPONSE_V21.md](REVIEW_RESPONSE_V21.md), [AI_INDEX.md](AI_INDEX.md), and [CITATION_AUDIT_V21.md](CITATION_AUDIT_V21.md) first. The supplied novelty report is a historical agent report; its older numerical-certification and branch-disappearance wording is superseded by this revision. Private review prose and the original problem prompt are not redistributed. Current environment and scoped replay logs are retained separately from original logs. The eight snapshots now use the physical 384-point horizontal grid; `x_base` preserves the 256-point spectral grid.

A short current check is `cd code && python3 revision_checks.py && python3 verify_numbers.py`. Both work under `python3 -O`. This checks stored arrays and arithmetic, not a fresh full simulation campaign. The full reproduction script remains a development workflow that regenerates outputs. No decimal prefactor is certified as a uniform bound.
