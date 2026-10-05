# Where each result comes from

Each row names a result in the paper, the script that produces it, the stored data it reads or writes, and the log
that records the run. `code/verify_numbers.py` checks every number in the right-hand columns against the text and
exits with a nonzero status if any differs (log: `logs/verify_numbers.log`).

## Theory and constants

| Result in the paper | Script (in `code/`) | Data written or read | Log |
|---|---|---|---|
| Green's functions, Appendix A | `slab_green.py`, `slab_green_sf.py`, `check_symbols.py` | `slab_corr.pkl`, `slab_corr_sf.pkl` | `logs/check_symbols.log` |
| Constants ν₁, ν₂, ν₃, Appendix B | `verify_numbers.py` (recomputed in 30-digit arithmetic) | — | `logs/verify_numbers.log` |
| Lemma 3.6(i): low-frequency expansions, cusp −(s/2)κ, endpoint identity μ_sf(κ;0,1) = −½κe^{−κ} | `lowfreq_check.py`; rechecked in `verify_numbers.py` | — | `logs/lowfreq_check.log` |
| Slab-correction suprema C_ns, C_sf (Section 5.1) | `corr_constants2.py 1`, `corr_constants2.py 2` (uses `kernels.py`) | `corr_constants2_n1.json`, `corr_constants2_n2.json` | `logs/corr_constants2_n{1,2}.log` |
| Truncation radius doubling and independent mpmath evaluation | `corr_convergence.py` (then `corr_convergence.py n2`, see note) | `corr_convergence.json` | `logs/corr_convergence.log`, `logs/corr_convergence_n2.log` |
| Carlson majorants of Lemma 3.6(iii) | `carlson_sup.py` | `carlson_sup.json` | `logs/carlson_sup.log` |
| Constants A, B, κ₀, c_ns, c_sf of Theorem 2.1 (both routes); 𝒯_max of Theorem 2.3 | `constants.py` | `paper/facts_constants.json` | `logs/constants.log` |
| Wall functions Φ, Lemma 3.4, Figure 1 | `phi_profile2.py` (via `run_profiles.sh`), `kernel_norm.py` | `phi_ns_n{1,2}_26.json`, `phi_sf_n{1,2}_26.json` | `code/prof_*.log` |
| Wall constant c_w, Proposition 2.2 | `cw_compute.py` | `cw_values.json` | `logs/cw_compute.log` |
| Bound functional 𝔅, Table 2, Figure 2 | `bound_curve.py`, `make_tables.py`, `fig_theory.py` | `paper/facts_bounds.json`, `paper/tab_bound.tex` | `logs/make_tables.log`, `logs/fig_theory.log` |
| Interval-arithmetic cross-check of the symbols (not used for any quoted number) | `cert/symbols_arb.py`, `cert/test_symbols_arb.py` | — | `logs/test_symbols_arb.log` |

## Simulations

| Result in the paper | Script | Data | Log |
|---|---|---|---|
| Time-stepper and its validation (linear thresholds, influence matrix) | `ipc2d.py`, `test_ipc.py` | — | `logs/test_ipc.log` |
| Simulations at R = 10⁵, 10⁶ (both boundary conditions) | `run_dns.py` (queued by `queue_dns.sh`) | `data/dns/ih_*_cont.json`, `_hist.npz`, `_prof.npz`, `_state.npz`, `_snap.npz`; the earlier windows without `_cont` | `data/dns/*.log` |
| Table 3: physical-time averages, block means, heat balance, ratios, finite-window bound; flow character (plateaus, bursts) | `check_dns.py` | `data/dns/dns_checks.json` | `logs/check_dns.log` |
| Spectral resolution of the simulations | `dns_resolution.py` (reads the snapshots in the companion archive) | `data/dns/resolution.json` | `logs/dns_resolution.log` |
| Profile hypothesis of Chen et al. (full windows, sub-blocks, steady states) | `check_chen.py`; sub-block statistics recomputed in `verify_numbers.py` | `data/dns/chen_hypothesis_checks.json` | `logs/check_chen.log` |
| Figure 3 | `fig_dns.py` | `figs/fig3_ih_ns_R1e+06_G4_256x64_cont.*` | `logs/fig_dns.log` |

## Steady rolls

| Result in the paper | Script | Data | Log |
|---|---|---|---|
| Newton solver and its validation (Jacobian, agreement with time-stepping) | `steady2d.py`, `test_steady.py` | — | `logs/test_steady.log` |
| Blankenbach benchmark | `test_blankenbach.py` | `data/steady/blankenbach.json` | `logs/test_blankenbach.log` |
| Continuation in Ra with optimisation over the period (Table 4) | `steady_scan.py` | `data/steady/ih_ns.json`, `rb_ns.json`, `*_X_*.npy` | `data/steady/*.log` |
| Largest maxima over the scanned periods; refined maxima at Ra = 10⁵ | `rb_global.py`, `refine_maxima.py` | `data/steady/rb_ns_best_found.json`, `maxima_1e5.json` | `logs/rb_global.log`, `data/steady/maxima_1e5.log` |
| Continuation in the period at fixed Ra; shoulder (Figure 4c,d) | `gamma_scan.py`, `steady_shoulder.py` | `data/steady/gscan_*.json` | `data/steady/gscan_*.log` |
| Resolution checks and spectral tails | `res_check.py`, `res_check2.py`, `res_check3.py`, `spectral_tails.py` | `data/steady/res_check*.json`, `spectral_tails.json` | `logs/spectral_tails.log` |
| Fits and comparisons (exponents, fit sensitivity, factor above the bound) | `make_tables.py` | `paper/facts_steady.json`, `paper/tab_steady.tex` | `logs/make_tables.log` |
| Figure 4 | `fig_steady.py` | `figs/fig4_steady.*` | `logs/fig_steady.log` |

## Paper

| Product | Script | Log |
|---|---|---|
| Reference list (APA 7 via citeproc); every DOI checked against its Crossref record | `make_refs.py`, `check_dois.py` | `logs/make_refs.log`, `logs/check_dois.log` |
| PDF | `pdflatex main.tex` (twice) in `paper/` | `logs/latex_build.log` |
| Markdown | `make_md.py` | `logs/make_md.log` |

**Note on `corr_convergence`.** In the first full run (`logs/corr_convergence.log`) the radius-80 values for the
two-dimensional kernels (n = 2) were wrong: the radial grid of the Hankel transform stopped at 40, and the spline
was extrapolated beyond it. `kernels.kernel_from_symbol` now takes the radial grid out to the truncation radius.
The n = 2 cases were rerun with `corr_convergence.py n2` (`logs/corr_convergence_n2.log`), which replaces those
entries in `corr_convergence.json`. A fresh full run gives the same file. The radius-40 values, which are the ones
used in the paper, were not affected.
