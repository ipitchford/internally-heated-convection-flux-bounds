# Novelty-gate report: pointwise lower bounds on the bottom heat flux in internally heated convection at infinite Prandtl number

**Search date:** 5 October 2026
**Researcher/agent:** Claude (agent run for Evidence Press)
**Corpus boundary:** Semantic Scholar citation graph (citing papers of Arslan & Rojas 2024/2025; Whitehead & Doering 2011; Goluskin & Spiegel 2012; Arslan 2024 rotating; Arslan et al. 2024 heating profiles); arXiv API full-listing search `abs:"internally heated" AND abs:bounds` (all dates, sorted by date, 40 most recent); general web search (US index); full texts inspected for arXiv 2403.14407, 2406.10975, 2604.01380, 2602.21860, 2402.19240, 2102.06458, 2208.04421, 0911.2294, 2008.08752, 2212.05407, 2504.11051, nlin/0004020. No access to MathSciNet or zbMATH full search; Springer page for Chen et al. (2024) rate-limited.
**Intended contribution units:** theorem + complete proof (bounds); verification artefact (numerics); computation (steady rolls).
**Revision:** rerun on 5 October 2026 for the revised manuscript (section 11), after the referee audit. The claims were narrowed (constants computed but not interval-certified; steady-roll statements restricted to the computed roll family); no new collision was found.

## 1. Frozen claim

### Exact statement
For classical solutions of uniform internally heated Boussinesq convection at infinite Prandtl number in a horizontally periodic slab with isothermal plates (T = 0) and nonnegative initial temperature, with long-time maximum temperature 𝒯 and flux Rayleigh number R:
(a) no-slip plates: liminf time-averaged F_B ≥ 𝔅[(R𝒯/2)Φ_ns] ≥ c (Λ ln Λ)^{-1/3}, Λ = R𝒯/2;
(b) stress-free plates: ≥ c Λ^{-1/2};
(c) for any divergence-free tangential drift, 𝒯 ≤ 𝒯_max = 1 + (e/2π²) max{2, ln(2^{d+1}|Ω̃|)}, where Ω̃ is the periodic cell enlarged to have every period at least 1;
hence F_B ≳ (R ln R)^{-1/3} and ≳ R^{-1/2} at fixed aspect ratio.

### Exclusions
Not claimed: aspect-ratio-uniform bounds; finite Prandtl number; attainment of the bound by solutions; rigorous construction of Nu ≥ c Ra^{1/3} solutions in Rayleigh–Bénard convection.

### Permitted claim language before the gate
"improves the exponents of the best bounds known to us".

## 2. Normalisation and translation table

| Original notation | Canonical notation | Transformation checked |
|---|---|---|
| F_B = ⟨∂_zT⟩_h at z=0 (Arslan & Rojas 2025) | same | identical nondimensionalisation (length d, time d²/κ, temperature d²H/κρc_p) |
| ⟨wT⟩ ≤ 1/2 − c f(R) | F_B ≥ c f(R) | F_B = 1/2 − ⟨wT⟩ (time-averaged) |
| "convective heat flux asymmetry" (Wang et al., 2604.08164) | ⟨wT⟩ | same quantity |
| expected exit time τ_u (Iyer et al. 2010) | steady T for heat equation with unit source | T = τ_u when u is steady |

## 3. Fingerprints

Non-computational target (inequalities). Distinctive formula fingerprints searched: exponents 1/3 and 1/2 for F_B; "(R ln R)^{-1/3}"; wall coefficient ν₂/2 = 3√3/(4π) ≈ 0.4135 (d = 2) and 0.5577 (d = 3); half-space clamped Green's function in the form [(1+κa)e^{-κa} − (1+κb+2κ²zs)e^{-κb}]/(4κ³); the identity ⟨wT⟩_h ≤ max|w| ⟨T⟩_h.

## 4. Alias map

| Analysis | Probability | Mathematical physics | Geophysics | Other |
|---|---|---|---|---|
| lower bound on boundary flux; maximum/minimum principle; log-Lipschitz Stokes estimate | expected exit time with incompressible drift; hot spots | auxiliary-functional / background method; upper bound theory | radiogenically/volumetrically heated mantle; downward heat flux; stable bottom layer | optimal cooling; heat-flux asymmetry |

## 5. Search log

| Date | Corpus/database | Query | Results inspected | Consequence |
|---|---|---|---|---|
| 2026-10-05 | Semantic Scholar | citations of arXiv:2403.14407 | 3 (2604.01380, 2602.21860, 2503.20515) | none prove F_B bounds by a new method |
| 2026-10-05 | Semantic Scholar | citations of Whitehead & Doering 2011 (JMP 52, 093101) | 31 | IHC bound papers all use background/auxiliary functionals |
| 2026-10-05 | Semantic Scholar | citations of Goluskin & Spiegel 2012 | 51 | as above; numerical/phenomenological papers |
| 2026-10-05 | arXiv API | abs:"internally heated" AND abs:bounds | 30 entries | two 2026 entries checked (2604.08164 phenomenology; 2608.14334 optimal cooling) — no collision |
| 2026-10-05 | full text | 2604.01380 (rotating QG IHC, 2026) | sections 3–4 | quadratic auxiliary functional; authors report failing to get ⟨wθ⟩ bound below 1/2 — no collision |
| 2026-10-05 | full text | 2406.10975 (Arslan 2024) | §4 | uses Yan (2004) per-Fourier-mode pointwise estimate within background method; best non-rotating F_B ≳ R^{-2/3} — no collision |
| 2026-10-05 | full text | 2208.04421 (Song, Fantuzzi, Tobasco) | abstract, refs | dissipation bounds for balanced sources/sinks; different object |
| 2026-10-05 | full text | 0911.2294 (Iyer et al.) | Thms 1.1–1.2 | steady L^∞ bound by ball of equal volume on bounded domains: related to Theorem 2.3 (time-dependent slab variant); cite |
| 2026-10-05 | web | "R^{-1/3}" internally heated bottom flux minimum principle pointwise | ~20 | no collision |
| 2026-10-05 | web | infinite Prandtl steady rolls Nusselt-maximising aspect ratio | ~10 | Wen, Goluskin & Doering 2022 (Pr = 1); Wen et al. 2022/2026 (stress-free, Pr ≤ 10^{3/2}); no infinite-Pr no-slip study found |
| 2026-10-05 (rerun) | Crossref abstracts; academia.edu full text | Sondak, Smith & Waleffe (2015), JFM 784, 565–595 | abstract; passages on Pr range, the two maxima, exponents | **near-collision** for the steady-roll unit (see §6b) |
| 2026-10-05 (rerun) | full text (arXiv 2212.05407v2) | Wen, Takla, Goluskin & Chini (2026) | §§1, 3, 4, 6 | their Ra^{3/10} asymptotics and the exponent 0.298 are for **stress-free** walls (Pr = 0.1 data); they speculate the same scaling holds for no-slip walls |
| 2026-10-05 (rerun) | full text (open access) | Motoki, Kawahara & Shimizu (2021), JFM 914, A14 | abstract, §§1, 3, 5 | 3D multi-scale steady solution at Pr = 1, Nu ∼ Ra^{0.31} up to Ra ∼ 10^7; relevant to the construction question; not infinite Pr |
| 2026-10-05 (rerun) | Crossref abstract | Goluskin & van der Poel (2016) | abstract | finite-Pr IH simulations; F_B decrease reverses in 2D turbulent regime; cite |
| 2026-10-05 (rerun) | arXiv 2606.14138v4 | Ding, Wen & Li (2026), marginal energy stability | abstract, introduction | finite-Pr RB, different object; no collision |
| 2026-10-05 (rerun) | Crossref | Roberts (1979); Jiménez & Zufiria (1987) | metadata, abstract | classical boundary-layer theory of steady infinite-Pr cells at fixed aspect ratio; no Nu-maximising-period study |

## 6. Verified collisions or bridges

None for the unconditional bounds. One **conditional near-collision** for the exponent:
- Chen, Gao, Liao, Wan & Liu (2024), Acta Mech. Sin. 40, 323630, §3.2–3.3, Eqs. (14)–(22). For uniform heating at *finite* Pr, they combine the rigorous lower bound ⟨T⟩ ≥ 4.421R^{-1/3} (Lu, Doering & Busse 2004) with the profile inequality (14), T̄(z) ≤ min(q_b z, (1−z)(1−q_b)) (their q_b = our F_B), to get ⟨wT⟩ ≤ (1/4 − 8.842R^{-1/3})^{1/2}, i.e. F_B ≥ 1/2 − (1/4 − 8.842R^{-1/3})^{1/2} ≥ 8.842R^{-1/3}. Inequality (14) is introduced by "assuming that the horizontally averaged temperature profile satisfies the following inequality", justified only "in an immediately neighbouring domain of the wall" and by numerical results; their conclusions list it as a limitation. The exponent 1/3 therefore appears in the literature **conditionally on an unproven hypothesis**. Translation check: (14) ⟺ ∫_0^z J ≤ z²/2 and ∫_z^1 J ≥ −(1−z)²/2 for all z, with J = ⟨wT⟩_h. (Their text says the bound approaches 1/2 "at a rate of R^{-1/6}"; since 1/2 − √(1/4 − x) ≥ x, the implied rate is R^{-1/3}.) Consequence: our contribution is the first *unconditional* bound with exponent 1/3 found in this search (with a logarithm and at infinite Pr); the paper must cite Chen et al. and state the comparison, and can test hypothesis (14) numerically.

Related (to be cited, not colliding):
- Constantin & Doering (1999): maximum principle plus logarithmic estimate of the Stokes kernel for Rayleigh–Bénard convection at infinite Prandtl number — same family of ideas, different problem; the IHC flux bound and the integrating-factor argument are not there.
- Iyer, Novikov, Ryzhik & Zlatoš (2010); Berestycki et al. (2009): L^∞ bounds for drift-diffusion with incompressible drift (steady, bounded domains).
- Wen, Goluskin & Doering (2022): numerically Nu ∼ Ra^{1/3} for Nu-maximising steady rolls at Pr = 1.

## 6b. Rerun for the steady-roll unit (triggered by the computed Nu(Γ) shoulder)

- **Sondak, Smith & Waleffe (2015)** computed 2D no-slip steady rolls that maximise Nu over the horizontal wavenumber for Pr = 1, 4, 7, 10, 100 and Ra ≲ 10^9. Load-bearing passages (full text via academia.edu): "For Pr = 1, Pr = 4 and Pr = 7 the second local maximum does not overtake the first maximum, whereas for Pr = 10 and Pr = 100 the second maximum eventually becomes the global maximum"; "The optimal solutions exhibit a scaling of the form Nu ∼ Ra^{0.31}"; the second (higher-wavenumber) maximum has wavenumber ≈ 0.257 Ra^{0.256}. Infinite Pr was not computed. Translation: their "second local maximum" = our narrow-roll (Γ ∝ Ra^{-1/4}) branch; their "first" = the wide-roll branch (Γ ∝ Ra^{-1/5}).
- Consequence: the statement "at large Pr the narrow-roll maximum becomes the global optimum" is **known** (COLLISION, cite Sondak et al.). What remains ours, within the bounded search: (i) the computation at Pr = ∞ itself; (ii) the finding that at Pr = ∞ the wide-roll maximum, which coexists with the narrow one at Ra = 10^5 (and is still global there, by 0.013%), has disappeared by Ra = 1.8 × 10^5, leaving a shoulder at ΓRa^{1/5} ≈ 10.0–10.4 whose slope grows in magnitude from 0.012 to 0.115 between Ra = 1.8 × 10^5 and 10^7; (iii) the exponent 0.298 over 10^6 ≤ Ra ≤ 10^8 for no-slip walls, matching the stress-free Ra^{3/10} asymptotics of Wen et al. (2026) and supporting their conjecture that the scaling carries over to no-slip walls. Outcome for (i)–(iii): BOUNDED UNCERTAINTY. Sondak et al.'s Pr = 100 curves may already show the disappearance of the wide maximum; the summary of their figure on ResearchGate says so, but we did not find that statement in the passages of their text that we inspected.
- Wording permitted: "consistent with and extending to infinite Prandtl number the finding of Sondak et al. (2015)"; "not found in the bounded search". Forbidden: "first to show that narrow rolls are optimal at large Pr".

## 7. Negative evidence and missingness audit

- Chen, Gao, Liao, Wan & Liu (2024), Acta Mech. Sin. 40, 323630: full text supplied by the user and inspected (all 11 pages); see section 6.
- No MathSciNet/zbMATH search. Russian/Chinese-language literature not searched.
- Geodynamics literature on steady cells at infinite Prandtl number (rigid boundaries) is large; a Nu-maximising-aspect-ratio study may exist under other names.
- No author enquiry made (and, by standing instruction, none suggested).

## 8. Outcome by contribution unit

| Unit | Outcome | Confidence | Allowed wording | Forbidden wording |
|---|---|---:|---|---|
| Statement (F_B exponents 1/3, 1/2) | CLEAR TO PROVE (bounded), unconditional version only | moderate–high | "improves the best unconditional bounds known to us"; "the exponent 1/3 was previously obtained only under the profile hypothesis of Chen et al. (2024)" | "first R^{-1/3} bound", "previously unknown exponent" |
| Complete proof | CLEAR TO PROVE (bounded) | moderate–high | as above | — |
| L^∞ lemma (Theorem 2.3) | BRIDGE / incremental | high | "a time-dependent variant of known L^∞ bounds for incompressible drifts" | "new" without qualification |
| Near-wall kernel lemma | incremental (log is classical) | high | "explicit form with sharp wall coefficient" | "new logarithmic estimate" |
| Numerical steady rolls at Pr = ∞ | BOUNDED UNCERTAINTY (after rerun; narrow-roll optimality at large Pr is a COLLISION with Sondak et al. 2015) | moderate | "extends to infinite Prandtl number the large-Pr trend found by Sondak et al. (2015)"; "we did not find a computation of Nu-maximising no-slip rolls at infinite Prandtl number" | "first computation"; "first to find narrow-roll optima" |

## 9. Decision

**CLEAR TO PROVE (within the recorded boundary)** for the main bounds; **BRIDGE** for the L^∞ lemma; **BOUNDED UNCERTAINTY** for the steady-roll computations.

## 10. Rerun triggers

- If the log factor is removed (e.g. via a Chanillo–Malchiodi-type bilinear estimate), rerun against the Rayleigh–Bénard log-removal literature.
- If a Γ-uniform temperature bound is proved, rerun against the exit-time/hot-spot literature (Iyer et al.; Marcotte et al. 2018).
- If steady rolls are found to attain Nu ∼ Ra^{1/3} at infinite Pr, rerun against geodynamics steady-cell literature.

## 11. Revision rerun (5 October 2026, revised manuscript)

Trigger: the referee audit of version 1 and the resulting changes (repaired stress-free lemma, corrected constants, unconditional Figure 2 curve, narrowed steady-roll wording). The frozen claim in section 1 is unchanged except for the explicit temperature ceiling 𝒯_max.

| Date | Corpus/database | Query | Results inspected | Consequence |
|---|---|---|---|---|
| 2026-10-05 | arXiv API, sorted by submission date | abs:"internally heated convection" | 15 most recent (to 2604.08164) | no new bound paper since the first search |
| 2026-10-05 | arXiv API | abs:"internally heated" AND abs:"Prandtl" | 15 most recent; new entry 2608.28756 (Nicoski et al., Aug 2026) | rotating quasi-geostrophic simulations; no bound; no collision |
| 2026-10-05 | arXiv API | abs:"heat flux" AND abs:"minimum principle" | 1 (2110.10344, Kumar et al.) | background method with minimum principle; exponentially small bound; already cited |
| 2026-10-05 | arXiv API | abs:"internal heating" AND abs:"bound" AND abs:"convection" | 15 (incl. 2211.02587, Bouillaut et al.) | upper bounds on Nu for sources and sinks; different object; no collision |
| 2026-10-05 | arXiv abstracts | 2205.03175 (Arslan et al. 2023), 2110.10344 (Kumar et al. 2022) | abstracts | background method plus Hardy–Rellich inequalities; R^{-2} and exponentially small bounds; cited |
| 2026-10-05 | web | internally heated convection infinite Prandtl lower bound bottom heat flux 2026 | ~10 | no new result |
| 2026-10-05 | web | "internally heated convection" bound "R^{-1/3}" minimum principle proof | ~10 | no collision |
| 2026-10-05 | web | follow-ups to Chen et al. (2024) proving their profile hypothesis | ~10 | none found |

Equation-level comparison with Chen et al. (2024), from the full text: their Eq. (14) is the triangle hypothesis T̄ ≤ T_st(z, z₀) with z₀ = 1 − q_b; Eq. (15), q_b − 1 ≤ dT̄/dz ≤ q_b, is the sufficient condition from which they derive it; Eq. (21) is ⟨T⟩ ≤ 1/8 − ½⟨wT⟩²; Eq. (22) combines (21) with ⟨T⟩ ≥ 4.421R^{-1/3} (their Table 1, from Lu, Doering & Busse 2004). Hypothesis (14) is supported by a Taylor expansion at the walls and by numerical results, and the use of (14) is listed among the limitations in their conclusions. Translation unchanged from section 6.

Arb (interval) arithmetic was attempted for the slab-correction constants and abandoned for this version; the constants are reported as computed values with error control (paper, Section 5.1), and the paper's Table 1 (status) and claim language say so.

Outcome by unit: unchanged from section 8. The permitted wording now also includes "with constants computed numerically, without interval certification" wherever numerical values of c_ns, c_sf or C_ns, C_sf are quoted.
