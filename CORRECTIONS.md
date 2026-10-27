# Corrections

## 2026-10-27 — errata to internally-heated-convection-flux-bounds 2.1.0-candidate (DOI 10.5281/zenodo.23172948)

The archived PDF, the manuscript source, the Zenodo deposit and the tag `v2.1.0-candidate` are unchanged. These errata were confirmed by the Evidence Press archive audit and are also shown at <https://evidencepress.org/releases/internally-heated-convection-flux-bounds/>.

1. Corrected the bound quoted in the abstract, the introduction and the discussion. The archived text states T ≤ C(1 + ln |Ω|). Theorem 2.3 and its proof use the enlarged torus Ω̃, with each period enlarged to ⌈1/Γ_i⌉Γ_i ≥ 1 (so Ω̃ = Ω when every Γ_i ≥ 1). These three statements should read C(1 + ln |Ω̃|); as printed, the bound would be negative for |Ω| < 1/e.

   Theorem 2.3 itself is stated and proved correctly; only the summaries are affected. No theorem, certificate, computation or stated result is affected.

2. Corrected a resolution threshold in the numerical section. The archived text says the resolution was raised to (96, 96) “above Ra = 3 × 10⁷”. Table 4 uses (80, 72) at Ra = 3.16 × 10⁷ and first uses (96, 96) at Ra = 5.62 × 10⁷, so the text should say “and to (96, 96) from Ra = 5.62 × 10⁷”.

   The table and the fitted prefactor are correct. No theorem, certificate, computation or stated result is affected.
