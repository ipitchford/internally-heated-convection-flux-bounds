# Superseded and exploratory files

Kept for the record; nothing in the paper depends on them.

| File | What it was | Superseded by |
|---|---|---|
| `corr_constants.py`, `corr_constants_n{1,2}.json`, `corr_constants_n2.log` | Version-1 grid suprema of the slab-correction norms. The norms were truncated at radius 40 without analytic tails, which put C_sf about 0.008 (d = 2) and 0.012 (d = 3) below the exact lower bound ν₁/2 found by the referee. | `../corr_constants2.py` (analytic cusp and \|κ\|³ tails; endpoint identity as a test) |
| `phi_profile.py`, `phi_profile_n{1,2}.json`, `phi_profile_n2.log`, `phi_sf_n1_10.json` | Early, coarser wall-function profiles. | `../phi_profile2.py` and `../phi_*_26.json` |
| `carlson_probe.py` | Carlson majorants at six parameter points. | `../carlson_sup.py` (suprema over the parameter rectangle) |
| `expo_form.py` | Probe showing that a termwise exponential bound on the symbols fails because of cancellations. | the proof of Lemma 3.6(ii) (Neumann series in a decaying basis) |
