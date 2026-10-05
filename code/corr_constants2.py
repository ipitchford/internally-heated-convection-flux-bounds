"""Slab-correction constants C_ns = sup ||F^{-1}[kappa^2 d^2_zeta c_ns]||_1 and C_sf = sup ||F^{-1}[kappa^2 d_zeta c_sf]||_1
over 0 <= zeta <= 1/2, 0 <= s <= 1, for n = d-1 horizontal dimensions (estimates, not interval-certified).
Coarse grid, then two levels of local refinement around the grid maximum. Also records the norm of the
cusp-free remainder (for sf) and the endpoint value (zeta, s) = (0, 1), which is exactly nu_1/2 for sf."""
import numpy as np, sys, json, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import L1_corr_deriv
H = os.path.dirname(os.path.abspath(__file__))
n = int(sys.argv[1]); out = {}
for bc, d in [('ns', 2), ('sf', 1)]:
    t0 = time.time(); vals = {}
    f = lambda zeta, s: vals.setdefault((round(zeta, 6), round(s, 6)), L1_corr_deriv(zeta, s, n, bc, d))
    Z, S = np.linspace(0, 0.5, 11), np.linspace(0, 1, 21)
    for zeta in Z:
        for s in S: f(zeta, s)
    for h in [0.025, 0.00625]:                    # local refinement around the current maximum
        (zb, sb), vb = max(vals.items(), key=lambda kv: kv[1])
        for zeta in np.clip(zb + h*np.arange(-2, 3), 0, 0.5):
            for s in np.clip(sb + h*np.arange(-2, 3), 0, 1): f(float(zeta), float(s))
    (zb, sb), vb = max(vals.items(), key=lambda kv: kv[1])
    out[bc] = dict(sup=vb, argmax=[zb, sb], n_evals=len(vals), seconds=time.time()-t0,
                   endpoint_0_1=f(0.0, 1.0), grid={f'{k[0]:.5f},{k[1]:.5f}': v for k, v in vals.items()})
    print(n, bc, 'sup = %.7f at (zeta, s) = (%.5f, %.5f); endpoint (0,1): %.7f; %d evaluations, %.0f s'
          % (vb, zb, sb, out[bc]['endpoint_0_1'], len(vals), time.time()-t0), flush=True)
json.dump(out, open(os.path.join(H, f'corr_constants2_n{n}.json'), 'w'), indent=1)
