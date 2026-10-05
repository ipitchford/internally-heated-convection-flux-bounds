"""Convergence and independent checks of the slab-correction norms at the computed maximisers.
(1) the production routine with the truncation radius Y doubled;
(2) n = 1: an independent evaluation, K(y) = (1/pi) int_0^inf mu(kappa) cos(kappa y) dkappa by mpmath
    quadrature (cusp part in closed form), and ||K||_1 by quadrature between the sign changes of K."""
import numpy as np, mpmath as mp, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import L1_corr_deriv, corr_symbol
mp.mp.dps = 20
JS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'corr_convergence.json')
cases = [('ns', 1, 0.0, 0.9125), ('sf', 1, 0.0, 1.0), ('sf', 1, 0.0, 0.95), ('ns', 2, 0.0, 0.9), ('sf', 2, 0.0, 1.0), ('ns', 1, 0.2, 0.6)]
# optional: "python corr_convergence.py n2" recomputes only the n = 2 cases and merges them into the stored file
only = sys.argv[1] if len(sys.argv) > 1 else None
out = json.load(open(JS)) if (only and os.path.exists(JS)) else {}
for bc, n, z, s in cases:
    if only == 'n2' and n != 2: continue
    d = 2 if bc == 'ns' else 1
    v40 = L1_corr_deriv(z, s, n, bc, d, Y=40.0); v80 = L1_corr_deriv(z, s, n, bc, d, Y=80.0)
    rec = dict(bc=bc, n=n, zeta=z, s=s, Y40=v40, Y80=v80)
    if n == 1:
        cusp = -s/2 if bc == 'sf' else 0.0
        rem = lambda k: float(corr_symbol(np.array([float(k)]), z, s, bc, d)[0]) - cusp*float(k)*np.exp(-float(k))
        Kc = lambda y: mp.quad(lambda k: rem(k)*mp.cos(k*y), mp.linspace(0, 80, 41))/mp.pi + cusp*(1-y*y)/(mp.pi*(1+y*y)**2)
        ys = np.concatenate([np.linspace(0, 6, 241), np.linspace(6.05, 40, 120)])
        Kv = [float(Kc(y)) for y in ys]
        # locate sign changes, integrate |K| piecewise with Gauss-Legendre
        roots = [mp.findroot(Kc, (ys[i], ys[i+1]), solver='anderson') for i in range(len(ys)-1) if Kv[i]*Kv[i+1] < 0]
        pts = [0] + [float(r) for r in roots] + [40.0]
        L1 = 2*sum(abs(mp.quad(Kc, [pts[i], pts[i+1]])) for i in range(len(pts)-1))
        # analytic tails beyond y = 40: the cusp part exactly, and the |kappa|^3 term of the low-frequency
        # expansion, whose kernel decays like 6 c3/(pi y^4) (c3 from lowfreq_check.py)
        tail = abs(cusp)*2*40/(np.pi*(1+40**2))
        c3 = s*s if bc == 'ns' else s*(1-s*s-3*z*z)/4
        tail3 = abs(c3)*4/(np.pi*40.0**3)
        rec.update(independent=float(L1)+tail+tail3, independent_tail_cusp=tail, independent_tail_k3=tail3,
                   sign_changes=len(roots))
    out[f'{bc}_n{n}_{z}_{s}'] = rec
    print(json.dumps(rec), flush=True)
json.dump(out, open(JS, 'w'), indent=1)
