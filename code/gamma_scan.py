"""Continuation in aspect ratio Gamma at fixed Rayleigh number, to map Nu(Gamma) (RB) or F_B(Gamma) (IH)
and locate all local optima. Starts from a stored steady state of the Ra-continuation scan."""
import os, numpy as np, json, sys, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady2d import Steady
from steady_scan import solve
p = argparse.ArgumentParser()
p.add_argument('--mode'); p.add_argument('--bc', default='ns'); p.add_argument('--Ra', type=float)
p.add_argument('--start'); p.add_argument('--G0', type=float); p.add_argument('--Gmin', type=float); p.add_argument('--Gmax', type=float)
p.add_argument('--Nz', type=int); p.add_argument('--modes_per_unit', type=float); p.add_argument('--fac', type=float, default=1.04)
p.add_argument('--out')
a = p.parse_args()
X0 = np.load(a.start)
def nx_for(G): return int(max(32, 8*np.ceil(a.modes_per_unit*G/8)))
recs = []
for direction in (+1, -1):
    G, X = a.G0, X0.copy()
    while a.Gmin <= G <= a.Gmax:
        Nx = nx_for(G)
        st = Steady(R=a.Ra, Gamma=G, Nx=Nx, Nz=a.Nz, bc=a.bc, rb=(a.mode == 'rb'))
        if X.shape != st.shape: X = Steady.interp(X, X.shape[0]-1, a.Nz, st.K)
        t0 = time.time(); Xn, ok, r = solve(st, X)
        if not ok:
            print(f'G={G:.4f} failed (|G|={r:.1e}); stopping this direction', flush=True); break
        d = st.diagnostics(Xn); X = Xn
        rec = dict(Ra=a.Ra, Gamma=G, Nx=Nx, Nz=a.Nz, wall=time.time()-t0, **{k: float(v) for k, v in d.items()})
        recs.append(rec)
        print(json.dumps({k: rec[k] for k in ['Gamma', 'Nx', 'Nu', 'FB', 'Tmax', 'wall']}), flush=True)
        json.dump(sorted(recs, key=lambda q: q['Gamma']), open(a.out, 'w'), indent=1)
        G = G*a.fac if direction > 0 else G/a.fac
