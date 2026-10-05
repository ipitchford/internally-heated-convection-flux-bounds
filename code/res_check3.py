"""Resolution convergence check for steady states: recompute selected scan points at higher resolution."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady2d import Steady
from steady_scan import solve
cases = [('rb', 1e8, '../data/steady/rb_ns.json', '../data/steady/rb_ns_X_1.000e+08.npy', [(128, 128)])]
out = []
for mode, Ra, jf, xf, resl in cases:
    rec = [r for r in json.load(open(jf)) if abs(r['Ra']/Ra-1) < 1e-6][0]
    X0 = np.load(xf); G = rec['Gamma']
    for Nz, Nx in resl:
        st = Steady(R=Ra, Gamma=G, Nx=Nx, Nz=Nz, bc='ns', rb=(mode == 'rb'))
        X = Steady.interp(X0, X0.shape[0]-1, Nz, st.K)
        X, ok, r = solve(st, X); d = st.diagnostics(X)
        # spectral tail: max |coefficient| in last 10% of Fourier modes relative to max
        tail = np.abs(X[:, int(0.9*st.K):]).max()/np.abs(X[:, 1:]).max()
        o = dict(mode=mode, Ra=Ra, Gamma=G, Nz=Nz, Nx=Nx, ok=bool(ok), Nu=d['Nu'], FB=d['FB'], FT=d['FT'], fourier_tail=float(tail))
        out.append(o); print(json.dumps(o), flush=True)
json.dump(out, open((ROOT + '/data/steady/res_check3.json'), 'w'), indent=1)
