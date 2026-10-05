"""Benchmark against Blankenbach et al. (1989), cases 1a-1c: isoviscous infinite-Prandtl convection in a unit
square with free-slip walls, Ra = 1e4, 1e5, 1e6. The square cell is half of our mirror-symmetric roll pair of
period 2 between stress-free plates. Reference Nusselt numbers (cases 1a-1c): 4.884409, 10.534095, 21.972465."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json
from steady2d import Steady
from steady_scan import solve
ref = {1e4: 4.884409, 1e5: 10.534095, 1e6: 21.972465}
X = None; out = []
for Ra, (Nz, Nx) in [(1e4, (48, 64)), (1e5, (64, 96)), (1e6, (96, 128)), (1e6, (128, 160))]:
    st = Steady(R=Ra, Gamma=2.0, Nx=Nx, Nz=Nz, bc='sf', rb=True)
    if X is None:
        z = st.s.z; X = np.zeros(st.shape); X[:, 0] = 1-z; X[:, 1] = 0.1*np.sin(np.pi*z)
    elif X.shape != st.shape:
        X = Steady.interp(X, X.shape[0]-1, Nz, st.K)
    X, ok, r = solve(st, X); d = st.diagnostics(X)
    rel = d['Nu']/ref[Ra]-1
    out.append(dict(Ra=Ra, Nz=Nz, Nx=Nx, ok=bool(ok), Nu=d['Nu'], ref=ref[Ra], rel_err=rel))
    print(f"Ra={Ra:.0e} (Nz,Nx)=({Nz},{Nx}) converged={ok} Nu={d['Nu']:.8f} ref={ref[Ra]} rel.err={rel:+.2e}", flush=True)
json.dump(out, open((ROOT + '/data/steady/blankenbach.json'), 'w'), indent=1)
