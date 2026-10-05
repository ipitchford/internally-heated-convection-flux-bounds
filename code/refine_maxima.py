"""Refine both local maxima of Nu(Gamma) at Ra = 1e5 (Rayleigh-Benard, no-slip, Pr = infinity) at two resolutions,
by continuation in Gamma from the stored narrow-roll solution and golden-section search on each maximum."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady2d import Steady
from steady_scan import solve
Ra = 1e5; X0 = np.load((ROOT + '/data/steady/rb_ns_X_1.000e+05.npy')); G0 = 0.731508
out = []
for Nz, mpu in [(48, 100), (64, 140)]:
    cache = {}
    def nu_at(G, Xs):
        Nx = int(max(32, 8*np.ceil(mpu*G/8))); st = Steady(R=Ra, Gamma=G, Nx=Nx, Nz=Nz, bc='ns', rb=True)
        X = Steady.interp(Xs, Xs.shape[0]-1, Nz, st.K) if Xs.shape != st.shape else Xs
        X, ok, r = solve(st, X); assert ok, (G, r); return st.diagnostics(X)['Nu'], X
    # continuation along the family, storing solutions
    Gs = np.exp(np.linspace(np.log(0.70), np.log(1.10), 15)); X = X0; sols = {}
    for G in Gs: nu, X = nu_at(G, X); sols[G] = (nu, X)
    for lo, hi in [(0.68, 0.82), (0.95, 1.10)]:
        a, b = np.log(lo), np.log(hi); gr = (np.sqrt(5)-1)/2
        def f(lg):
            G = float(np.exp(lg)); Gn = min(sols, key=lambda g: abs(np.log(g/G)))
            nu, X = nu_at(G, sols[Gn][1]); sols[G] = (nu, X); return nu
        c, d = b-gr*(b-a), a+gr*(b-a); fc, fd = f(c), f(d)
        while b-a > 2e-3:
            if fc >= fd: b, d, fd = d, c, fc; c = b-gr*(b-a); fc = f(c)
            else: a, c, fc = c, d, fd; d = a+gr*(b-a); fd = f(d)
        Gm = float(np.exp(c if fc >= fd else d)); o = dict(Ra=Ra, Nz=Nz, modes_per_unit=mpu, branch='narrow' if hi < 0.9 else 'wide', Gamma=Gm, Nu=max(fc, fd))
        out.append(o); print(json.dumps(o), flush=True)
json.dump(out, open((ROOT + '/data/steady/maxima_1e5.json'), 'w'), indent=1)
