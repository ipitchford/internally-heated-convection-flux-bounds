"""Continuation of steady 2D infinite-Pr rolls in Rayleigh number with aspect-ratio optimisation.
mode 'rb': maximise Nu over Gamma (Rayleigh-Benard, no-slip by default)
mode 'ih': minimise F_B over Gamma (internally heated, isothermal walls)
Each solve: Newton with LU reuse (chord iterations), refactorising when convergence stalls."""
import os, numpy as np, scipy.linalg as sla, json, time, sys, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady2d import Steady

def res_norm(st, X): return np.linalg.norm(st.residual(X))/np.sqrt(st.n)

def solve(st, X, tol=1e-9, maxit=40, verbose=False):
    lu = None; r = res_norm(st, X); it = 0; nfac = 0
    while it < maxit and r > tol:
        if lu is None:
            lu = sla.lu_factor(st.jacobian(X), check_finite=False); nfac += 1
        G = st.residual(X)
        dX = -sla.lu_solve(lu, G.ravel(), check_finite=False).reshape(st.shape)
        lam = 1.0; ok = False
        while lam > 1e-3:
            Xn = X + lam*dX; rn = res_norm(st, Xn)
            if np.isfinite(rn) and rn < r*(1-0.25*lam): ok = True; break
            lam *= 0.5
        if not ok:
            if nfac and lu is not None and it > 0:
                lu = None; it += 1; continue        # refactorise at current X
            return X, False, r
        if rn > 0.3*r: lu = None                     # slow contraction: refactorise next time
        X, r = Xn, rn; it += 1
        if verbose: print(f'    it {it} |G|={r:.2e} lam={lam:.2f}', flush=True)
    return X, r <= tol, r

def resolution(Ra, mode):
    lr = np.log10(Ra)
    if lr <= 4.5: return 32, 32
    if lr <= 5.5: return 48, 40
    if lr <= 6.5: return 64, 56
    if lr <= 7.5: return 80, 72
    if lr <= 8.5: return 96, 96
    return 128, 128

class Scanner:
    def __init__(self, mode, bc, out):
        self.mode, self.bc, self.out = mode, bc, out
        self.cache = {}; self.records = []

    def evaluate(self, Ra, G, X0, Nz, Nx):
        st = Steady(R=Ra, Gamma=G, Nx=Nx, Nz=Nz, bc=self.bc, rb=(self.mode == 'rb'))
        if X0 is not None and X0.shape != st.shape:
            X0 = Steady.interp(X0, X0.shape[0]-1, Nz, st.K)
        X, ok, r = solve(st, X0)
        d = st.diagnostics(X) if ok else None
        obj = None
        if ok: obj = d['Nu'] if self.mode == 'rb' else -d['FB']     # maximise obj
        return X, ok, obj, d

    def optimise(self, Ra, Gc, X0, Nz, Nx, width=0.25, tolG=0.01):
        """golden-section maximisation of obj over log Gamma in [Gc e^-width, Gc e^width]"""
        evals = {}
        def f(lg, Xs):
            G = float(np.exp(lg))
            X, ok, obj, d = self.evaluate(Ra, G, Xs, Nz, Nx)
            if not ok: return -np.inf, Xs, None
            evals[G] = (obj, d); return obj, X, d
        a, b = np.log(Gc)-width, np.log(Gc)+width; gr = (np.sqrt(5)-1)/2
        c = b - gr*(b-a); dpt = a + gr*(b-a)
        fc, Xc, _ = f(c, X0); fd, Xd, _ = f(dpt, X0 if not np.isfinite(fc) else Xc)
        # expand bracket if optimum at an edge
        for _ in range(6):
            if fc >= fd and np.isfinite(fc):
                fa, Xa, _ = f(a, Xc)
                if fa > fc: b, dpt, fd, Xd = dpt, c, fc, Xc; c = a; fc, Xc = fa, Xa; a = a - width; continue
            break
        while b - a > tolG:
            if fc >= fd:
                b, dpt, fd, Xd = dpt, c, fc, Xc; c = b - gr*(b-a); fc, Xc, _ = f(c, Xd)
            else:
                a, c, fc, Xc = c, dpt, fd, Xd; dpt = a + gr*(b-a); fd, Xd, _ = f(dpt, Xc)
        if fc >= fd: return float(np.exp(c)), Xc, evals
        return float(np.exp(dpt)), Xd, evals

    def run(self, Ra_list, G0, X0):
        Gc, X = G0, X0
        for Ra in Ra_list:
            Nz, Nx = resolution(Ra, self.mode); t0 = time.time()
            Gs, Xs, evals = self.optimise(Ra, Gc, X, Nz, Nx)
            obj, d = evals[Gs]
            rec = dict(Ra=Ra, Gamma=Gs, Nz=Nz, Nx=Nx, wall=time.time()-t0, **{k: float(v) for k, v in d.items()},
                       scan={f'{g:.5f}': float(v[0]) for g, v in sorted(evals.items())})
            self.records.append(rec)
            print(json.dumps({k: rec[k] for k in ['Ra', 'Gamma', 'Nz', 'Nx', 'Nu', 'FB', 'FT', 'Tmax', 'wall']}), flush=True)
            json.dump(self.records, open(self.out, 'w'), indent=1)
            np.save(self.out.replace('.json', f'_X_{Ra:.3e}.npy'), Xs)
            Gc, X = Gs, Xs

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--mode'); p.add_argument('--bc', default='ns')
    p.add_argument('--lo', type=float); p.add_argument('--hi', type=float); p.add_argument('--step', type=float, default=0.25)
    p.add_argument('--G0', type=float); p.add_argument('--out')
    a = p.parse_args()
    Ras = 10**np.arange(np.log10(a.lo), np.log10(a.hi)+1e-9, a.step)
    sc = Scanner(a.mode, a.bc, a.out)
    # initial state near onset
    Nz, Nx = resolution(Ras[0], a.mode)
    st = Steady(R=Ras[0], Gamma=a.G0, Nx=Nx, Nz=Nz, bc=a.bc, rb=(a.mode == 'rb'))
    z = st.s.z; X0 = np.zeros(st.shape)
    X0[:, 0] = (1-z) if a.mode == 'rb' else z*(1-z)/2
    X0[:, 1] = (0.1 if a.mode == 'rb' else 0.02)*np.sin(np.pi*z)
    sc.run(Ras, a.G0, X0)
