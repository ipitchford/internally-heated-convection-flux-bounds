"""Evaluate the lower bound  F_B >= B(Lambda) = int_0^1 z e^{-Lambda P(z)} dz / int_0^1 e^{-Lambda P(z)} dz,
P(z) = int_0^z Phi(zeta) dzeta,  Lambda = R*Tsup/2,  for no-slip (ns) and stress-free (sf) walls."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys
from scipy.interpolate import CubicSpline
from scipy.special import gamma as G
nu2 = {1: 3*np.sqrt(3)/(2*np.pi), 2: 1.115419}

def load_phi(bc, n):
    d = json.load(open(f'{ROOT}/code/phi_{bc}_n{n}_26.json'))
    z = np.array(sorted(float(k) for k in d)); P = np.array([d[k] for k in sorted(d, key=float)])
    return z, P

def make_Phi(bc, n):
    z, P = load_phi(bc, n)
    if bc == 'ns':
        a = nu2[n]/2; g = P/z**2 - a*np.log(1/z)          # smooth remainder
        b0 = g[0]
        spl = CubicSpline(np.log(z), g)
        def Phi(zz):
            zz = np.maximum(np.minimum(zz, 1-zz), 1e-300); out = np.empty_like(zz)
            lo = zz < z[0]; hi = ~lo
            out[lo] = zz[lo]**2*(a*np.log(1/zz[lo]) + b0)
            out[hi] = zz[hi]**2*(a*np.log(1/zz[hi]) + spl(np.log(zz[hi])))
            return out
        return Phi, a, b0
    else:
        c0 = P[0]/z[0]; spl = CubicSpline(np.log(z), P/z)
        def Phi(zz):
            zz = np.minimum(zz, 1-zz); out = np.empty_like(zz)
            lo = zz < z[0]; out[lo] = c0*zz[lo]; out[~lo] = zz[~lo]*spl(np.log(zz[~lo]))
            return out
        return Phi, c0, None

def B_of_Lambda(Phi, Lams):
    zg = np.unique(np.concatenate([np.geomspace(1e-9, 0.5, 20000), np.linspace(0.5, 1, 4001)]))
    ph = Phi(zg)
    P = np.concatenate([[0], np.cumsum(0.5*(ph[1:]+ph[:-1])*np.diff(zg))])
    P += ph[0]*zg[0]/3     # tiny initial piece (Phi ~ z^2 or z near 0; negligible)
    out = []
    for L in Lams:
        E = np.exp(-L*P)
        out.append(np.trapezoid(zg*E, zg)/np.trapezoid(E, zg))
    return np.array(out)

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    Lams = 10.0**np.arange(1, 12.01, 0.5)
    res = {}
    for bc in ['ns', 'sf']:
        Phi, A, b = make_Phi(bc, n)
        B = B_of_Lambda(Phi, Lams)
        if bc == 'ns':
            asym = (G(2/3)/G(1/3))*(9/(A*Lams*np.log(Lams)))**(1/3)
        else:
            asym = 1/np.sqrt(np.pi*Lams*A/2)
        res[bc] = dict(Lambda=Lams.tolist(), B=B.tolist(), asym=asym.tolist(), A=A, b=b)
        print(f'--- n={n} bc={bc}  coefficient A={A:.5f}' + (f'  b={b:.4f}' if b is not None else ''))
        for L, x, y in zip(Lams[::2], B[::2], asym[::2]):
            extra = (x*(L*np.log(L))**(1/3)) if bc=='ns' else x*np.sqrt(L)
            print(f'Lambda=1e{np.log10(L):4.1f}  B={x:.4e}  asym={y:.4e}  ratio={x/y:.3f}  '
                  + ('B*(L lnL)^(1/3)' if bc=='ns' else 'B*sqrt(L)') + f'={extra:.3f}')
    json.dump(res, open(f'{ROOT}/code/bound_curve_n{n}.json','w'), indent=1)
