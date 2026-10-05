"""Phi_bc,n(z) = int_0^1 ||K_{z,s}||_1 ds on a grid of z (uses kernels.py)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, sys, json
from kernels import L1_norm
from numpy.polynomial.legendre import leggauss
n, bc = int(sys.argv[1]), sys.argv[2]
xg, wg = leggauss(24)
def gl(f, a, b, log=False):
    if log:
        la, lb = np.log(a), np.log(b); s = np.exp(0.5*(lb-la)*xg + 0.5*(lb+la))
        return 0.5*(lb-la)*np.sum(wg*s*np.array([f(x) for x in s]))
    s = 0.5*(b-a)*xg + 0.5*(b+a); return 0.5*(b-a)*np.sum(wg*np.array([f(x) for x in s]))
def Phi(z):
    f = lambda s: L1_norm(z, s, n, bc)
    tot = gl(f, z*1e-4, z*0.5, log=True) + gl(f, 0.5*z, z) + gl(f, z, min(1.5*z, 1.0))
    if 1.5*z < 1: tot += gl(f, 1.5*z, 1.0, log=True)
    return tot
zs = [float(x) for x in sys.argv[3].split(',')]
res = {}
for z in zs:
    P = Phi(z); res[z] = P
    print(n, bc, z, '%.6e'%P, 'Phi/z=%.5f'%(P/z), 'Phi/z^2=%.5f'%(P/z**2), flush=True)
json.dump(res, open(f'{ROOT}/code/phi_{bc}_n{n}_{len(zs)}.json','w'))
