import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, sys, json
from kernel_norm import L1_norm
from numpy.polynomial.legendre import leggauss
n = int(sys.argv[1])
xg, wg = leggauss(24)
def gl(f, a, b, log=False):
    if log:
        la, lb = np.log(a), np.log(b); t = 0.5*(lb-la)*xg + 0.5*(lb+la); s = np.exp(t)
        return 0.5*(lb-la)*np.sum(wg*s*np.array([f(x) for x in s]))
    s = 0.5*(b-a)*xg + 0.5*(b+a); return 0.5*(b-a)*np.sum(wg*np.array([f(x) for x in s]))
def Phi(z):
    f = lambda s: L1_norm(z, s, n)
    # pieces: (0,z/100) log, (z/100, z) log, (z, min(2z,1)) , (2z, 1) log
    tot = gl(f, z*1e-4, z*0.5, log=True) + gl(f, 0.5*z, z) + gl(f, z, min(1.5*z, 1.0))
    if 1.5*z < 1: tot += gl(f, 1.5*z, 1.0, log=True)
    # tiny s<z*1e-4 piece: ||K|| ~ c s for s<<z, negligible
    return tot
zs = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
coef = {1: 3*np.sqrt(3)/(4*np.pi), 2: (3/(2))*(8/3)*(0.6)**2.5/2}[n]
res = {}
for z in zs:
    P = Phi(z); res[z] = P
    print(n, z, '%.5e'%P, ' Phi/z^2=%.4f'%(P/z**2), ' (Phi/z^2)/ln(1/z)=%.4f'%(P/z**2/np.log(1/z)), ' asym coef %.4f'%coef, flush=True)
json.dump(res, open(f'{ROOT}/code/phi_profile_n{n}.json','w'))
