"""The stress-free wall coefficient c_w = lim Phi_sf(z)/z = int_0^1 ||F^{-1}[kappa^2 d_z g_sf(kappa; 0, s)]||_1 ds.
d = 2 (n = 1), method A: closed-form physical-space kernel from the sine series,
    k(y; s) = (1/2) sum_j sin(j pi s) (1 - j pi |y|) e^{-j pi |y|}, summed in closed form with q = e^{-pi|y|}.
d = 2 and 3, method B: half-space part in closed form plus the slab correction from kernels.py (cusp
    subtracted analytically), integrated over |y| <= 40 (the total kernel decays exponentially).
Gauss-Legendre in s on [0, 1] with a graded split near s = 0."""
import numpy as np, sys, os, json, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT
from kernels import kernel_from_symbol, corr_symbol
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(40)
def gl(f, a, b):
    s = 0.5*(b-a)*xg + 0.5*(b+a); return 0.5*(b-a)*np.sum(wg*np.array([f(x) for x in s]))
def s_integral(f):
    return gl(f, 0, 0.02) + gl(f, 0.02, 0.2) + gl(f, 0.2, 0.6) + gl(f, 0.6, 1.0)
# ---- method A (d = 2)
def kA(y, s):
    q = np.exp(-np.pi*np.abs(y)); th = np.pi*s
    D = 1 - 2*q*np.cos(th) + q*q
    S0 = q*np.sin(th)/D
    S1 = q*np.sin(th)*(1 - q*q)/D**2           # sum_j j sin(j th) q^j = q dS0/dq
    return 0.5*(S0 - np.pi*np.abs(y)*S1)
def L1A(s):
    y = np.unique(np.concatenate([np.geomspace(1e-7, 40, 40000), np.linspace(0, 40, 40001)]))
    return 2*np.trapezoid(np.abs(kA(y, s)), y)
# ---- method B (d = 2, 3)
def L1B(s, n):
    Kc = kernel_from_symbol(lambda k: corr_symbol(k, 0.0, s, 'sf', 1), n, 2-s, cusp=-s/2)
    if n == 1:
        khs = lambda y: (s/2)*(s*s - y*y)/(np.pi*(s*s + y*y)**2)
    else:
        khs = lambda r: (s/2)*(2*s*s - r*r)/(2*np.pi*(s*s + r*r)**2.5)
    y = np.unique(np.concatenate([np.geomspace(max(s, 1e-6)*1e-5, 40, 40000), np.linspace(0, 40, 40001)]))
    f = np.abs(khs(y) + Kc(y))*(y if n == 2 else 1.0)
    return (2*np.pi if n == 2 else 2.0)*np.trapezoid(f, y)
out = {}
cwA = s_integral(L1A); out['cw_d2_methodA'] = cwA; print('c_w (d=2), closed-form sine-series kernel: %.6f' % cwA, flush=True)
cwB2 = s_integral(lambda s: L1B(s, 1)); out['cw_d2_methodB'] = cwB2; print('c_w (d=2), half-space + correction:   %.6f' % cwB2, flush=True)
cwB3 = s_integral(lambda s: L1B(s, 2)); out['cw_d3_methodB'] = cwB3; print('c_w (d=3), half-space + correction:   %.6f' % cwB3, flush=True)
json.dump(out, open(os.path.join(ROOT, 'code', 'cw_values.json'), 'w'), indent=1)
