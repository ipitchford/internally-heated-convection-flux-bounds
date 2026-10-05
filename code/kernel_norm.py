"""Numerical evaluation of Phi_n(z) = int_0^1 ||K_{z,s}||_{L^1(R^n)} ds for the no-slip slab.
w(x,z) = R int_0^1 (K_{z,s} * T(.,s))(x) ds, symbol of K_{z,s}: m(kappa)=kappa^2 g_kappa(z,s).
K = K_hs (closed form, half-space) + K_c (slab correction; smooth symbol, FFT/Hankel quadrature)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, mpmath as mp, pickle, sympy as sp, sys
from scipy.special import j0
from scipy.interpolate import CubicSpline
mp.mp.dps = 50
kk, zz_, ss_ = sp.symbols('kappa z s', positive=True)
a_sol, b_sol = pickle.load(open((ROOT + '/code/slab_corr.pkl'),'rb'))
Aexp = kk*zz_*sp.cosh(kk*zz_) - sp.sinh(kk*zz_); Bexp = kk*zz_*sp.sinh(kk*zz_)
corr_mp = sp.lambdify((kk, zz_, ss_), kk**2*(a_sol*Aexp + b_sol*Bexp), 'mpmath')

def corr_symbol_double(k, z, s):
    """kappa^2 * (g_slab - g_hs)(z,s) via 4x4 solve in decaying exponential basis; k array (k>=1)."""
    e0 = np.exp(-k); 
    # basis values/derivs at z=0 and z=1
    M = np.zeros((len(k),4,4))
    # phi1=e^{-kz}, phi2=kz e^{-kz}, phi3=e^{-k(1-z)}, phi4=k(1-z)e^{-k(1-z)}
    M[:,0,:] = np.stack([np.ones_like(k), 0*k, e0, k*e0],1)
    M[:,1,:] = np.stack([-k, k, k*e0, k*(-1+k)*e0],1)
    M[:,2,:] = np.stack([e0, k*e0, np.ones_like(k), 0*k],1)
    M[:,3,:] = np.stack([-k*e0, k*(1-k)*e0, k, -k],1)
    a, b = 1-s, 1+s
    g1 = ((1+k*a)*np.exp(-k*a) - (1+k*b+2*k*k*s)*np.exp(-k*b))/(4*k**3)
    dg1 = (-k*k*a*np.exp(-k*a) + k*k*(1-s+2*k*s)*np.exp(-k*b))/(4*k**3)
    rhs = np.stack([0*k, 0*k, -g1, -dg1],1)
    c = np.linalg.solve(M, rhs[...,None])[...,0]
    val = c[:,0]*np.exp(-k*z) + c[:,1]*k*z*np.exp(-k*z) + c[:,2]*np.exp(-k*(1-z)) + c[:,3]*k*(1-z)*np.exp(-k*(1-z))
    return k*k*val

def corr_symbol(k, z, s):
    out = np.empty(len(k))
    lo = k < 1.0
    for i in np.where(lo)[0]:
        out[i] = 0.0 if k[i] < 1e-8 else float(corr_mp(mp.mpf(k[i]), mp.mpf(z), mp.mpf(s)))
    if (~lo).any(): out[~lo] = corr_symbol_double(k[~lo], z, s)
    return out

def K_hs_1d(y, z, s):
    a = abs(z-s); b = z+s
    with np.errstate(divide='ignore'):
        L = np.log((b*b+y*y)/(a*a+y*y))
    return (L/(8*np.pi) + (a*a/(a*a+y*y) - b*b/(b*b+y*y))/(4*np.pi)
            - z*s*(b*b-y*y)/(b*b+y*y)**2/(2*np.pi))

def K_hs_2d(r, z, s):
    a = abs(z-s); b = z+s
    with np.errstate(divide='ignore'):
        t1 = 0.25*((a*a+r*r)**-0.5 - (b*b+r*r)**-0.5)
        t2 = 0.25*(a*a*(a*a+r*r)**-1.5 - b*b*(b*b+r*r)**-1.5)
    t3 = -0.5*z*s*(2*b*b-r*r)*(b*b+r*r)**-2.5
    return (t1+t2+t3)/(2*np.pi)

def Kc_on_grid(z, s, n):
    decay = 2-z-s; kmax = 80.0/decay
    if n == 1:
        P = 160.0; N = 2**int(np.ceil(np.log2(2*kmax*P/(2*np.pi))))
        kf = 2*np.pi*np.fft.rfftfreq(N, d=P/N)
        cs = corr_symbol(np.maximum(kf,1e-12), z, s); cs[0] = 0.0
        Ky = np.fft.irfft(cs, n=N) * N / P   # (1/2pi) sum m e^{iky} dk, dk=2pi/P
        y = np.arange(N)*P/N; m = y < P/2
        return CubicSpline(y[m], Ky[m])
    else:
        dk = 0.01; k = np.arange(dk/2, kmax, dk)
        cs = corr_symbol(k, z, s)
        r = np.linspace(0, 40, 2001)
        Kr = (j0(np.outer(r, k)) @ (cs*k))*dk/(2*np.pi)
        return CubicSpline(r, Kr)

def L1_norm(z, s, n=1):
    Kc = Kc_on_grid(z, s, n)
    sc = min(z, s, abs(z-s)) if abs(z-s) > 0 else min(z,s)
    yg = np.unique(np.concatenate([[0.0], np.geomspace(max(sc,1e-9)*1e-4, 40.0, 6000),
                                   np.linspace(0, 20*max(z,s), 3000)]))
    if n == 1:
        K = K_hs_1d(yg, z, s) + Kc(yg)
        f = np.abs(K); f[~np.isfinite(f)] = 0
        return 2*np.trapezoid(f, yg)
    else:
        K = K_hs_2d(yg, z, s) + Kc(yg)
        f = np.abs(K)*yg; f[~np.isfinite(f)] = 0
        return 2*np.pi*np.trapezoid(f, yg)

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv)>1 else 1
    import time; t0=time.time()
    for z in [0.2, 0.1, 0.05]:
        out = [L1_norm(z, s, n) for s in [z/4, z/2, 2*z, 4*z, 0.5, 0.9]]
        print(n, z, ['%.4e'%v for v in out], ' 4s/z^2 at s=4z: %.3f'%(out[3]*4*z/z**2), '%.1fs'%(time.time()-t0), flush=True)
