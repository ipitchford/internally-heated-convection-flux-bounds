"""Horizontal velocity kernels for infinite-Prandtl Stokes flow in the slab 0<z<1.
w(x,z) = R int_0^1 (K_{z,s} * T(.,s))(x) ds,  K^_{z,s}(kappa) = kappa^2 g_kappa(z,s).
bc='ns' (no-slip, clamped: g=g'=0) or 'sf' (stress-free: g=g''=0) at both walls.
Decomposition K = K_hs (half-space part, closed form) + K_c (slab correction, smooth symbol)."""
import numpy as np, mpmath as mp, pickle, sympy as sp
from scipy.special import j0
from scipy.interpolate import CubicSpline
mp.mp.dps = 50
import os
D = os.path.dirname(os.path.abspath(__file__)) + os.sep
_k, _z, _s = sp.symbols('kappa z s', positive=True)
_a, _b = pickle.load(open(D+'slab_corr.pkl','rb'))
_ns_c = _a*(_k*_z*sp.cosh(_k*_z) - sp.sinh(_k*_z)) + _b*_k*_z*sp.sinh(_k*_z)
_a2, _b2 = pickle.load(open(D+'slab_corr_sf.pkl','rb'))
_sf_c = _a2*sp.sinh(_k*_z) + _b2*_k*_z*sp.cosh(_k*_z)
_mpf = {}
for bc, cexpr in [('ns', _ns_c), ('sf', _sf_c)]:
    for d in range(3):
        _mpf[(bc,d)] = sp.lambdify((_k,_z,_s), _k**2*sp.diff(cexpr, _z, d), 'mpmath')

def _G0(x, k):  # free-space Green's function and derivatives, vectorised in k
    ax = abs(x); e = np.exp(-k*ax)
    return (1+k*ax)*e/(4*k**3), -x*e/(4*k), (k*ax-1)*e/(4*k)

def _basis(k, z, d):
    """decaying basis phi1=e^{-kz}, phi2=kz e^{-kz}, phi3=e^{-k(1-z)}, phi4=k(1-z)e^{-k(1-z)} and derivs"""
    e1, e3 = np.exp(-k*z), np.exp(-k*(1-z)); y = 1-z
    if d == 0: return np.stack([e1, k*z*e1, e3, k*y*e3], -1)
    if d == 1: return np.stack([-k*e1, k*(1-k*z)*e1, k*e3, k*(-1+k*y)*e3], -1)
    if d == 2: return np.stack([k*k*e1, k*k*(k*z-2)*e1, k*k*e3, k*k*(k*y-2)*e3], -1)

def _corr_double(k, z, s, bc, d):
    """kappa^2 d^d/dz^d (g_slab - g_hs)(z,s), k>=1, via 4x4 solve in the decaying basis"""
    rows = (0, 1) if bc == 'ns' else (0, 2)
    M = np.zeros((len(k),4,4)); rhs = np.zeros((len(k),4))
    M[:,0,:] = _basis(k, 0.0, rows[0]); M[:,1,:] = _basis(k, 0.0, rows[1])
    M[:,2,:] = _basis(k, 1.0, rows[0]); M[:,3,:] = _basis(k, 1.0, rows[1])
    if bc == 'ns':   # g_hs satisfies the bottom BCs; mismatch at z=1 (z>s)
        a, b = 1-s, 1+s
        g1 = ((1+k*a)*np.exp(-k*a) - (1+k*b+2*k*k*s)*np.exp(-k*b))/(4*k**3)
        dg1 = (-k*k*a*np.exp(-k*a) + k*k*(1-s+2*k*s)*np.exp(-k*b))/(4*k**3)
        rhs[:,2], rhs[:,3] = -g1, -dg1
    else:
        Gm, Gp = _G0(1-s, k), _G0(1+s, k)
        rhs[:,2], rhs[:,3] = -(Gm[0]-Gp[0]), -(Gm[2]-Gp[2])
    c = np.linalg.solve(M, rhs[...,None])[...,0]
    return k*k*np.sum(c*_basis(k, z, d), -1)

def corr_symbol(k, z, s, bc='ns', d=0):
    k = np.asarray(k, float); out = np.empty(len(k)); lo = k < 1.0
    f = _mpf[(bc,d)]
    for i in np.where(lo)[0]:
        out[i] = 0.0 if k[i] < 1e-9 else float(f(mp.mpf(k[i]), mp.mpf(z), mp.mpf(s)))
    if (~lo).any(): out[~lo] = _corr_double(k[~lo], z, s, bc, d)
    return out

def K_hs(y, z, s, n, bc):
    a = abs(z-s); b = z+s
    with np.errstate(divide='ignore', invalid='ignore'):
        if n == 1:
            K = np.log((b*b+y*y)/(a*a+y*y))/(8*np.pi) + (a*a/(a*a+y*y) - b*b/(b*b+y*y))/(4*np.pi)
            if bc == 'ns': K = K - z*s*(b*b-y*y)/(b*b+y*y)**2/(2*np.pi)
        else:
            K = 0.25*((a*a+y*y)**-0.5 - (b*b+y*y)**-0.5) + 0.25*(a*a*(a*a+y*y)**-1.5 - b*b*(b*b+y*y)**-1.5)
            if bc == 'ns': K = K - 0.5*z*s*(2*b*b-y*y)*(b*b+y*y)**-2.5
            K = K/(2*np.pi)
    return K

def kernel_from_symbol(symfun, n, decay, cusp=None, rmax=40.0):
    """inverse Fourier (n=1) / Hankel (n=2) transform of a smooth even symbol decaying like e^{-decay*k}.
    A |k| cusp c1*|k| at k=0 (slow 1/y^2 tail) is removed as c1*|k|e^{-|k|} and added back in closed form.
    The returned kernel is valid for |y| < 80 (n = 1, FFT period 160) or |y| <= rmax (n = 2, radial grid)."""
    kmax = 80.0/decay
    if cusp is None:
        k0 = 1e-5; cusp = float(symfun(np.array([k0]))[0])/k0
    sym = lambda k: symfun(k) - cusp*k*np.exp(-k)
    if n == 1:
        P = 160.0; N = 2**int(np.ceil(np.log2(2*kmax*P/(2*np.pi))))
        kf = 2*np.pi*np.fft.rfftfreq(N, d=P/N)
        cs = sym(np.maximum(kf, 1e-12)); cs[0] = 0.0
        Ky = np.fft.irfft(cs, n=N)*N/P
        y = np.arange(N)*P/N; m = y < P/2
        sp_ = CubicSpline(y[m], Ky[m])
        return lambda yy: sp_(yy) + cusp*(1-yy**2)/(np.pi*(1+yy**2)**2)
    dk = 0.01; k = np.arange(dk/2, kmax, dk); cs = sym(k)
    r = np.linspace(0, rmax, int(round(50*rmax))+1)
    sp_ = CubicSpline(r, (j0(np.outer(r, k)) @ (cs*k))*dk/(2*np.pi))
    return lambda rr: sp_(rr) + cusp*(2-rr**2)*(1+rr**2)**-2.5/(2*np.pi)

def L1(Kfun, n, scales):
    sc = max(min(scales), 1e-9)
    yg = np.unique(np.concatenate([[0.0], np.geomspace(sc*1e-4, 40.0, 6000), np.linspace(0, 20*max(scales), 3000)]))
    f = np.abs(Kfun(yg)) * (yg if n == 2 else 1.0); f[~np.isfinite(f)] = 0
    return (2*np.pi if n == 2 else 2.0)*np.trapezoid(f, yg)

def L1_norm(z, s, n=1, bc='ns'):
    Kc = kernel_from_symbol(lambda k: corr_symbol(k, z, s, bc, 0), n, 2-z-s)
    sc = [z, s] + ([abs(z-s)] if abs(z-s) > 0 else [])
    return L1(lambda y: K_hs(y, z, s, n, bc) + Kc(y), n, sc)

def cusp_tail(n, Y):
    """exact L1 norm of F^{-1}[kappa e^{-kappa}] outside |y| <= Y (the kernel has one sign there)"""
    return 2*Y/(np.pi*(1+Y*Y)) if n == 1 else Y*Y/(1+Y*Y)**1.5

def L1_corr_deriv(z, s, n=1, bc='ns', d=2, Y=40.0, return_parts=False):
    """||F^{-1}[kappa^2 d^d/dz^d c(kappa; z, s)]||_{L^1(R^n)} for the slab correction c.
    The |kappa| cusp c1*kappa (present for bc='sf', where c1 = -s/2) is split off as c1*kappa*e^{-kappa}; the
    kernel is integrated over |y| <= Y and the exact tail of the cusp kernel outside Y is added, which is
    correct when |c1| dominates the remainder kernel there (checked by tail_check below). Version 1 of this
    routine truncated at |y| = 40 without the tail."""
    cusp = -s/2 if (bc == 'sf' and d == 1) else 0.0
    Kc = kernel_from_symbol(lambda k: corr_symbol(k, z, s, bc, d), n, 2-z-s, cusp=cusp, rmax=max(40.0, Y))
    sc = [1.0]
    yg = np.unique(np.concatenate([[0.0], np.geomspace(1e-5, Y, 8000), np.linspace(0, Y, 8001)]))
    f = np.abs(Kc(yg))*(yg if n == 2 else 1.0)
    core = (2*np.pi if n == 2 else 2.0)*np.trapezoid(f, yg)
    # |kappa|^3 term of the remainder symbol: kernel ~ c3 * (6/(pi y^4)) (n=1) or c3 * 1.4324/r^5 (n=2) at large y
    c3 = s*s if (bc == 'ns' and d == 2) else (s*(1 - s*s - 3*z*z)/4 if (bc == 'sf' and d == 1) else 0.0)
    tail3 = abs(c3)*(4/(np.pi*Y**3) if n == 1 else 2*np.pi*(8/np.pi)*0.75*np.sqrt(np.pi)/(4*np.sqrt(np.pi)/3)/(3*Y**3))
    tail = abs(cusp)*cusp_tail(n, Y) + tail3
    if return_parts:
        return core + tail, core, tail, float(Kc(np.array([Y]))[0])
    return core + tail
