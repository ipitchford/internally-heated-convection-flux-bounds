import numpy as np, mpmath as mp, scipy.linalg as sl, pickle, sympy as sp
from ipc2d import IPC2D
# 1) velocity map vs explicit Green's function: T_hat(z) = sin(pi z) z^2, mode k
k_, z_, s_ = sp.symbols('kappa z s', positive=True)
a_sol, b_sol = pickle.load(open('slab_corr.pkl','rb'))
def ghs(zz, ss):
    a = abs(zz-ss); b = zz+ss
    return ((1+k_*a)*sp.exp(-k_*a) - (1+k_*b+2*k_**2*zz*ss)*sp.exp(-k_*b))/(4*k_**3)
corr = a_sol*(k_*z_*sp.cosh(k_*z_) - sp.sinh(k_*z_)) + b_sol*k_*z_*sp.sinh(k_*z_)
cf = sp.lambdify((k_, z_, s_), corr, 'mpmath')
mp.mp.dps = 30
def g_slab(kk, zz, ss):
    a = abs(zz-ss); b = zz+ss; kk = mp.mpf(kk)
    hs = ((1+kk*a)*mp.exp(-kk*a) - (1+kk*b+2*kk**2*zz*ss)*mp.exp(-kk*b))/(4*kk**3)
    lo, hi = min(zz,ss), max(zz,ss)
    return hs + cf(kk, mp.mpf(lo), mp.mpf(hi))   # correction formula derived for z<s; symmetric in (z,s)
sim = IPC2D(R=1.0, Gamma=2*np.pi/3.0, Nx=8, Nz=48, bc='ns')
j = 1; kk = sim.k[j]; f = np.sin(np.pi*sim.z)*sim.z**2
w_num = sim.Aw[j] @ f
for zi in [3, 10, 24, 40]:
    zz = float(sim.z[zi])
    w_ex = kk**2*mp.quad(lambda ss: g_slab(kk, zz, ss)*mp.sin(mp.pi*ss)*ss**2, [0, zz, 1])
    print('k=%.3f z=%.4f  w_num=%.10e  w_exact=%.10e' % (kk, zz, w_num[zi], float(w_ex)))
# 2) linear stability: RB infinite Pr no-slip (Ra_c=1707.76 at k=3.117), IH no-slip and stress-free
def growth(R, k, bc, rb, Nz=48):
    s = IPC2D(R=R, Gamma=2*np.pi/k, Nx=8, Nz=Nz, bc=bc, rb=rb)
    N = s.Nz; I = np.eye(N+1)
    dTc = -np.ones(N+1) if rb else (0.5 - s.z)
    Lop = s.D2 - k*k*I - np.diag(dTc) @ s.Aw[1]
    # Dirichlet: restrict to interior
    Li = Lop[1:N, 1:N]
    ev = np.linalg.eigvals(Li)
    return ev.real.max()
from scipy.optimize import brentq
for bc, rb, k0, Rg in [('ns', True, 3.117, (1000, 3000)), ('ns', False, 4.0, (20000, 60000)), ('sf', False, 3.0, (3000, 40000))]:
    def crit(k):
        return brentq(lambda R: growth(R, k, bc, rb), *Rg)
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(crit, bracket=(k0*0.8, k0, k0*1.2), tol=1e-6)
    print(bc, 'RB' if rb else 'IH', 'R_c = %.2f at k = %.4f' % (res.fun, res.x))
