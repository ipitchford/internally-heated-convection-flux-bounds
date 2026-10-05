"""Stress-free slab: (D^2-k^2)^2 g = delta(z-s), g = g'' = 0 at z=0,1.
g_sf = g_hs_sf + c, g_hs_sf = G0(z-s) - G0(z+s) (odd image), c = alpha sinh(kz) + beta kz cosh(kz)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import sympy as sp, pickle
k, z, s = sp.symbols('kappa z s', positive=True)
G0 = lambda x: (1+k*sp.Abs(x))*sp.exp(-k*sp.Abs(x))/(4*k**3)
# z>s branch at z=1
ghs_top = (1+k*(z-s))*sp.exp(-k*(z-s))/(4*k**3) - (1+k*(z+s))*sp.exp(-k*(z+s))/(4*k**3)
al, be = sp.symbols('alpha beta')
c = al*sp.sinh(k*z) + be*k*z*sp.cosh(k*z)
eqs = [(c+ghs_top).subs(z,1), sp.diff(c+ghs_top, z, 2).subs(z,1)]
sol = sp.solve(eqs, [al, be], dict=True)[0]
a_sol, b_sol = sp.simplify(sol[al]), sp.simplify(sol[be])
print('alpha =', a_sol); print('beta =', b_sol)
# verify full slab function satisfies BCs and matches the sine series numerically
gl = (1+k*(s-z))*sp.exp(-k*(s-z))/(4*k**3) - (1+k*(z+s))*sp.exp(-k*(z+s))/(4*k**3)  # z<s
full_low = gl + c.subs({al:a_sol, be:b_sol})
f = sp.lambdify((k,z,s), full_low, 'mpmath')
import mpmath as mp
mp.mp.dps = 30
for (kk, zz, ss) in [(0.7, 0.2, 0.6), (3.0, 0.1, 0.35), (12.0, 0.05, 0.4)]:
    series = sum(2*mp.sin(j*mp.pi*zz)*mp.sin(j*mp.pi*ss)/((j*mp.pi)**2+kk**2)**2 for j in range(1, 40000))
    print(kk, zz, ss, mp.nstr(f(kk,zz,ss), 12), mp.nstr(series, 12))
pickle.dump((a_sol, b_sol), open((ROOT + '/code/slab_corr_sf.pkl'),'wb'))
