"""Slab Green's function for (D^2-k^2)^2 g = delta(z-s) on (0,1), clamped: g=g'=0 at z=0,1.
Written as g = g_hs (half-space, wall at 0) + correction c = a*A(z) + b*B(z),
A = k z cosh kz - sinh kz, B = k z sinh kz (both vanish to 2nd order at z=0)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import sympy as sp
k, z, s = sp.symbols('kappa z s', positive=True)
# half-space Green's function, z>s branch and z<s branch
def ghs(zz, ss, branch):
    a = (zz-ss) if branch=='z>s' else (ss-zz)
    b = zz+ss
    return ((1+k*a)*sp.exp(-k*a) - (1+k*b+2*k**2*zz*ss)*sp.exp(-k*b))/(4*k**3)
# check ODE and BCs
for br in ['z>s','z<s']:
    G = ghs(z, s, br)
    L = sp.simplify(sp.diff(G,z,4) - 2*k**2*sp.diff(G,z,2) + k**4*G)
    print(br, 'L g =', L)
G0 = ghs(z,s,'z<s')
print('g(0)=', sp.simplify(G0.subs(z,0)), " g'(0)=", sp.simplify(sp.diff(G0,z).subs(z,0)))
# jump conditions at z=s: g'' continuous, g''' jumps by 1
Gp, Gm = ghs(z,s,'z>s'), ghs(z,s,'z<s')
for j in range(4):
    print('jump d^%d:'%j, sp.simplify((sp.diff(Gp,z,j)-sp.diff(Gm,z,j)).subs(z,s)))
# correction
A = k*z*sp.cosh(k*z) - sp.sinh(k*z); B = k*z*sp.sinh(k*z)
a_, b_ = sp.symbols('a b')
Gtop = ghs(z,s,'z>s')  # at z=1, z>s
eqs = [ (a_*A + b_*B + Gtop).subs(z,1), sp.diff(a_*A + b_*B + Gtop, z).subs(z,1) ]
sol = sp.solve(eqs, [a_, b_], dict=True)[0]
a_sol = sp.simplify(sol[a_]); b_sol = sp.simplify(sol[b_])
print('a =', a_sol); print('b =', b_sol)
import pickle
pickle.dump((a_sol,b_sol), open((ROOT + '/code/slab_corr.pkl'),'wb'))
