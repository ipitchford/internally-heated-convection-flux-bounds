"""Probe: Carlson-type upper bounds for the slab-correction kernels at selected (zeta, s) (not certified)."""
import numpy as np, mpmath as mp, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import corr_symbol, L1_corr_deriv
def sym(bc, d, zeta, s):
    cusp = -s/2 if bc == 'sf' else 0.0
    return lambda k: corr_symbol(np.atleast_1d(k), zeta, s, bc, d) - cusp*np.atleast_1d(k)*np.exp(-np.atleast_1d(k))
def carlson(bc, n, zeta, s):
    d = 2 if bc == 'ns' else 1
    k = np.linspace(1e-6, 120, 240001); h = k[1]-k[0]
    f = sym(bc, d, zeta, s)(k)
    f1 = np.gradient(f, h); f2 = np.gradient(f1, h)
    if n == 1:
        I0, I1 = np.trapezoid(f*f, k), np.trapezoid(f1*f1, k)
        return np.sqrt(2)*(I0*I1)**0.25
    J0 = np.trapezoid(f*f*k, k); J1 = np.trapezoid(f1*f1*k, k); J2 = np.trapezoid((f2 + f1/k)**2*k, k)
    return np.sqrt(np.sqrt(J0*J2) + J1)
for bc, n, z, s in [('ns', 1, 0.0, 0.9125), ('ns', 1, 0.25, 0.5), ('ns', 2, 0.0, 0.9), ('sf', 1, 0.0, 0.8), ('sf', 1, 0.3, 0.5), ('sf', 2, 0.0, 0.8)]:
    d = 2 if bc == 'ns' else 1
    true = L1_corr_deriv(z, s, n, bc, d)
    cusp = (s/2)*(2/np.pi if n == 1 else 4/(3*np.sqrt(3))) if bc == 'sf' else 0.0
    print(f'{bc} n={n} (zeta,s)=({z},{s}): computed L1 = {true:.5f}; cusp part = {cusp:.5f}; Carlson bound on remainder = {carlson(bc, n, z, s):.5f}')
