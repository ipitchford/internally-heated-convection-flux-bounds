"""Second route to the slab-correction constants: the Carlson-type majorants of Lemma 3.6(iii),
    C_ns <= sup_Q C[mu_ns],   C_sf <= nu_1/2 + sup_Q C[r_sf],
evaluated on the same (zeta, s) grid as corr_constants2.py, with two levels of refinement around the maximiser.
The functionals C[phi] involve only integrals of squares of the symbol and its first two derivatives (no
oscillatory inverse transform), so they are less exposed to truncation error than the direct L1 norms.
Derivatives are taken by second-order finite differences on a uniform grid (h = 5e-4 on [1e-6, 120]); the
result is checked against h = 2.5e-4 at the maximiser. Not interval-certified.
Usage: python carlson_sup.py  -> code/carlson_sup.json"""
import numpy as np, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import corr_symbol

NU1 = {1: 2/np.pi, 2: 4/(3*np.sqrt(3))}


def symbol_fast(k, zeta, s, bc, d, deg=48):
    """phi on a fine grid: on kappa < 1 the symbol is evaluated in 50-digit arithmetic at deg+1 Chebyshev
    points of [0, 1] and interpolated (it is analytic there, Lemma 3.6(i)); on kappa >= 1 it is evaluated
    directly in double precision."""
    out = np.empty_like(k); lo = k < 1.0
    x = 0.5*(1 - np.cos(np.pi*(np.arange(deg+1) + 0.5)/(deg+1)))
    fx = corr_symbol(x, zeta, s, bc, d)
    cheb = np.polynomial.chebyshev.Chebyshev.fit(x, fx, deg, domain=[0, 1])
    out[lo] = cheb(k[lo]); out[~lo] = corr_symbol(k[~lo], zeta, s, bc, d)
    return out


def carlson(bc, n, zeta, s, h=5e-4, kmax=120.0):
    d = 2 if bc == 'ns' else 1
    k = np.arange(1e-6, kmax, h)
    f = symbol_fast(k, zeta, s, bc, d)
    if bc == 'sf':
        f = f + (s/2)*k*np.exp(-k)          # r_sf = mu_sf + (s/2) kappa e^{-kappa}
    f1 = np.gradient(f, h, edge_order=2); f2 = np.gradient(f1, h, edge_order=2)
    if n == 1:
        I0, I1 = np.trapezoid(f*f, k), np.trapezoid(f1*f1, k)
        return float(np.sqrt(2)*(I0*I1)**0.25)
    J0 = np.trapezoid(f*f*k, k); J1 = np.trapezoid(f1*f1*k, k); J2 = np.trapezoid((f2 + f1/k)**2*k, k)
    return float(np.sqrt(np.sqrt(J0*J2) + J1))


def sup_over_Q(bc, n):
    grid = {}
    def ev(z, s):
        key = (round(z, 6), round(s, 6))
        if key not in grid:
            grid[key] = carlson(bc, n, z, s)
        return grid[key]
    zs, ss = np.linspace(0, 0.5, 11), np.linspace(0, 1, 21)
    for z in zs:
        for s in ss:
            ev(z, s)
    hz, hs = 0.05, 0.05
    for level in range(2):
        (z0, s0) = max(grid, key=grid.get)
        hz, hs = hz/4, hs/4
        for z in np.clip(z0 + hz*np.arange(-4, 5), 0, 0.5):
            for s in np.clip(s0 + hs*np.arange(-4, 5), 0, 1):
                ev(z, s)
    (z0, s0) = max(grid, key=grid.get)
    return grid[(z0, s0)], (z0, s0), len(grid)


if __name__ == '__main__':
    out = {}
    for n in (1, 2):
        for bc in ('ns', 'sf'):
            t = time.time()
            sup, arg, ne = sup_over_Q(bc, n)
            check = carlson(bc, n, arg[0], arg[1], h=2.5e-4)
            C = sup if bc == 'ns' else NU1[n]/2 + sup
            out[f'{bc}_n{n}'] = dict(sup_carlson=sup, argmax=arg, n_evals=ne, half_step_check=check,
                                     C_majorant=C, seconds=time.time()-t)
            print(f'{bc} n={n}: sup C[.] = {sup:.6f} at {arg} (h/2: {check:.6f}); C_majorant = {C:.6f}; '
                  f'{ne} evaluations', flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'carlson_sup.json'), 'w'),
              indent=1)
