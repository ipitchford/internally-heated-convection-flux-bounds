import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb
from symbols_arb import r_small, r_large
def ball(a, b): return arb((a+b)/2, (b-a)/2)
for bc in ['ns', 'sf']:
    for k in [0.5, 1.5, 2.9, 6.0]:
        f = r_small if k < 3 else r_large
        for lab, K, Z, S in [('point k, box 0.005', arb(k), ball(0, 0.005), ball(0.9, 0.905)),
                             ('k-ball 0.02, point params', ball(k, k+0.02), arb(0.0025), arb(0.9025)),
                             ('k-ball 0.002, point params', ball(k, k+0.002), arb(0.0025), arb(0.9025))]:
            c = f(K, Z, S, bc).coeffs()
            print(bc, k, lab.ljust(28), ' '.join(f'[{float(x.mid()):+.4e} +/- {float(x.rad()):.1e}]' for x in c))
