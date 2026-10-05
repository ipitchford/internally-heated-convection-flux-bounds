import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb
from symbols_arb import r_small, r_large
def ball(a, b): return arb((a+b)/2, (b-a)/2)
for bc in ['ns', 'sf']:
    for (k0, k1) in [(0.0, 0.02), (0.5, 0.52), (1.5, 1.52), (2.9, 2.92), (3.0, 3.02), (6.0, 6.05), (20, 20.2), (60, 61)]:
        for (h) in [0.02, 0.005]:
            Z, S = ball(0.0, h), ball(0.9, 0.9+h)
            t = time.time()
            f = r_small if k1 <= 3.05 else r_large
            c = f(ball(k0, k1), Z, S, bc).coeffs()
            dt = time.time()-t
            print(bc, f.__name__, (k0, k1), 'box', h, ' '.join(f'[{float(x.mid()):+.4e} +/- {float(x.rad()):.1e}]' for x in c), f'{dt*1e3:.0f} ms')
