import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from flint import arb
from symbols_arb import r_small, r_large
from kernels import corr_symbol
for bc in ['ns', 'sf']:
    d = 2 if bc == 'ns' else 1
    for (z, s) in [(0.0, 0.9), (0.25, 0.5), (0.5, 1.0), (0.1, 0.05)]:
        for k in [0.05, 0.5, 1.5, 3.0, 5.0, 12.0]:
            ref = corr_symbol(np.array([k]), z, s, bc, d)[0] + (s/2*k*np.exp(-k) if bc == 'sf' else 0)
            out = []
            for f in ([r_small] if k <= 6 else []) + ([r_large] if k >= 1.5 else []):
                v = f(arb(k), arb(z), arb(s), bc).coeffs()[0]
                out.append(f'{f.__name__}={float(v.mid()):+.12e} (rad {float(v.rad()):.1e})')
            print(bc, z, s, k, f'ref={ref:+.12e}', ' '.join(out))
