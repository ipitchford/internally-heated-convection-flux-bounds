import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, sys, json
from kernels import L1_corr_deriv
n = int(sys.argv[1]); out = {}
for bc, d in [('ns', 2), ('sf', 1)]:
    best = (0, None)
    vals = {}
    for zeta in np.linspace(0, 0.5, 6):
        for s in np.linspace(0, 1, 11):
            v = L1_corr_deriv(zeta, s, n, bc, d); vals[f'{zeta:.2f},{s:.2f}'] = v
            if v > best[0]: best = (v, (zeta, s))
    print(n, bc, 'd=%d'%d, 'sup ||F^-1[k^2 d^%d c]||_1 = %.5f at (zeta,s)=%s'%(d, best[0], best[1]), flush=True)
    out[bc] = {'sup': best[0], 'argmax': best[1], 'grid': vals}
json.dump(out, open(f'{ROOT}/code/corr_constants_n{n}.json','w'), indent=1)
