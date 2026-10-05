"""Spectral tails of every stored steady solution: largest coefficient in the last 10% of Fourier modes and of
Chebyshev modes, relative to the largest coefficient. X has shape (Nz+1, K): Chebyshev-Lobatto values in z,
cosine coefficients in x."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, glob, json, re
from scipy.fft import dct
out = []
for f in sorted(glob.glob((ROOT + '/data/steady/[ri][bh]_ns_X_*.npy'))):
    X = np.load(f); N = X.shape[0]-1; K = X.shape[1]
    C = dct(X, type=1, axis=0)/N; C[0] /= 2; C[-1] /= 2        # Chebyshev coefficients along z
    mx = np.abs(C[:, 1:]).max()
    ft = np.abs(C[:, int(0.9*K):]).max()/mx; ct = np.abs(C[int(0.9*N):, :]).max()/mx
    Ra = float(re.search(r'X_([0-9.e+]+)\.npy', f).group(1)); mode = 'rb' if '/rb_' in f else 'ih'
    out.append(dict(mode=mode, Ra=Ra, Nz=N, K=K, fourier_tail=float(ft), cheb_tail=float(ct)))
    print(f'{mode} Ra={Ra:.3e} Nz={N} K={K} fourier_tail={ft:.1e} cheb_tail={ct:.1e}')
json.dump(out, open((ROOT + '/data/steady/spectral_tails.json'), 'w'), indent=1)
m = [o for o in out if o['Ra'] >= 1e6]
print('max over Ra>=1e6: fourier %.1e, cheb %.1e' % (max(o['fourier_tail'] for o in m), max(o['cheb_tail'] for o in m)))
