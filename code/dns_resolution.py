"""Spectral resolution of the simulations: for the stored snapshots and the final state of each continuation run,
the largest Fourier coefficient of T in the last 10% of the retained modes and the largest Chebyshev coefficient
in the last 10% of the Chebyshev modes, each relative to the largest coefficient. Writes data/dns/resolution.json."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, glob, json
from scipy.fft import dct


def tails(Tphys_or_hat, spectral=False, Nk=None):
    """Tphys: (Nz+1, M) physical values on Chebyshev-Lobatto x uniform grid; or Th: (Nz+1, Nk) Fourier coefficients."""
    Th = Tphys_or_hat if spectral else np.fft.rfft(Tphys_or_hat, axis=1)/Tphys_or_hat.shape[1]
    if Nk is not None: Th = Th[:, :Nk]
    N = Th.shape[0]-1; K = Th.shape[1]
    C = dct(Th.real, type=1, axis=0)/N + 1j*dct(Th.imag, type=1, axis=0)/N; C[0] /= 2; C[-1] /= 2
    mx = np.abs(C[:, 1:]).max() if K > 1 else np.abs(C).max()
    return float(np.abs(C[:, int(0.9*K):]).max()/mx), float(np.abs(C[int(0.9*N):, :]).max()/mx)


out = []
for f in sorted(glob.glob(ROOT + '/data/dns/ih_*_cont.json')):
    info = json.load(open(f)); st = np.load(f.replace('.json', '_state.npz')); sn = np.load(f.replace('.json', '_snap.npz'))
    Nk = st['Th'].shape[1]
    ft, ct = tails(st['Th'], spectral=True)
    rec = dict(name=info['name'], final_state=dict(t=float(st['t']), fourier_tail=ft, cheb_tail=ct), snapshots=[])
    for i, t in enumerate(sn['t']):
        a, b = tails(sn['T'][i], Nk=Nk)
        rec['snapshots'].append(dict(t=float(t), fourier_tail=a, cheb_tail=b))
    rec['max_fourier_tail'] = max([ft] + [s['fourier_tail'] for s in rec['snapshots']])
    rec['max_cheb_tail'] = max([ct] + [s['cheb_tail'] for s in rec['snapshots']])
    out.append(rec)
    print(f"{info['name']}: max Fourier tail {rec['max_fourier_tail']:.1e}, max Chebyshev tail {rec['max_cheb_tail']:.1e} "
          f"(final state {ft:.1e}, {ct:.1e}; {len(sn['t'])} snapshots at t = {', '.join(f'{t:.2f}' for t in sn['t'])})")
json.dump(out, open(ROOT + '/data/dns/resolution.json', 'w'), indent=1)
