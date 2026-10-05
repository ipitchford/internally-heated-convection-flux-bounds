"""Test the profile hypothesis (14) of Chen et al. (2024): <T>_h(z) <= min(F_B z, (1-z)(1-F_B)),
on time-averaged DNS profiles and on steady states (where it must hold pointwise in time)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, glob, sys
from profile_extrema import inspect_profile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
res = []
def test(z, Theta, FB, label):
    tri = np.minimum(FB*z, (1-z)*(1-FB)); m = (z > 1e-4) & (z < 1-1e-4)
    r = Theta[m]/tri[m]; i = np.argmax(r)
    rec = dict(case=label, FB=float(FB), max_ratio=float(r[i]), z_at_max=float(z[m][i]),
               **inspect_profile(z, Theta, FB))
    res.append(rec); print(json.dumps(rec))
for jf in sorted(glob.glob(ROOT + '/data/dns/ih_*_cont.json')):
    info = json.load(open(jf))
    d = np.load(jf.replace('.json', '_prof.npz'))
    test(d['z'], d['Tbar'], info['FB'], 'DNS ' + info['name'])
    # block-to-block spread of the ratio (each block profile against its own block-mean F_B)
    rb = []
    for Tbk, B in zip(d['Tbar_blocks'], info['blocks']):
        z = d['z']; tri = np.minimum(B['FB']*z, (1-z)*(1-B['FB'])); m = (z > 1e-4) & (z < 1-1e-4)
        rb.append(float((Tbk[m]/tri[m]).max()))
    res[-1]['block_max_ratios'] = rb
    print('   blocks:', ' '.join(f'{x:.5f}' for x in rb))
from steady2d import Steady
for jf in sorted(glob.glob((ROOT + '/data/steady/ih_*.json'))):
    recs = json.load(open(jf))
    for r in recs:
        X = np.load(jf.replace('.json', f"_X_{r['Ra']:.3e}.npy"))
        Nz = X.shape[0]-1; z = (1+np.cos(np.pi*np.arange(Nz+1)/Nz))/2
        test(z, X[:, 0], r['FB'], f"steady {jf.split('/')[-1][:-5]} R={r['Ra']:.3g} G={r['Gamma']:.3f}")
json.dump(res, open((ROOT + '/data/dns/chen_hypothesis_checks.json'), 'w'), indent=1)
