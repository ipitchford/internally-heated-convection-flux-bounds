"""Global Nu-maximising rolls for the Rayleigh-Benard scan. The Ra-continuation scan uses a local golden-section
search; where a fixed-Ra Gamma scan shows a higher local maximum elsewhere on the roll family, that maximum
replaces the scan value. Writes steady/rb_ns_best_found.json."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import json, glob, os, re, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steady_shoulder import local_maxima
S = (ROOT + '/data/steady/')
rb = json.load(open(S+'rb_ns.json')); out = []
ref = json.load(open(S+'maxima_1e5.json')) if os.path.exists(S+'maxima_1e5.json') else []
gs = {}
for f in glob.glob(S+'gscan_rb_*.json'):
    d = json.load(open(f)); gs[round(np.log10(d[0]['Ra']), 2)] = f
for r in rb:
    q = dict(Ra=r['Ra'], Gamma=r['Gamma'], Nu=r['Nu'], Nz=r['Nz'], Nx=r['Nx'], source='Ra-continuation scan', local_maxima=[])
    k = round(np.log10(r['Ra']), 2)
    if k in gs:
        lm = local_maxima(gs[k]); q['local_maxima'] = lm
        rr = [x for x in ref if abs(x['Ra']/r['Ra']-1) < 1e-9 and x['Nz'] == max(y['Nz'] for y in ref)]
        if rr:
            best = max(rr, key=lambda x: x['Nu'])
            if best['Nu'] > q['Nu']*(1+1e-6):
                q.update(Gamma=best['Gamma'], Nu=best['Nu'], Nz=best['Nz'], Nx=int(max(32, 8*np.ceil(best['modes_per_unit']*best['Gamma']/8))),
                         source=f"refined maximum ({best['branch']} branch, Nz={best['Nz']})")
            q['local_maxima'] = [dict(Gamma=x['Gamma'], Nu=x['Nu'], branch=x['branch']) for x in rr]
        else:
            for m in lm:
                if m['Nu'] > q['Nu']*(1+1e-5) and abs(np.log(m['Gamma']/q['Gamma'])) > 0.1:
                    q.update(Gamma=m['Gamma'], Nu=m['Nu'], source='Gamma scan')
    out.append(q)
json.dump(out, open(S+'rb_ns_best_found.json', 'w'), indent=1)
for q in out: print('%.3g %.4f %.5f %s %s' % (q['Ra'], q['Gamma'], q['Nu'], q['source'], [(round(m['Gamma'], 3), round(m['Nu'], 5)) for m in q['local_maxima']]))
