"""Regression checks for the October 2026 referee corrections. Works under -O."""
import json
from pathlib import Path
import numpy as np
from profile_extrema import inspect_profile

root = Path(__file__).resolve().parents[1]
results = []
for path in sorted((root/'data/dns').glob('*_cont.json')):
    info = json.loads(path.read_text())
    with np.load(path.with_name(path.stem+'_prof.npz')) as d:
        record = inspect_profile(d['z'], d['Tbar'], info['FB'], info['FT'])
    record['case'] = path.stem
    results.append(record)
bad = next(x for x in results if x['case'] == 'ih_ns_R1e+06_G4_256x64_cont')
if not (bad['grid_pass'] and not bad['slope_compatible'] and not bad['interpolant_pass']):
    raise RuntimeError('Known off-grid cooling violation not detected')
if abs(bad['top_limit_ratio'] - 1.0001190511321745) > 1e-8:
    raise RuntimeError('Unexpected boundary-slope regression')
z = (1+np.cos(np.arange(65)*np.pi/64))/2
good = inspect_profile(z, z*(1-z)/2, .5, .5)
bad_control = inspect_profile(z, .501*z*(1-z), .5, .5)
if not good['interpolant_pass'] or bad_control['slope_compatible'] or bad_control['interpolant_pass']:
    raise RuntimeError('Positive/negative profile control failed')
snapshots = sorted((root/'data/dns').glob('*_snap.npz'))
if len(snapshots) != 8:
    raise RuntimeError('Expected all eight snapshots')
for path in snapshots:
    with np.load(path) as d:
        if d['T'].shape != d['w'].shape or d['T'].shape[-2:] != (len(d['z']),len(d['x'])):
            raise RuntimeError('Snapshot coordinate mismatch: '+path.name)
print(json.dumps({'profiles': results, 'snapshots_checked':len(snapshots),
                  'controls':'conductive pass; perturbed profile rejected'}, indent=2))
