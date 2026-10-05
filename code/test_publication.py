"""Fail-closed supporting replay and semantic controls; no continuum proof claim."""
from pathlib import Path
import subprocess,sys,json,hashlib,tempfile,shutil
import numpy as np
r=Path(__file__).resolve().parents[1]
for script in ['revision_checks.py','verify_numbers.py']:
 for flags in [[],['-O']]:
  subprocess.run([sys.executable,*flags,str(r/'code'/script)],cwd=r/'code',check=True)
from profile_extrema import inspect_profile
z=(1+np.cos(np.arange(65)*np.pi/64))/2
for amplitude in [.501,.51]:
 d=inspect_profile(z,amplitude*z*(1-z),.5,.5)
 if d['interpolant_pass'] or d['slope_compatible']:raise RuntimeError('Corrupt profile accepted')
with tempfile.TemporaryDirectory(prefix='ihc-negative-') as tmp:
 root=Path(tmp);shutil.copytree(r/'code',root/'code',ignore=shutil.ignore_patterns('__pycache__'))
 shutil.copytree(r/'data/dns',root/'data/dns')
 p=next((root/'data/dns').glob('*_snap.npz'))
 with np.load(p) as d:data={k:d[k] for k in d.files}
 data['x']=data['x'][:-1];np.savez_compressed(p,**data)
 result=subprocess.run([sys.executable,'-O','revision_checks.py'],cwd=root/'code',capture_output=True,text=True)
 if result.returncode==0 or 'Snapshot coordinate mismatch' not in result.stderr:raise RuntimeError('Corrupt snapshot not rejected for intended reason')
print('PUBLICATION SUPPORTING REPLAY PASSED: normal, optimized, two profile controls, snapshot corruption')
