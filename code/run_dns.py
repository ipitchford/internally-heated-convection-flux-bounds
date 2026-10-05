"""Driver: run one IH (or RB) infinite-Pr DNS case and save statistics needed to test the bound."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import sys, json, numpy as np, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipc2d import IPC2D
import argparse
p = argparse.ArgumentParser()
p.add_argument('--R', type=float); p.add_argument('--bc', default='ns'); p.add_argument('--G', type=float, default=4.0)
p.add_argument('--Nx', type=int, default=256); p.add_argument('--Nz', type=int, default=64)
p.add_argument('--tend', type=float); p.add_argument('--tavg', type=float); p.add_argument('--rb', action='store_true')
p.add_argument('--wall', type=float, default=None); p.add_argument('--tag', default='')
p.add_argument('--restart', default=None); p.add_argument('--blocks', type=int, default=8)
a = p.parse_args()
name = f"{'rb' if a.rb else 'ih'}_{a.bc}_R{a.R:.0e}_G{a.G:g}_{a.Nx}x{a.Nz}{a.tag}"
s = IPC2D(R=a.R, Gamma=a.G, Nx=a.Nx, Nz=a.Nz, bc=a.bc, Q=0.0 if a.rb else 1.0, rb=a.rb)
if a.restart:
    d = np.load(a.restart); s.Th = d['Th']; s.t = float(d['t'])
t0 = time.time()
hist, acc, snaps = s.run(t_end=a.tend, t_avg_start=a.tavg, diag_every=1, log_every=20000, max_wall=a.wall,
                         snapshot_every=20000, n_blocks=a.blocks)
out = dict(name=name, R=a.R, bc=a.bc, G=a.G, Nx=a.Nx, Nz=a.Nz, rb=a.rb, wall=time.time()-t0, t_final=s.t,
           restart=a.restart, averaging='physical-time trapezoidal; diagnostics from a single time level')
if acc is not None:
    Bs = acc['blocks']; W = sum(B['W'] for B in Bs)
    tot = {k: sum(B[k] for B in Bs)/W for k in ('FB', 'FT', 'wT', 'Tm')}
    blk = [dict(ta=B['ta'], tb=B['tb'], W=B['W'], **{k: B[k]/B['W'] for k in ('FB', 'FT', 'wT', 'Tm')}) for B in Bs if B['W'] > 0]
    out.update(FB=tot['FB'], FT=tot['FT'], wT=tot['wT'], Tmean=tot['Tm'], Tmax=acc['Tmax'], Tmin=acc['Tmin'],
               t_avg=[acc['t0'], acc['t1']], window=W, nsteps=acc['nsteps'], Tm_start=acc['Tm0'], Tm_end=acc['Tm1'],
               imbalance=1-tot['FB']-tot['FT'], dTm_dt=(acc['Tm1']-acc['Tm0'])/W, blocks=blk)
    np.savez(f'{ROOT}/data/dns/{name}_prof.npz', z=s.z, Tbar=sum(B['Tbar'] for B in Bs)/W, J=sum(B['J'] for B in Bs)/W,
             Wmax=acc['Wmax'], Tbar_blocks=np.array([B['Tbar']/B['W'] for B in Bs if B['W'] > 0]),
             J_blocks=np.array([B['J']/B['W'] for B in Bs if B['W'] > 0]))
np.savez(f'{ROOT}/data/dns/{name}_hist.npz', hist=hist)
np.savez(f'{ROOT}/data/dns/{name}_state.npz', Th=s.Th, t=s.t)
if snaps:
    fields_T = np.array([x[1] for x in snaps])
    fields_w = np.array([x[2] for x in snaps])
    x_padded = np.arange(fields_T.shape[-1]) * a.G / fields_T.shape[-1]
    if fields_T.shape != fields_w.shape or fields_T.shape[-2] != len(s.z):
        raise ValueError('Snapshot field/coordinate shapes disagree')
    np.savez_compressed(f'{ROOT}/data/dns/{name}_snap.npz', t=np.array([x[0] for x in snaps]),
                        T=fields_T, w=fields_w, x=x_padded, x_base=s.x, z=s.z)
json.dump(out, open(f'{ROOT}/data/dns/{name}.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
