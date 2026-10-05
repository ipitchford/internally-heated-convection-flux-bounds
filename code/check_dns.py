"""Tests of the proof's inequalities on time-averaged simulation statistics, with window uncertainties.

Primary data: the continuation runs (suffix _cont), whose averages are physical-time (trapezoidal)
integrals over every step with all diagnostics taken from one time level, split into 8 equal blocks.
Secondary data: the version-1 windows, re-averaged here with trapezoidal physical-time weights from their
stored scalar histories (their stored profiles are sample averages and are not used)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, glob
from bound_curve import make_Phi, B_of_Lambda
out = []
for jf in sorted(glob.glob(ROOT + '/data/dns/ih_*_cont.json')):
    info = json.load(open(jf))
    d = np.load(jf.replace('.json', '_prof.npz')); z, Tb, J, W = d['z'], d['Tbar'], d['J'], d['Wmax']
    Phi, a, b = make_Phi(info['bc'], 1)
    R, Ts, tau = info['R'], info['Tmax'], info['window']
    m = (z > 1e-3) & (z < 1-1e-3); lo = m & (z < 0.2)
    r1 = (J/np.maximum(W*Tb, 1e-300))[m]; r2 = (W/(0.5*R*Ts*Phi(z)))[m]
    blocks = np.array([[B['FB'], B['FT'], B['wT']] for B in info['blocks']])
    nb = len(blocks); se = blocks.std(0, ddof=1)/np.sqrt(nb)
    Bw = float(B_of_Lambda(Phi, [0.5*R*Ts])[0])
    # the version-1 window of the same run, re-averaged with physical-time (trapezoidal) weights
    v1 = jf.replace('_cont.json', '.json'); v1rec = None
    if os.path.exists(v1):
        i1 = json.load(open(v1)); h = np.load(v1.replace('.json', '_hist.npz'))['hist']
        t = h[:, 0]; w = (t >= i1['t_avg'][0]-1e-12) & (t <= i1['t_avg'][1]+1e-12)
        tw, L = t[w], t[w][-1]-t[w][0]
        tr = lambda col: float(np.trapezoid(h[w, col], tw)/L)
        v1rec = dict(window=[float(tw[0]), float(tw[-1])], FB_trapezoid=tr(1), FT_trapezoid=tr(2), FB_sample_mean=i1['FB'],
                     FT_sample_mean=i1['FT'], dTm_dt=float((h[w, 3][-1]-h[w, 3][0])/L), imbalance=1-tr(1)-tr(2))
    # character of the flow over the window: time spent at a steady plateau, and number of excursions from it
    hh = np.load(jf.replace('.json', '_hist.npz'))['hist']; th, fb = hh[:, 0], hh[:, 1]
    ww = (th >= info['t_avg'][0]-1e-12) & (th <= info['t_avg'][1]+1e-12); th, fb = th[ww], fb[ww]
    dtw = np.gradient(th); vals, inv = np.unique(np.round(fb, 5), return_inverse=True)
    plateau = float(vals[np.argmax(np.bincount(inv, weights=dtw))])
    at = np.abs(fb - plateau) < 1e-4; frac = float(np.sum(dtw[at])/np.sum(dtw))
    # bursts: departures by more than 1e-3, merging departures less than 0.2 time units apart
    far = np.abs(fb - plateau) > 1e-3; starts = th[1:][np.diff(far.astype(int)) == 1]
    if far[0]: starts = np.concatenate([[th[0]], starts])
    exc = int(len(starts) - np.sum(np.diff(starts) < 0.2)) if len(starts) else 0
    character = dict(plateau_FB=plateau if frac > 0.2 else None, time_fraction_at_plateau=frac,
                     excursions=exc if frac > 0.2 else None, FB_min=float(fb.min()), FB_max=float(fb.max()))
    rec = dict(name=info['name'], R=R, bc=info['bc'], window=info['t_avg'], tau=tau, nsteps=info['nsteps'], character=character,
               FB=info['FB'], FT=info['FT'], wT=info['wT'], Tmean=info['Tmean'],
               FB_block_se=float(se[0]), FT_block_se=float(se[1]), FB_block_range=[float(blocks[:, 0].min()), float(blocks[:, 0].max())],
               imbalance=info['imbalance'], dTm_dt=info['dTm_dt'], identity_residual=info['imbalance']-info['dTm_dt'],
               Tsup=Ts, Tmin=info['Tmin'], max_J_over_WTheta=float(r1.max()), max_W_over_bound=float(r2.max()),
               max_W_over_bound_z_lt_0p2=float((W/(0.5*R*Ts*Phi(z)))[lo].max()), min_J_z_lt_0p2=float(J[lo].min()),
               finite_window_bound=Bw - Ts/tau, bound_window_T=Bw, bound_T_eq_1_8=float(B_of_Lambda(Phi, [0.5*R*0.125])[0]),
               v1_window=v1rec)
    out.append(rec); print(json.dumps(rec))
json.dump(out, open(ROOT + '/data/dns/dns_checks.json', 'w'), indent=1)
