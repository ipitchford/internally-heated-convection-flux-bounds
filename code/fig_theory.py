"""Figures 1-2: wall functions and lower-bound curves."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from bound_curve import make_Phi, B_of_Lambda
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'axes.spines.top': False, 'axes.spines.right': False})
nu2 = {1: 3*np.sqrt(3)/(2*np.pi), 2: 1.115419}
have2sf = os.path.exists((ROOT + '/code/phi_sf_n2_26.json'))
# ---------- Figure 1: wall functions
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
z = np.geomspace(1e-4, 0.5, 300)
for n, ls in [(1, '-'), (2, '--')]:
    Phi, a, b = make_Phi('ns', n)
    d = json.load(open(f'{ROOT}/code/phi_ns_n{n}_26.json')); zz = np.array(sorted(map(float, d))); P = np.array([d[k] for k in sorted(d, key=float)])
    ax[0].plot(np.log(1/zz), P/zz**2, 'o', ms=2.5, color='C0' if n == 1 else 'C1')
    ax[0].plot(np.log(1/z), a*np.log(1/z) + b, ls, lw=0.8, color='C0' if n == 1 else 'C1',
               label=f'$d={n+1}$: $\\frac{{\\nu_2}}{{2}}\\ln(1/z){b:+.2f}$')
ax[0].set_xlabel('$\\ln(1/z)$'); ax[0].set_ylabel('$\\Phi_{\\rm ns}(z)/z^2$'); ax[0].legend(frameon=False, fontsize=7)
ax[0].set_title('(a) no-slip', fontsize=9)
for n in [1, 2]:
    if n == 2 and not have2sf: continue
    d = json.load(open(f'{ROOT}/code/phi_sf_n{n}_26.json')); zz = np.array(sorted(map(float, d))); P = np.array([d[k] for k in sorted(d, key=float)])
    ax[1].semilogx(zz, P/zz, 'o-', ms=2.5, lw=0.8, label=f'$d={n+1}$')
ax[1].set_xlabel('$z$'); ax[1].set_ylabel('$\\Phi_{\\rm sf}(z)/z$'); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title('(b) stress-free', fontsize=9)
fig.tight_layout(); fig.savefig((ROOT + '/figs/fig1_wall_functions.pdf')); fig.savefig((ROOT + '/figs/fig1_wall_functions.png'), dpi=200)
# ---------- Figure 2: bounds vs R (d=2), with Arslan & Rojas (2025) and simulations
R = np.geomspace(1e4, 1e12, 200)
fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
AR = {'ns': lambda R: 0.1555*R**(-2/3) - 0.7104*R**(-0.5)*np.log(1-1.5334*R**(-1/3)),
      'sf': lambda R: 0.0131*R**(-40/29) - 0.0312*R**(-35/29)*np.log(1-0.6301*R**(-10/29))}
ARdata = {'ns': [(5e4,.0782),(1e5,.0571),(5e5,.162),(1e6,.133),(5e6,.183),(1e7,.213),(5e7,.242),(1e8,.268),(5e8,.275),(1e9,.313),(5e9,.329)],
          'sf': [(5e4,.0782),(1e5,.108),(5e5,.161),(1e6,.180),(5e6,.241),(1e7,.257),(5e7,.299),(1e8,.312),(5e8,.345),(1e9,.356)]}
checks = json.load(open((ROOT + '/data/dns/dns_checks.json'))) if os.path.exists((ROOT + '/data/dns/dns_checks.json')) else []
for k, bc in enumerate(['ns', 'sf']):
    Phi, a, b = make_Phi(bc, 1)
    Tmax = json.load(open(ROOT + '/paper/facts_constants.json'))['Tmax_d2_period4']
    ax[k].loglog(R, B_of_Lambda(Phi, 0.5*R*Tmax), color='C3', lw=1.6,
                 label=f'Corollary 2.4, $\\mathcal{{T}}_{{\\max}}={Tmax:.3f}$ (estimate)')
    for Ts, col, ls in [(1.0, 'C2', '-.'), (1/8, 'C0', ':')]:
        ax[k].loglog(R, B_of_Lambda(Phi, 0.5*R*Ts), color=col, ls=ls, lw=1.2,
                     label=f'Theorem 2.1 if $\\mathcal{{T}}\\leq{"1/8" if Ts < 1 else "1"}$ (conditional)')
    ax[k].loglog(R, AR[bc](R), 'k--', lw=1, label='Arslan & Rojas (2025)')
    Rd, wT = np.array(ARdata[bc]).T; ax[k].loglog(Rd, 0.5-wT, 'ks', ms=3, mfc='none', label='DNS, Arslan & Rojas')
    ours = [c for c in checks if c['bc'] == bc]
    if ours:
        ax[k].loglog([c['R'] for c in ours], [c['FB'] for c in ours], '^', color='C4', ms=4, label='DNS, this work')
    ax[k].set_xlabel('$R$'); ax[k].set_ylabel('$\\overline{F_B}$'); ax[k].set_ylim(1e-9, 1)
    ax[k].set_title('(a) no-slip' if bc == 'ns' else '(b) stress-free', fontsize=9)
ax[0].legend(frameon=False, fontsize=6.5, loc='lower left')
fig.tight_layout(); fig.savefig((ROOT + '/figs/fig2_bounds.pdf')); fig.savefig((ROOT + '/figs/fig2_bounds.png'), dpi=200)
print('ok')
