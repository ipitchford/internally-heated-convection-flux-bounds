import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from bound_curve import make_Phi
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'axes.spines.top': False, 'axes.spines.right': False})
name = sys.argv[1] if len(sys.argv) > 1 else 'ih_ns_R1e+06_G4_256x64_cont'
info = json.load(open(f'{ROOT}/data/dns/{name}.json')); d = np.load(f'{ROOT}/data/dns/{name}_prof.npz')
z, Tb, J, W = d['z'], d['Tbar'], d['J'], d['Wmax']; o = np.argsort(z); z, Tb, J, W = z[o], Tb[o], J[o], W[o]
Phi, a, b = make_Phi(info['bc'], 1); R, Ts = info['R'], info['Tmax']
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
ax[0].plot(Tb, z, 'k', label=r'$\Theta$'); FB = info['FB']
ax[0].plot(np.minimum(FB*z, (1-z)*(1-FB)), z, 'C3--', lw=0.8, label='Chen et al. triangle')
ax[0].set_xlabel('temperature'); ax[0].set_ylabel('$z$'); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title('(a) mean profile', fontsize=9)
ax[1].plot(J, z, 'C0', label=r'$J=\langle wT\rangle_h$'); ax[1].plot(W*Tb, z, 'C1', label=r'$W\Theta$'); ax[1].plot(z - FB, z, 'C2:', lw=0.8, label=r'$z-F_B$')
ax[1].axvline(0, color='0.6', lw=0.5); ax[1].set_xlabel('flux'); ax[1].set_xlim(-0.05, 0.8); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title('(b) fluxes', fontsize=9)
ax[2].semilogx(W, z, 'C0', label=r'$W=\max_{x,t}|w|$'); ax[2].semilogx(0.5*R*Ts*Phi(z), z, 'k--', label=r"$R\mathcal{T}\Phi(z)/2$")
ax[2].set_xlabel('vertical speed'); ax[2].set_xlim(1e-2, None); ax[2].legend(frameon=False, fontsize=7); ax[2].set_title('(c) velocity bound', fontsize=9)
fig.suptitle(f"{'no-slip' if info['bc']=='ns' else 'stress-free'}, $R=10^{{{int(round(np.log10(R)))}}}$, time-averaged", fontsize=9)
fig.tight_layout(); fig.savefig(f'{ROOT}/figs/fig3_{name}.pdf'); fig.savefig(f'{ROOT}/figs/fig3_{name}.png', dpi=200)
print('ok', J[(z>0.005)&(z<0.3)].min(), J[(z>0.005)&(z<0.3)].max())
