import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.ticker; import matplotlib.pyplot as plt
from bound_curve import make_Phi, B_of_Lambda
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'axes.spines.top': False, 'axes.spines.right': False})
ih = json.load(open((ROOT + '/data/steady/ih_ns.json'))); rb = json.load(open((ROOT + '/data/steady/rb_ns_best_found.json')))
fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.4)); ax = axs.ravel()
# (a) IH
R = np.array([r['Ra'] for r in ih]); FB = np.array([r['FB'] for r in ih]); Tm = np.array([r['Tmax'] for r in ih])
ARd = np.array([(5e4,.0782),(1e5,.0571),(5e5,.162),(1e6,.133),(5e6,.183),(1e7,.213),(5e7,.242),(1e8,.268),(5e8,.275),(1e9,.313),(5e9,.329)])
ax[0].loglog(R, FB, 'o-', ms=3, lw=0.8, label='steady, $F_B$-minimising')
ax[0].loglog(ARd[:, 0], 0.5-ARd[:, 1], 'ks', ms=3, mfc='none', label='DNS (Arslan & Rojas)')
Phi, a, b = make_Phi('ns', 1); Rg = np.geomspace(5e4, 5e9, 60)
ax[0].loglog(Rg, B_of_Lambda(Phi, Rg/16), 'C2', lw=1, label=r'bound, $\mathcal{T}=1/8$')
ax[0].loglog(R, B_of_Lambda(Phi, 0.5*R*Tm), 'C3--', lw=1, label=r'bound, steady $\mathcal{T}$')
ax[0].set_xlabel('$R$'); ax[0].set_ylabel('$F_B$'); ax[0].legend(frameon=False, fontsize=6); ax[0].set_title('(a) internally heated', fontsize=9)
# (b) RB compensated
Ra = np.array([r['Ra'] for r in rb]); Nu = np.array([r['Nu'] for r in rb])
wgd_glob = np.array([(1e6, 10.39171), (10**6.5, 14.62274), (1e7, 20.64616), (10**7.5, 29.23512), (1e8, 41.48855), (10**8.5, 58.99612), (1e9, 84.40976)])
wgd_G2 = np.array([(1e6, 8.896615), (10**6.5, 11.71065), (1e7, 15.34493), (10**7.5, 20.47798), (1e8, 27.50669), (1e9, 48.94284)])
ax[1].semilogx(Ra, Nu/Ra**(1/3), 'o-', ms=3, lw=0.8, label=r'Pr$=\infty$, Nu-maximising (this work)')
ax[1].semilogx(wgd_glob[:, 0], wgd_glob[:, 1]/wgd_glob[:, 0]**(1/3), 'k^', ms=3, mfc='none', label='Pr$=1$, global max (WGD 2022)')
ax[1].semilogx(wgd_G2[:, 0], wgd_G2[:, 1]/wgd_G2[:, 0]**(1/3), 'kv', ms=3, mfc='none', label=r'Pr$=1$, $\Gamma=2$ (WGD 2022)')
ax[1].set_xlabel('Ra'); ax[1].set_ylabel(r'Nu$\,$Ra$^{-1/3}$'); ax[1].legend(frameon=False, fontsize=6); ax[1].set_title('(b) Rayleigh–Bénard', fontsize=9)
# (c) Nu(Gamma)
SC = [('gscan_rb_1e5.json', 1e5, r'$10^5$'), ('gscan_rb_1.78e5.json', 10**5.25, r'$1.8\times10^5$'), ('gscan_rb_3.16e5.json', 10**5.5, r'$3.2\times10^5$'),
      ('gscan_rb_5.62e5.json', 10**5.75, r'$5.6\times10^5$'), ('gscan_rb_1e6.json', 1e6, r'$10^6$'), ('gscan_rb_1e7.json', 1e7, r'$10^7$')]
SC = [x for x in SC if os.path.exists((ROOT + '/data/steady/') + x[0])]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(SC)))
for (f, Rv, lab), col in zip(SC, cols):
    d = sorted(json.load(open((ROOT + '/data/steady/') + f)), key=lambda x: x['Gamma']); G = np.array([x['Gamma'] for x in d]); N = np.array([x['Nu'] for x in d])
    ax[2].semilogx(G, N/N.max(), 'o-', ms=2, lw=0.8, color=col, label='Ra$=$'+lab)
ax[2].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g')); ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[2].set_xticks([0.2, 0.5, 1.0]); ax[2].set_xlabel(r'$\Gamma$'); ax[2].set_ylabel(r'Nu$/\max_\Gamma$Nu'); ax[2].legend(frameon=False, fontsize=6, loc='lower center'); ax[2].set_title(r'(c) Pr$=\infty$, Nu($\Gamma$)', fontsize=9)
# (d) local slope d ln Nu / d ln Gamma against Gamma Ra^(1/5)
from steady_shoulder import slopes
for (f, Rv, lab), col in zip(SC, cols):
    Gm, sl = slopes((ROOT + '/data/steady/') + f); ax[3].semilogx(Gm*Rv**0.2, sl, 'o-', ms=2, lw=0.8, color=col, label='Ra$=$'+lab)
ax[3].axhline(0, color='k', lw=0.5); ax[3].axvline(13.7, color='0.5', lw=0.8, ls=':'); ax[3].axvline(10.4, color='C3', lw=0.8, ls='--')
ax[3].text(14.0, 0.47, 'Pr$=1$\noptimum', fontsize=6, color='0.4'); ax[3].text(7.2, 0.42, 'wide-roll\noptimum,\nRa$\\leq10^5$', fontsize=6, color='C3')
ax[3].set_ylim(-0.4, 0.6); ax[3].set_xlabel(r'$\Gamma\,$Ra$^{1/5}$'); ax[3].set_ylabel(r'$d\ln$Nu$/d\ln\Gamma$')
ax[3].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g')); ax[3].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[3].set_xticks([3, 5, 10, 20])
ax[3].legend(frameon=False, fontsize=5.5, loc='upper right', bbox_to_anchor=(1.0, 0.78), ncol=2, columnspacing=0.8, handlelength=1.5); ax[3].set_title(r'(d) local slope along the roll family', fontsize=9)
fig.tight_layout(); fig.savefig((ROOT + '/figs/fig4_steady.pdf')); fig.savefig((ROOT + '/figs/fig4_steady.png'), dpi=200); print('ok')
