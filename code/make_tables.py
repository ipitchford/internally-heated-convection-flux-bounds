import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bound_curve import make_Phi, B_of_Lambda
AR = {'ns': lambda R: 0.1555*R**(-2/3) - 0.7104*R**(-0.5)*np.log(1-1.5334*R**(-1/3)),
      'sf': lambda R: 0.0131*R**(-40/29) - 0.0312*R**(-35/29)*np.log(1-0.6301*R**(-10/29))}
Phi = {(bc, n): make_Phi(bc, n)[0] for bc in ['ns', 'sf'] for n in [1, 2]}
facts = {}
# range of B*(L lnL)^(1/3) over Lambda >= e^3
L = np.geomspace(np.e**3, 1e12, 200)
for n in [1, 2]:
    v = B_of_Lambda(Phi[('ns', n)], L)*(L*np.log(L))**(1/3); facts[f'ns_n{n}_range'] = [float(v.min()), float(v.max())]
    v = B_of_Lambda(Phi[('sf', n)], L)*np.sqrt(L); facts[f'sf_n{n}_range'] = [float(v.min()), float(v.max())]
Tmax = json.load(open(ROOT + '/paper/facts_constants.json'))['Tmax_d2_period4']
facts['Tmax_d2_period4'] = Tmax
for R in [1e5, 1e6, 1e9, 1e12]:
    for bc in ['ns', 'sf']:
        facts[f'{bc}_R{R:.0e}_T1/8'] = float(B_of_Lambda(Phi[(bc, 1)], [R/16])[0])
        facts[f'{bc}_R{R:.0e}_Tmax'] = float(B_of_Lambda(Phi[(bc, 1)], [0.5*R*Tmax])[0])
        facts[f'{bc}_R{R:.0e}_AR'] = float(AR[bc](R))
print(json.dumps(facts, indent=1))
# Table: B for Lambda values, d=2,3, both bc
rows = []
for lg in [3, 4, 6, 8, 10, 12]:
    Lm = 10.0**lg; r = [f'$10^{{{lg}}}$']
    for bc in ['ns', 'sf']:
        for n in [1, 2]:
            r.append(f'{B_of_Lambda(Phi[(bc, n)], [Lm])[0]:.3e}'.replace('e-0', '\\times10^{-').replace('e-', '\\times10^{-') + '}')
    rows.append(' & '.join(['$'+x+'$' if i else x for i, x in enumerate(r)]) + r' \\')
tab = r"""\begin{table}[t]\centering
\caption{Uncertified numerical evaluation of $\mathfrak B[\Lambda\Phi_{\rm bc}]$ of Theorem~\ref{thm:main} ($\Lambda=\frac12R\mathcal T$).}\label{tab:bound}
\begin{tabular}{lcccc}\hline
$\Lambda$ & no-slip, $d=2$ & no-slip, $d=3$ & stress-free, $d=2$ & stress-free, $d=3$\\\hline
""" + '\n'.join(rows) + r"""
\hline\end{tabular}\end{table}
"""
open((ROOT + '/paper/tab_bound.tex'), 'w').write(tab)
# DNS table (continuation runs; physical-time averages)
ch = json.load(open((ROOT + '/data/dns/dns_checks.json')))
def tex_sci(x):
    if x == 0: return '$0$'
    e = int(np.floor(np.log10(abs(x)))); m = x/10**e
    return f'${"-" if m < 0 else ""}{abs(m):.1f}\\times10^{{{e}}}$'
def sgn(x, f='.4f'):
    return ('$-$' if x < 0 else '') + format(abs(x), f)
rows = []
for c in ch:
    imb = c['FB'] + c['FT'] - 1
    lo, hi = c['FB_block_range']
    rows.append(' & '.join([
        'no-slip' if c['bc'] == 'ns' else 'stress-free', f"$10^{{{int(round(np.log10(c['R'])))}}}$",
        f"$[{c['window'][0]:.0f},{c['window'][1]:.0f}]$", f"{c['FB']:.4f}", f"{lo:.4f}--{hi:.4f}", tex_sci(imb),
        f"{c['Tsup']:.4f}", f"{c['max_J_over_WTheta']:.3f}", f"{c['max_W_over_bound']:.3f}",
        sgn(c['min_J_z_lt_0p2']), f"{c['bound_window_T']:.4f}", sgn(c['finite_window_bound'])]) + r' \\')
tab = r"""\begin{table}[t]\centering\footnotesize
\caption{Simulations ($d=2$, period 4, $256\times64$ modes), averaged in physical time over the window
$[t_0,t_0+\tau]$ shown, and the descriptive range over eight sub-blocks (not a long-time confidence interval). The sixth column is the
window's heat imbalance, which equals $-\Delta\langle T\rangle/\tau$ exactly. $\mathcal T_\tau$ is the
maximum temperature over the window. The two ratios must not exceed one (Lemmas~\ref{lem:IF} and \ref{lem:rep}); $\min J$ is
over $z<0.2$. The last two columns are $\mathfrak B[\frac12R\mathcal T_\tau\Phi]$ and the finite-window value
$\mathfrak B[\frac12R\mathcal T_\tau\Phi]-\mathcal T_\tau/\tau$ whose exact formula occurs in Lemma~\ref{lem:IF}; these decimal evaluations are uncertified.}\label{tab:dns}
\setlength{\tabcolsep}{3pt}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lccccccccccc}\hline
walls & $R$ & window & $\overline{F_B}$ & block means & $\overline{F_B}+\overline{F_T}-1$ & $\mathcal T_\tau$ & $\max\frac{J}{W\Theta}$ & $\max\frac{W}{\frac12R\mathcal T_\tau\Phi}$ & $\min J$ & $\mathfrak B$ & $\mathfrak B-\frac{\mathcal T_\tau}{\tau}$\\\hline
""" + '\n'.join(rows) + r"""
\hline\end{tabular}}\end{table}
"""
open((ROOT + '/paper/tab_dns.tex'), 'w').write(tab)
json.dump(facts, open((ROOT + '/paper/facts_bounds.json'), 'w'), indent=1)

# ---------------- steady-state table and facts
ih = json.load(open((ROOT + '/data/steady/ih_ns.json'))); rb = json.load(open((ROOT + '/data/steady/rb_ns.json')))
ARd = np.array([(5e4,.0782),(1e5,.0571),(5e5,.162),(1e6,.133),(5e6,.183),(1e7,.213),(5e7,.242),(1e8,.268),(5e8,.275),(1e9,.313),(5e9,.329)])
lr, lf = np.log(ARd[:, 0]), np.log(0.5-ARd[:, 1])
sf = {}
cmp = []
for r in ih:
    if r['Ra'] < 1e5: continue
    turb = np.exp(np.interp(np.log(r['Ra']), lr, lf)); cmp.append((r['Ra'], r['FB'], turb, r['FB']/turb-1))
sf['ih_vs_turb'] = cmp
last = ih[-1]; sf['ih_last'] = dict(R=last['Ra'], FB=last['FB'], bound_T18=float(B_of_Lambda(Phi[('ns',1)], [last['Ra']/16])[0]))
sf['ih_last']['factor'] = sf['ih_last']['FB']/sf['ih_last']['bound_T18']
# unconditional comparison: Corollary 2.4 with T_max for the roll's own period (enlarged to at least 1)
import math
Gt = math.ceil(1/last['Gamma'])*last['Gamma']; Tm_roll = 1 + np.e/(2*np.pi**2)*max(2, np.log(8*Gt))
sf['ih_last'].update(Tmax_roll=Tm_roll, bound_Tmax=float(B_of_Lambda(Phi[('ns',1)], [0.5*last['Ra']*Tm_roll])[0]))
sf['ih_last']['factor_unconditional'] = sf['ih_last']['FB']/sf['ih_last']['bound_Tmax']
Rs = np.array([r['Ra'] for r in ih]); Ts = np.array([r['Tmax'] for r in ih])
Bm = B_of_Lambda(Phi[('ns',1)], 0.5*Rs*Ts); m = Rs >= 1e6
sf['bound_measuredT_fit_exponent'] = float(np.polyfit(np.log(Rs[m]), np.log(Bm[m]), 1)[0])
sf['Tmax_fit_exponent'] = float(np.polyfit(np.log(Rs[m]), np.log(Ts[m]), 1)[0])
# sensitivity of the power-law fits to the lower cutoff of the fitted range
rbg0 = json.load(open((ROOT + '/data/steady/rb_ns_best_found.json')))
sf['fit_sensitivity'] = {}
for name, data, key in (('ih', ih, 'FB'), ('rb', rbg0, 'Nu')):
    Rv = np.array([r['Ra'] for r in data]); V = np.array([r[key] for r in data])
    for cut in (1e5, 1e6, 1e7):
        mm = Rv >= cut*0.999; pf = np.polyfit(np.log(Rv[mm]), np.log(V[mm]), 1)
        sf['fit_sensitivity'][f'{name}_from_{cut:.0e}'] = dict(n=int(mm.sum()), exponent=float(pf[0]), prefactor=float(np.exp(pf[1])),
                                                              range=[float(Rv[mm].min()), float(Rv[mm].max())])
print(json.dumps(sf, indent=1))
json.dump(sf, open((ROOT + '/paper/facts_steady.json'), 'w'), indent=1)
# two-part table: internally heated (left) and Rayleigh-Benard (right), one row per Rayleigh number
rbg = json.load(open((ROOT + '/data/steady/rb_ns_best_found.json')))
def sci(x):
    e = int(np.floor(np.log10(x)+1e-9)); m = x/10**e
    return f'$10^{{{e}}}$' if abs(m-1) < 1e-6 else f'${m:.2f}\\times10^{{{e}}}$'
L = [f"{sci(r['Ra'])} & {r['Nz']}, {r['Nx']} & {r['Gamma']:.4f} & {r['FB']:.5f} & {r['Tmax']:.4f}" for r in ih]
Rr = []
for q in rbg:
    tag = '$^\\dagger$' if q['source'].startswith('refined') else ''
    Rr.append(f"{sci(q['Ra'])} & {q['Nz']}, {q['Nx']} & {q['Gamma']:.4f}{tag} & {q['Nu']:.5f}")
n = max(len(L), len(Rr)); L += [' & & & & ']*(n-len(L)); Rr += [' & & & ']*(n-len(Rr))
rows = [a + ' & ' + b + r' \\' for a, b in zip(L, Rr)]
import glob
rc = [c for f in sorted(glob.glob((ROOT + '/data/steady/res_check*.json'))) for c in json.load(open(f))]
rcl = [c for c in rc if c['mode'] == 'ih']; rcr = sorted([c for c in rc if c['mode'] == 'rb'], key=lambda c: (c['Ra'], c['Nz']))
L2 = [f"{sci(c['Ra'])} & {c['Nz']}, {c['Nx']} & & {c['FB']:.8f} & " for c in rcl]
R2 = [f"{sci(c['Ra'])} & {c['Nz']}, {c['Nx']} & & {c['Nu']:.6f}" for c in rcr]
n = max(len(L2), len(R2)); L2 += [' & & & & ']*(n-len(L2)); R2 += [' & & & ']*(n-len(R2))
resrows = '\\hline\n\\multicolumn{5}{l|}{Resolution check} & \\multicolumn{4}{l}{Resolution check}\\\\\n' + '\n'.join(a + ' & ' + b + r' \\' for a, b in zip(L2, R2))
tab = r"""\begin{table}[t]\centering\footnotesize
\caption{Extremal steady rolls at infinite Prandtl number, no-slip walls ($d=2$). Left: internally heated,
$F_B$-minimising period $\Gamma^*$. Right: Rayleigh--B\'enard, Nu-maximising period. Resolution $N_z,N_x$ is per
period. $^\dagger$Largest maximum found by continuation in $\Gamma$ (Section~\ref{sec:steady}); the narrow-roll
local maximum at this Rayleigh number is $\mathrm{Nu}=4.79376$ at $\Gamma=0.7440$.}\label{tab:steady}
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{lcccc|lccc}\hline
\multicolumn{5}{c|}{internally heated} & \multicolumn{4}{c}{Rayleigh--B\'enard}\\
$R$ & $N_z,N_x$ & $\Gamma^*$ & $F_B$ & $\mathcal T$ & Ra & $N_z,N_x$ & $\Gamma^*$ & Nu\\\hline
""" + '\n'.join(rows) + '\n' + resrows + r"""
\hline\end{tabular}\end{table}
"""
open((ROOT + '/paper/tab_steady.tex'), 'w').write(tab)
