"""Explicit constants of Theorems 2.1 and 2.3 and Proposition 2.2, from exact nu_j values and the computed
(not interval-certified) slab-correction suprema C_ns, C_sf. Writes paper/facts_constants.json."""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT
from scipy.special import gamma as G
from scipy.optimize import brentq
nu = {2: dict(nu1=2/np.pi, nu2=3*np.sqrt(3)/(2*np.pi), nu3=6/np.pi),          # n = 1 (d = 2)
      3: dict(nu1=4/(3*np.sqrt(3)), nu2=1.1154192037, nu3=2.7693233725)}     # n = 2 (d = 3)
cw = json.load(open(os.path.join(ROOT, 'code', 'cw_values.json')))
cwd = {2: cw['cw_d2_methodA'], 3: cw['cw_d3_methodB']}
r = (np.e-2)/(np.e+2)
out = {}
for d in (2, 3):
    n = d-1; v = nu[d]
    cc = json.load(open(os.path.join(ROOT, 'code', f'corr_constants2_n{n}.json')))
    Cns, Csf = cc['ns']['sup'], cc['sf']['sup']
    A = (5*v['nu2'] + 3*v['nu1'] + v['nu3'])/8; B = 1.5*v['nu1'] + Cns/2
    # kappa_0: largest kappa with kappa^3 [7A + 3B + 3A ln^+(1/(2 kappa))] <= 27
    F = lambda k: k**3*(7*A + 3*B + 3*A*max(0.0, np.log(1/(2*k)))) - 27
    k0 = (27/(7*A+3*B))**(1/3) if 7*A+3*B <= 216 else brentq(F, 1e-6, 0.5)
    Csfp = v['nu1']/2 + Csf
    # second route: the Carlson majorants of Lemma 3.6(iii) (carlson_sup.py), which bound every norm from above
    cs = json.load(open(os.path.join(ROOT, 'code', 'carlson_sup.json')))
    CnsC, CsfC = cs[f'ns_n{n}']['C_majorant'], cs[f'sf_n{n}']['C_majorant']
    BC = 1.5*v['nu1'] + CnsC/2; k0C = (27/(7*A+3*BC))**(1/3); CsfpC = v['nu1']/2 + CsfC
    carl = dict(C_ns=CnsC, C_sf=CsfC, B=BC, kappa0=k0C, c_ns=r*k0C, C_sf_prime=CsfpC, c_sf=r*np.sqrt(2/CsfpC),
                sf_threshold=8/CsfpC, seven_A_3B_le_216=bool(7*A+3*BC <= 216))
    out[f'd{d}'] = dict(nu1=v['nu1'], nu2=v['nu2'], nu3=v['nu3'], A=A, C_ns=Cns, B=B, kappa0=k0, c_ns=r*k0,
                        C_sf=Csf, C_sf_exact_lower=v['nu1']/2, C_sf_prime=Csfp, c_sf=r*np.sqrt(2/Csfp),
                        c_sf_if_Csf_eq_nu1_half=r*np.sqrt(2/v['nu1']), sf_threshold=8/Csfp, c_w=cwd[d],
                        prop22_ns=G(2/3)/G(1/3)*(18/v['nu2'])**(1/3), prop22_sf=np.sqrt(2/(np.pi*cwd[d])),
                        N_d=2**(d+1), wall_coeff=v['nu2']/2, carlson_route=carl,
                        checks=dict(two_B_ge_3nu1_ge_A=bool(2*B >= 3*v['nu1'] >= A),
                                    A3_plus_B_gt_24e3=bool(A/3 + B > 24*np.exp(-3)),
                                    kappa0_cubed_lt_4=bool(k0**3 < 4), seven_A_3B_le_216=bool(7*A+3*B <= 216)))
# explicit temperature ceiling for the simulated geometry (d = 2, period 4)
Tmax = lambda d, vol: 1 + np.e/(2*np.pi**2)*max(2, np.log(2**(d+1)*vol))
out['Tmax_d2_period4'] = Tmax(2, 4.0)
out['Tmax_d2_period1'] = Tmax(2, 1.0)
json.dump(out, open(os.path.join(ROOT, 'paper', 'facts_constants.json'), 'w'), indent=1)
for k, v in out.items():
    print(k, {kk: (round(float(vv), 6) if isinstance(vv, (float, np.floating)) else vv) for kk, vv in v.items()} if isinstance(v, dict) else round(v, 6))
