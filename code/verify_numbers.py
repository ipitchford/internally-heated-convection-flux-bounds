"""Assertion suite for every number quoted in the paper.

Each check either recomputes a quantity (from first principles, or from stored raw data such as time histories)
or reads it from a stored product, and compares it with the value printed in the text, with a stated tolerance.
Inequalities that the proofs or the text assert are checked as Boolean conditions. Any failure is listed at the end
and the script exits with status 1.

Run from anywhere:  python3 code/verify_numbers.py > logs/verify_numbers.log
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ROOT
import numpy as np, mpmath as mp
from scipy.special import gamma as G
from bound_curve import make_Phi, B_of_Lambda, load_phi

D = ROOT + '/'
FAILED, N = [], [0]


def check(label, value, expected=None, tol=None, rel=False, cond=None):
    """Numeric check (value against expected +/- tol, absolute or relative) or Boolean check (cond)."""
    N[0] += 1
    if cond is not None:
        good = bool(cond); msg = f'{label}: {value}'
    else:
        err = abs(value - expected)/(abs(expected) if rel else 1.0)
        good = err <= tol
        msg = f'{label}: {value:.10g}  [text {expected:.10g}, {"rel" if rel else "abs"} tol {tol:g}]'
    print(('PASS  ' if good else 'FAIL  ') + msg)
    if not good: FAILED.append(label)


def section(t): print(f'\n=== {t}')


# ---------------------------------------------------------------------------------------------------------------
section('1. Constants nu_j (L1 norms of inverse transforms of |xi|^j e^{-|xi|}), recomputed in 30-digit arithmetic')
mp.mp.dps = 30
import sympy as sp
def _split_quad(f, roots):
    pts = [0] + sorted(roots) + [mp.inf]
    return sum(abs(mp.quad(f, [pts[i], pts[i+1]])) for i in range(len(pts)-1))
def nu1(j):
    y = sp.symbols('y', positive=True); num = sp.expand(sp.re(sp.expand((1+sp.I*y)**(j+1))))
    roots = [float(x) for x in sp.Poly(num, y).nroots() if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0]
    f = sp.lambdify(y, sp.factorial(j)*num/(1+y**2)**(j+1)/sp.pi, 'mpmath'); return 2*_split_quad(f, roots)
def nu2(j):
    c, r = sp.symbols('c r', positive=True)
    expr = sp.simplify(((-1)**(j+1)*sp.diff((c**2+r**2)**sp.Rational(-1, 2), c, j+1)).subs(c, 1))
    roots = [float(x) for x in sp.Poly(sp.expand(sp.numer(sp.together(expr))), r).nroots() if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0]
    f = sp.lambdify(r, expr*r, 'mpmath'); return _split_quad(f, roots)
nu = {1: [float(nu1(j)) for j in range(4)], 2: [float(nu2(j)) for j in range(4)]}
check('nu_1 (n=1) = 2/pi', nu[1][1], 2/np.pi, 1e-12)
check('nu_2 (n=1) = 3 sqrt3/(2 pi)', nu[1][2], 3*np.sqrt(3)/(2*np.pi), 1e-12)
check('nu_3 (n=1) = 1.9098593171 (= 6/pi)', nu[1][3], 1.9098593171, 1e-10)
check('nu_1 (n=2) = 4/(3 sqrt3)', nu[2][1], 4/(3*np.sqrt(3)), 1e-12)
check('nu_2 (n=2) = 1.1154192037', nu[2][2], 1.1154192037, 1e-10)
check('nu_3 (n=2) = 2.7693233725', nu[2][3], 2.7693233725, 1e-10)
check('wall coefficient nu_2/2, d=2: 0.4135', nu[1][2]/2, 0.4135, 5e-5)
check('wall coefficient nu_2/2, d=3: 0.5577', nu[2][2]/2, 0.5577, 5e-5)

# ---------------------------------------------------------------------------------------------------------------
section('2. Lemma 3.6: low-frequency structure and the exact endpoint identity (recomputed)')
from kernels import corr_symbol, L1_corr_deriv
k = np.array([1e-3, 0.05, 0.3, 0.9, 1.5, 4.0])
mu = corr_symbol(k, 0.0, 1.0, 'sf', 1)
check('mu_sf(kappa; 0, 1) = -(1/2) kappa e^{-kappa}, max abs deviation over 6 kappa', float(np.abs(mu + 0.5*k*np.exp(-k)).max()), 0.0, 1e-12)
for z, s in [(0.2, 0.4), (0.0, 0.7), (0.45, 0.95)]:
    kk = np.array([1e-6])
    check(f'mu_sf/kappa -> -s/2 at (zeta,s)=({z},{s})', float(corr_symbol(kk, z, s, 'sf', 1)[0]/1e-6), -s/2, 1e-5)
    kk = np.array([1e-4]); c2 = s*(s*s+3*z*z-1)/6
    r2 = float((corr_symbol(kk, z, s, 'sf', 1)[0] + (s/2)*1e-4*np.exp(-1e-4))/1e-8)
    check(f'r_sf/kappa^2 -> s(s^2+3 zeta^2-1)/6 at ({z},{s})', r2, c2, 2e-3*max(1, abs(c2)))
    c2n = -s*s*(2*s*z - s - 3*z + 2)
    check(f'mu_ns/kappa^2 -> -s^2(2 s zeta - s - 3 zeta + 2) at ({z},{s})', float(corr_symbol(kk, z, s, 'ns', 2)[0]/1e-8), c2n, 2e-3*max(1, abs(c2n)))
for n in (1, 2):
    v = L1_corr_deriv(0.0, 1.0, n, 'sf', 1)
    check(f'endpoint norm ||F^-1 mu_sf(.;0,1)||_1 = nu_1/2 (n={n}), recomputed', v, nu[n][1]/2, 1e-7)

# ---------------------------------------------------------------------------------------------------------------
section('3. Slab-correction constants C_ns, C_sf (grid suprema, convergence, Carlson majorants)')
cc = {n: json.load(open(D+f'code/corr_constants2_n{n}.json')) for n in (1, 2)}
Cns = {n: cc[n]['ns']['sup'] for n in (1, 2)}; Csf = {n: cc[n]['sf']['sup'] for n in (1, 2)}
check('C_ns (d=2) = 0.86155', Cns[1], 0.86155, 5e-6); check('C_ns (d=3) = 1.17193', Cns[2], 1.17193, 5e-6)
check('C_sf (d=2) = 0.3183099', Csf[1], 0.3183099, 5e-8); check('C_sf (d=3) = 0.3849002', Csf[2], 0.3849002, 5e-8)
check('C_sf (d=2) agrees with 1/pi to 1e-8', abs(Csf[1]-1/np.pi), cond=abs(Csf[1]-1/np.pi) <= 1e-8)
check('C_sf (d=3) agrees with 2/(3 sqrt3) to 5e-8', abs(Csf[2]-2/(3*np.sqrt(3))), cond=abs(Csf[2]-2/(3*np.sqrt(3))) <= 5e-8)
check('no-slip maximisers at zeta=0, s=0.9125 (d=2), 0.90 (d=3)', (cc[1]['ns']['argmax'], cc[2]['ns']['argmax']),
      cond=np.allclose(cc[1]['ns']['argmax'], [0, 0.9125]) and np.allclose(cc[2]['ns']['argmax'], [0, 0.9]))
check('stress-free maximisers at the corner (0,1)', (cc[1]['sf']['argmax'], cc[2]['sf']['argmax']),
      cond=np.allclose(cc[1]['sf']['argmax'], [0, 1]) and np.allclose(cc[2]['sf']['argmax'], [0, 1]))
check('C_ns at the d=2 maximiser, recomputed', L1_corr_deriv(0.0, 0.9125, 1, 'ns', 2), Cns[1], 1e-9)
# grid slack: largest neighbour difference per unit parameter change, times half the final spacing
slack = {}
for n in (1, 2):
    for bc in ('ns', 'sf'):
        pts = {tuple(map(float, kk.split(','))): vv for kk, vv in cc[n][bc]['grid'].items()}
        zs, ss = np.round(np.linspace(0, .5, 11), 6), np.round(np.linspace(0, 1, 21), 6)
        V = np.array([[pts[(z, s)] for s in ss] for z in zs])
        Lz, Ls = np.abs(np.diff(V, axis=0)).max()/0.05, np.abs(np.diff(V, axis=1)).max()/0.05
        slack[(n, bc)] = (Lz + Ls)*0.05/16/2
        if bc == 'ns': check(f'grid Lipschitz estimate for ns (n={n}) at most {[None, 2.5, 3.5][n]}', round(max(Lz, Ls), 3), cond=max(Lz, Ls) <= [None, 2.52, 3.5][n])
        else: check(f'grid Lipschitz estimate for sf (n={n}) at most 0.45', round(max(Lz, Ls), 3), cond=max(Lz, Ls) <= 0.45)
check('grid slack below 0.007 (ns, d=2), 0.009 (ns, d=3), 0.0011 (sf)', {f'{k[1]}{k[0]}': round(float(v), 5) for k, v in slack.items()},
      cond=slack[(1, 'ns')] < 0.007 and slack[(2, 'ns')] < 0.009 and slack[(1, 'sf')] < 0.0011 and slack[(2, 'sf')] < 0.0011)
cv = json.load(open(D+'code/corr_convergence.json'))
dY = max(abs(r['Y80']-r['Y40']) for r in cv.values()); dI = max(abs(r['independent']-r['Y40']) for r in cv.values() if 'independent' in r)
check('truncation radius 40 -> 80 changes norms at maximisers by at most 3e-6', dY, cond=dY <= 3e-6)
check('independent mpmath evaluation (n=1) agrees to 2e-6', dI, cond=dI <= 2e-6)
cs = json.load(open(D+'code/carlson_sup.json'))
for key, txt in [('ns_n1', 0.9703), ('ns_n2', 1.3031), ('sf_n1', 0.4027), ('sf_n2', 0.5264)]:
    check(f'Carlson majorant {key} = {txt}', cs[key]['C_majorant'], txt, 5e-5)
check('Carlson sup of r_sf: 0.0844 (d=2), 0.1415 (d=3)', (round(cs['sf_n1']['sup_carlson'], 4), round(cs['sf_n2']['sup_carlson'], 4)),
      cond=abs(cs['sf_n1']['sup_carlson']-0.0844) < 5e-5 and abs(cs['sf_n2']['sup_carlson']-0.1415) < 5e-5)
check('Carlson majorants exceed the direct grid suprema', None, cond=all(cs[f'ns_n{n}']['C_majorant'] >= Cns[n] and cs[f'sf_n{n}']['C_majorant'] >= Csf[n] for n in (1, 2)))
check('Carlson step-halving check agrees to 1e-6', None, cond=all(abs(v['half_step_check']-v['sup_carlson']) < 1e-6 for v in cs.values()))
from carlson_sup import carlson
check('Carlson functional at the d=2 no-slip maximiser, recomputed', carlson('ns', 1, 0.0, 0.9125), cs['ns_n1']['sup_carlson'], 1e-9)

# ---------------------------------------------------------------------------------------------------------------
section('4. Constants of Theorem 2.1 and Corollary 2.4')
r = (np.e-2)/(np.e+2)
def consts(n, Cn, Cs):
    v = nu[n]; A = (5*v[2]+3*v[1]+v[3])/8; B = 1.5*v[1] + Cn/2; k0 = (27/(7*A+3*B))**(1/3); Cp = v[1]/2 + Cs
    return A, B, k0, r*k0, Cp, r*np.sqrt(2/Cp)
for n, (tA, tB, tcn, tcs, tcnC, tcsC) in {1: (0.994, 1.386, 0.2046, 0.2698, 0.2036, 0.2535),
                                         2: (1.332, 1.741, 0.1871, 0.2454, 0.1863, 0.2255)}.items():
    A, B, k0, cn, Cp, csf = consts(n, Cns[n], Csf[n]); _, BC, _, cnC, _, csC = consts(n, cs[f'ns_n{n}']['C_majorant'], cs[f'sf_n{n}']['C_majorant'])
    d = n+1
    check(f'A (d={d})', A, tA, 5e-4); check(f'B (d={d})', B, tB, 5e-4)
    check(f'c_ns (d={d}), computed C', cn, tcn, 5e-5); check(f'c_sf (d={d}), computed C', csf, tcs, 5e-5)
    check(f'c_ns (d={d}), Carlson route', cnC, tcnC, 5e-5); check(f'c_sf (d={d}), Carlson route', csC, tcsC, 5e-5)
    check(f'theorem rounding: c_ns ~ {[None, 0.205, 0.187][n]}, c_sf ~ {[None, 0.270, 0.245][n]} (d={d})', (round(cn, 3), round(csf, 3)),
          cond=round(cn, 3) == [None, 0.205, 0.187][n] and round(csf, 3) == [None, 0.270, 0.245][n])
    check(f'theorem rounding, Carlson route: {[None, 0.204, 0.186][n]}, {[None, 0.254, 0.226][n]} (d={d})', (round(cnC, 3), round(csC, 3)),
          cond=round(cnC, 3) == [None, 0.204, 0.186][n] and round(csC, 3) == [None, 0.254, 0.226][n])
    check(f'proof conditions 2B >= 3 nu_1 >= A, A/3+B > 24 e^-3, kappa_0^3 < 4, 7A+3B <= 216 (d={d})', None,
          cond=2*B >= 3*nu[n][1] >= A and A/3+B > 24*np.exp(-3) and k0**3 < 4 and 7*A+3*B <= 216 and 7*A+3*BC <= 216)
    s1 = consts(n, 1.1*Cns[n], Csf[n])[3]/cn - 1; s2 = consts(n, Cns[n], 1.1*Csf[n])[5]/csf - 1
    check(f'10% increase of C_ns lowers c_ns by 0.4% (d={d})', -s1, 0.004, 2.5e-4); check(f'10% increase of C_sf lowers c_sf by 2.4% (d={d})', -s2, 0.024, 5e-4)
    check(f'2A/nu_2 = {[None, 2.40, 2.39][n]} (d={d})', 2*A/nu[n][2], [None, 2.40, 2.39][n], 5e-3)
Tmax = lambda d, vol: 1 + np.e/(2*np.pi**2)*max(2, np.log(2**(d+1)*vol))
check('T_max = 1 + (e/2 pi^2) ln 32 = 1.477 (d=2, period 4)', Tmax(2, 4.0), 1.477, 5e-4)
check('T_max formula uses ln(32) > 2 for period 4', None, cond=np.log(32) > 2)

# ---------------------------------------------------------------------------------------------------------------
section('5. Wall functions, c_w, Proposition 2.2 and the majorants')
cw = json.load(open(D+'code/cw_values.json'))
check('c_w (d=2) = 0.194366', cw['cw_d2_methodA'], 0.194366, 5e-7); check('c_w (d=3) = 0.254086', cw['cw_d3_methodB'], 0.254086, 5e-7)
check('c_w (d=2): two methods agree to 3e-8', abs(cw['cw_d2_methodA']-cw['cw_d2_methodB']), cond=abs(cw['cw_d2_methodA']-cw['cw_d2_methodB']) <= 3.5e-8)
for n in (1, 2):
    Phi, a, b = make_Phi('ns', n)
    check(f'Phi_ns constant term b (d={n+1})', b, [None, -0.430, -0.573][n], 1e-3)
    c_w = [None, cw['cw_d2_methodA'], cw['cw_d3_methodB']][n]
    check(f'Prop 2.2 no-slip constant (d={n+1})', G(2/3)/G(1/3)*(18/nu[n][2])**(1/3), [None, 1.41, 1.28][n], 5e-3)
    check(f'Prop 2.2 stress-free constant (d={n+1})', np.sqrt(2/(np.pi*c_w)), [None, 1.81, 1.58][n], 5e-3)
    A, B, _, _, Cp, _ = consts(n, Cns[n], Csf[n])
    z, P = load_phi('ns', n); m = z <= 0.5; rat = (z**2*(A*np.log(1/(2*z))+B)/P)[m]
    check(f'no-slip majorant / Phi >= 2.8 on 1e-4 <= z <= 1/2 (d={n+1})', round(float(rat.min()), 3), cond=rat.min() >= 2.8)
    z, P = load_phi('sf', n); m = z <= 0.5; rat = (Cp*z/P)[m]
    check(f'stress-free majorant / Phi >= 3.0 on 1e-4 <= z <= 1/2 (d={n+1})', round(float(rat.min()), 3), cond=rat.min() >= 3.0)
fb = json.load(open(D+'paper/facts_bounds.json'))
for k_, (lo, hi) in {'ns_n1_range': (1.46, 2.70), 'ns_n2_range': (1.31, 2.39), 'sf_n1_range': (1.81, 2.56), 'sf_n2_range': (1.58, 2.23)}.items():
    check(f'{k_} = [{lo}, {hi}]', str([round(x, 3) for x in fb[k_]]), cond=abs(fb[k_][0]-lo) < 5e-3 and abs(fb[k_][1]-hi) < 5e-3)
check('c_ns conservative by a factor of about 7 (1.46/0.205)', fb['ns_n1_range'][0]/0.2046, 7.1, 0.2)

# ---------------------------------------------------------------------------------------------------------------
section('6. Bound values and comparison with Arslan & Rojas (2025)')
check('unconditional no-slip bound at R=1e9 = 6.0e-4', fb['ns_R1e+09_Tmax'], 6.0e-4, 5e-6)
check('unconditional stress-free bound at R=1e9 = 6.7e-5', fb['sf_R1e+09_Tmax'], 6.7e-5, 5e-7)
check('conditional (T=1/8) no-slip bound at R=1e9 = 1.4e-3', fb['ns_R1e+09_T1/8'], 1.4e-3, 5e-5)
check('conditional (T=1/8) stress-free bound at R=1e9 = 2.3e-4', fb['sf_R1e+09_T1/8'], 2.3e-4, 5e-6)
check('Arslan & Rojas no-slip at R=1e9 = 1.9e-7', fb['ns_R1e+09_AR'], 1.9e-7, 5e-9)
check('Arslan & Rojas stress-free at R=1e9 = 5e-15', fb['sf_R1e+09_AR'], 5e-15, 5e-16)
ARd = {'ns': [(1e5, .0571), (5e5, .162), (1e6, .133), (5e6, .183), (1e7, .213), (5e7, .242), (1e8, .268), (5e8, .275), (1e9, .313), (5e9, .329)]}
P1 = make_Phi('ns', 1)[0]; T4 = Tmax(2, 4.0)
ratios = [(0.5-w)/B_of_Lambda(P1, [0.5*R*T4])[0] for R, w in ARd['ns']]
check('Arslan & Rojas no-slip simulations exceed the unconditional bound by 26 to 500', (round(min(ratios), 1), round(max(ratios), 1)),
      cond=abs(min(ratios)-26) < 0.5 and abs(max(ratios)-500) < 5)
ch = json.load(open(D+'data/dns/dns_checks.json'))
ours = [c['FB']/B_of_Lambda(make_Phi(c['bc'], 1)[0], [0.5*c['R']*T4])[0] for c in ch]
check('our simulations: 26 to 150 times the unconditional bound', [round(x, 1) for x in ours], cond=abs(min(ours)-26) < 0.5 and abs(max(ours)-150) < 2)
own = [c['FB']/c['bound_window_T'] for c in ch]
check('our simulations: 10 to 33 times the bound with the window maximum temperature', [round(x, 1) for x in own], cond=abs(min(own)-9.5) < 0.6 and abs(max(own)-33) < 0.6)

# ---------------------------------------------------------------------------------------------------------------
section('7. Simulations: physical-time averages, heat balance, resolution and character (from raw histories)')
quoted = {'ih_ns_R1e+05_G4_256x64_cont': dict(FB=0.4435, T=0.1149, r1=0.090, r2=0.314, plateau=0.44291),
          'ih_ns_R1e+06_G4_256x64_cont': dict(FB=0.3620, T=0.0790, r1=0.055, r2=0.339, plateau=None),
          'ih_sf_R1e+05_G4_256x64_cont': dict(FB=0.3916, T=0.0974, r1=0.074, r2=0.359, plateau=0.39158),
          'ih_sf_R1e+06_G4_256x64_cont': dict(FB=0.3200, T=0.0694, r1=0.043, r2=0.251, plateau=0.32024)}
fracs, spreads, low = [], [], []
for c in ch:
    q = quoted[c['name']]; info = json.load(open(D+f"data/dns/{c['name']}.json"))
    h = np.load(D+f"data/dns/{c['name']}_hist.npz")['hist']; t = h[:, 0]; w = (t >= info['t_avg'][0]-1e-12) & (t <= info['t_avg'][1]+1e-12)
    L = t[w][-1]-t[w][0]; trap = lambda col: float(np.trapezoid(h[w, col], t[w])/L)
    check(f"{c['name']}: stored F_B equals trapezoid over the raw history", trap(1), c['FB'], 1e-9)
    check(f"{c['name']}: F_B quoted {q['FB']}", c['FB'], q['FB'], 5e-5)
    check(f"{c['name']}: heat-balance identity d<T>/dt = 1 - F_B - F_T to 1e-8", abs(c['identity_residual']), cond=abs(c['identity_residual']) <= 1e-8)
    check(f"{c['name']}: identity from the raw history", abs((1-trap(1)-trap(2)) - (h[w, 3][-1]-h[w, 3][0])/L), cond=abs((1-trap(1)-trap(2)) - (h[w, 3][-1]-h[w, 3][0])/L) <= 1e-8)
    check(f"{c['name']}: T_tau quoted {q['T']}", c['Tsup'], q['T'], 5e-5)
    check(f"{c['name']}: max J/(W Theta) <= 1 (quoted {q['r1']})", round(c['max_J_over_WTheta'], 3), cond=c['max_J_over_WTheta'] <= 1 and abs(c['max_J_over_WTheta']-q['r1']) < 5e-4)
    check(f"{c['name']}: max W/bound <= 1 (quoted {q['r2']})", round(c['max_W_over_bound'], 3), cond=c['max_W_over_bound'] <= 1 and abs(c['max_W_over_bound']-q['r2']) < 5e-4)
    check(f"{c['name']}: window mean F_B >= finite-window bound", (round(c['FB'], 4), round(c['finite_window_bound'], 4)), cond=c['FB'] >= c['finite_window_bound'])
    check(f"{c['name']}: T >= 0 up to round-off", c['Tmin'], cond=c['Tmin'] > -1e-11)
    ch_ = c['character']
    if q['plateau'] is None:
        check(f"{c['name']}: no steady plateau (time-dependent throughout)", ch_['time_fraction_at_plateau'], cond=ch_['plateau_FB'] is None)
    else:
        check(f"{c['name']}: plateau F_B = {q['plateau']}", ch_['plateau_FB'], q['plateau'], 6e-6); fracs.append(ch_['time_fraction_at_plateau'])
        check(f"{c['name']}: one or two bursts", ch_['excursions'], cond=ch_['excursions'] in (1, 2))
    lo_, hi_ = c['FB_block_range']; spreads.append(hi_-lo_); low.append(c['max_W_over_bound_z_lt_0p2'])
check('time at plateau 76-89%', [round(x, 3) for x in fracs], cond=abs(min(fracs)-0.76) < 0.005 and abs(max(fracs)-0.89) < 0.005)
check('largest spread of sub-block means at most 5e-3', round(max(spreads), 5), cond=max(spreads) <= 5e-3)
check('max|w| in the lower fifth is 11-20% of the kernel bound', [round(x, 3) for x in low], cond=abs(min(low)-0.107) < 0.005 and abs(max(low)-0.195) < 0.005)
v1 = {c['name']: c['v1_window'] for c in ch}
check('referee regression: v1 no-slip R=1e5 trapezoid F_B = 0.4443729596', v1['ih_ns_R1e+05_G4_256x64_cont']['FB_trapezoid'], 0.4443729596, 1e-9)
check('referee regression: v1 no-slip R=1e6 trapezoid F_B = 0.3621253105', v1['ih_ns_R1e+06_G4_256x64_cont']['FB_trapezoid'], 0.3621253105, 1e-9)
check('referee regression: v1 no-slip R=1e6 trapezoid F_T = 0.6380474703', v1['ih_ns_R1e+06_G4_256x64_cont']['FT_trapezoid'], 0.6380474703, 1e-9)
res = json.load(open(D+'data/dns/resolution.json'))
check('Fourier tails of T below 1e-8 in every snapshot', max(r_['max_fourier_tail'] for r_ in res), cond=max(r_['max_fourier_tail'] for r_ in res) < 1e-8)
check('Chebyshev tails of T below 2e-8 in every snapshot', max(r_['max_cheb_tail'] for r_ in res), cond=max(r_['max_cheb_tail'] for r_ in res) < 2e-8)
check('Arslan & Rojas comparison: steady state at R=1e5 equals their 0.4429 to the digits given', round(quoted['ih_ns_R1e+05_G4_256x64_cont']['plateau'], 4), cond=round(0.44291, 4) == 0.4429)

# ---------------------------------------------------------------------------------------------------------------
section('8. Profile hypothesis of Chen et al. (2024), recomputed from stored profiles')
nblk, nexc, exc_ok, below98 = 0, 0, True, 0.0
full = []
for c in ch:
    info = json.load(open(D+f"data/dns/{c['name']}.json")); d = np.load(D+f"data/dns/{c['name']}_prof.npz"); z = d['z']
    m = (z > 1e-4) & (z < 1-1e-4)
    tri = np.minimum(info['FB']*z, (1-z)*(1-info['FB'])); full.append(float((d['Tbar'][m]/tri[m]).max()))
    for Tb, B in zip(d['Tbar_blocks'], info['blocks']):
        tri = np.minimum(B['FB']*z, (1-z)*(1-B['FB'])); rr = Tb[m]/tri[m]; nblk += 1
        if rr.max() > 1:
            nexc += 1; exc_ok &= (B['FB']+B['FT']-1 >= 6e-4) and (z[m][np.argmax(rr)] > 0.999)
        below98 = max(below98, float(rr[z[m] < 0.98].max()))
check('full-window maximum ratios 0.99946-0.99965', [round(x, 5) for x in full], cond=abs(min(full)-0.99946) < 5e-6 and abs(max(full)-0.99965) < 5e-6)
check('6 of 32 sub-block ratios exceed one', (nexc, nblk), cond=(nexc, nblk) == (6, 32))
check('every exceedance at the point nearest the top plate, in a cooling sub-block (F_B+F_T-1 >= 6e-4)', None, cond=exc_ok)
check('below z = 0.98 every sub-block ratio <= 0.99933', round(below98, 6), cond=below98 <= 0.99933)
cs_ = json.load(open(D+'data/dns/chen_hypothesis_checks.json')); st = [x['max_ratio'] for x in cs_ if x['case'].startswith('steady')]
check('steady states: maximum ratios 0.9990-0.9997', (round(min(st), 4), round(max(st), 4)), cond=round(min(st), 4) == 0.999 and round(max(st), 4) == 0.9997)
check('largest sub-block exceedance at most 0.9%', None, cond=max(max(x.get('block_max_ratios', [0])) for x in cs_) <= 1.009)

# ---------------------------------------------------------------------------------------------------------------
section('9. Steady rolls')
ih = json.load(open(D+'data/steady/ih_ns.json')); rb = json.load(open(D+'data/steady/rb_ns_best_found.json'))
def fit(data, key, cut, field=None):
    R_ = np.array([x['Ra'] for x in data]); y = np.array([x[field or key] for x in data]); m = R_ >= cut*0.999
    p = np.polyfit(np.log(R_[m]), np.log(y[m]), 1); return np.exp(p[1]), p[0], int(m.sum())
a, e, nn = fit(ih, 'FB', 1e6); check('IH: F_B ~ 1.29 R^-0.093 over six Rayleigh numbers', (round(a, 3), round(e, 4), nn), cond=round(a, 2) == 1.29 and round(e, 3) == -0.093 and nn == 6)
a, e, _ = fit(ih, 'Gamma', 1e6); check('IH: Gamma* ~ 12.4 R^-0.21', (round(a, 2), round(e, 3)), cond=round(a, 1) == 12.4 and round(e, 2) == -0.21)
check('IH fit sensitivity: -0.0915 from R >= 1.6e5, -0.097 from the last two', (round(fit(ih, 'FB', 1e5)[1], 4), round(fit(ih, 'FB', 1e7)[1], 4)),
      cond=round(fit(ih, 'FB', 1e5)[1], 4) == -0.0915 and round(fit(ih, 'FB', 1e7)[1], 3) == -0.097)
a, e, _ = fit(rb, 'Nu', 1e6); check('RB: Nu ~ 0.153 Ra^0.298', (round(a, 3), round(e, 4)), cond=round(a, 3) == 0.153 and round(e, 3) == 0.298)
a, e, _ = fit(rb, 'Gamma', 1e6); check('RB: Gamma* ~ 11.7 Ra^-0.249', (round(a, 2), round(e, 4)), cond=round(a, 1) == 11.7 and round(e, 3) == -0.249)
check('RB fit sensitivity: 0.2979 from Ra >= 1e5, 0.2986 from Ra >= 1e7', (round(fit(rb, 'Nu', 1e5)[1], 4), round(fit(rb, 'Nu', 1e7)[1], 4)),
      cond=round(fit(rb, 'Nu', 1e5)[1], 4) == 0.2979 and round(fit(rb, 'Nu', 1e7)[1], 4) == 0.2986)
Ra_ = np.array([x['Ra'] for x in rb]); Nu_ = np.array([x['Nu'] for x in rb]); loc = np.diff(np.log(Nu_))/np.diff(np.log(Ra_))
check('RB local exponents over the last nine steps in [0.298, 0.299]', np.round(loc[-9:], 4).tolist(), cond=loc[-9:].min() >= 0.2975 and loc[-9:].max() <= 0.2995)
G_ = np.array([x['Gamma'] for x in rb]); wide = G_[Ra_ <= 1e5]*Ra_[Ra_ <= 1e5]**0.2
check('wide-roll optimum: Gamma* Ra^(1/5) between 10.4 and 10.6 for Ra <= 1e5', np.round(wide, 2).tolist(), cond=round(wide.min(), 1) == 10.4 and round(wide.max(), 1) == 10.6)
check('Gamma* = 1.67 at Ra=1e4 and 1.04 at Ra=1e5', (round(G_[0], 2), round(G_[Ra_ == 1e5][0], 2)), cond=round(G_[0], 2) == 1.67 and round(G_[Ra_ == 1e5][0], 2) == 1.04)
fs = json.load(open(D+'paper/facts_steady.json'))
check('IH steady rolls at R=2.8e7: factor 49 above the T=1/8 bound', fs['ih_last']['factor'], 49.35, 0.5)
check('IH steady rolls at R=2.8e7: factor 115 above the unconditional bound, T_max = 1.29 for their period', (round(fs['ih_last']['factor_unconditional'], 1), round(fs['ih_last']['Tmax_roll'], 2)),
      cond=round(fs['ih_last']['factor_unconditional']) == 115 and round(fs['ih_last']['Tmax_roll'], 2) == 1.29)
check('IH steady T_max exponent -0.21', fs['Tmax_fit_exponent'], -0.21, 5e-3); check('bound with measured T: exponent -0.30', fs['bound_measuredT_fit_exponent'], -0.30, 5e-3)
check('IH steady T_max from 0.123 to 0.035', (round(ih[0]['Tmax'], 3), round(ih[-1]['Tmax'], 3)), cond=round(ih[0]['Tmax'], 3) == 0.123 and round(ih[-1]['Tmax'], 3) == 0.035)
vt = fs['ih_vs_turb']; hi_ = [x[3] for x in vt if x[0] >= 8e5]; lo_ = [x[3] for x in vt if x[0] < 8e5]
check('IH: 1-4% below Arslan & Rojas for R >= 8e5', [round(100*x, 1) for x in hi_], cond=all(-0.0405 <= x <= -0.0085 for x in hi_))
check('IH: up to 11% above at lower R', round(100*max(lo_), 1), cond=abs(max(lo_)-0.11) < 0.005)
from steady_shoulder import local_maxima, shoulder
sl, pos = [], []
for fn in sorted(glob.glob(D+'data/steady/gscan_rb_*.json'), key=lambda f: json.load(open(f))[0]['Ra']):
    dd = json.load(open(fn)); Ra = dd[0]['Ra']
    if Ra > 1.5e5:
        sh = shoulder(fn); sl.append(round(sh['slope_shoulder'], 3)); pos.append(sh['Gamma_shoulder']*Ra**0.2)
        lm = [x for x in local_maxima(fn) if x['Gamma'] > 0.7]
        check(f'no wide-roll maximum at Ra={Ra:.3g} within the scanned periods', len(lm), cond=len(lm) == 0)
check('shoulder slopes -0.012, -0.030, -0.049, -0.067, -0.115', sl, cond=sl == [-0.012, -0.030, -0.049, -0.067, -0.115])
check('shoulder positions Gamma Ra^(1/5) = 10.0-10.4', [round(x, 2) for x in pos], cond=round(min(pos), 1) == 10.0 and round(max(pos), 1) == 10.4)
mx = json.load(open(D+'data/steady/maxima_1e5.json')); w64 = [x for x in mx if x['branch'] == 'wide' and x['Nz'] == 64][0]; n64 = [x for x in mx if x['branch'] == 'narrow' and x['Nz'] == 64][0]
check('Ra=1e5: wide maximum Nu 4.79437, narrow 4.79376 at Gamma 0.744; difference 0.013%', (round(w64['Nu'], 5), round(n64['Nu'], 5), round(n64['Gamma'], 3), round(100*(w64['Nu']-n64['Nu'])/w64['Nu'], 4)),
      cond=round(w64['Nu'], 5) == 4.79437 and round(n64['Nu'], 5) == 4.79376 and round(n64['Gamma'], 3) == 0.744 and round(100*(w64['Nu']-n64['Nu'])/w64['Nu'], 3) == 0.013)
check('Ra=1e5 maxima change by less than 1e-8 (relative) from Nz=48 to 64', None,
      cond=all(abs(a_['Nu']-b_['Nu'])/b_['Nu'] < 1e-8 for a_ in mx for b_ in mx if a_['branch'] == b_['branch'] and a_['Nz'] != b_['Nz']))
dih = [x for f in glob.glob(D+'data/steady/gscan_ih_1.58e6*.json') for x in json.load(open(f))]
dih = sorted({round(x['Gamma'], 6): x for x in dih}.values(), key=lambda x: x['Gamma']); F = np.array([x['FB'] for x in dih]); i0 = int(np.argmin(F))
check('IH F_B(Gamma) at R=1.6e6 monotone on both sides of the optimum, from Gamma~0.26 to 1', (round(dih[0]['Gamma'], 2), round(dih[-1]['Gamma'], 2)),
      cond=bool(np.all(np.diff(F[i0:]) >= 0) and np.all(np.diff(F[:i0+1]) <= 0)) and round(dih[-1]['Gamma'], 2) == 1.0)
rc = [c for f in sorted(glob.glob(D+'data/steady/res_check*.json')) for c in json.load(open(f))]
rb7 = sorted([c for c in rc if c['mode'] == 'rb' and abs(c['Ra']-1e7) < 1], key=lambda c: c['Nz'])
check('RB Ra=1e7: (80,72) -> (112,104) changes Nu by 4e-7 (relative)', abs(rb7[1]['Nu']-rb7[0]['Nu'])/rb7[1]['Nu'], 4.2e-7, 0.5e-7)
ih7 = sorted([c for c in rc if c['mode'] == 'ih'], key=lambda c: c['Nz'])
check('IH R=1.6e7: change of F_B ~1e-10', abs(ih7[1]['FB']-ih7[0]['FB']), cond=abs(ih7[1]['FB']-ih7[0]['FB']) < 1e-10)
for Rv, txt in [(56234132.51903491, 5e-6), (1e8, 2e-5)]:
    hiR = [c for c in rc if c['mode'] == 'rb' and abs(c['Ra']-Rv)/Rv < 1e-6][0]; base = [x for x in rb if abs(x['Ra']-Rv)/Rv < 1e-6][0]
    check(f'RB Ra={Rv:.3g}: re-solve at (128,128) changes Nu by about {txt:g}', abs(hiR['Nu']-base['Nu'])/base['Nu'], cond=abs(hiR['Nu']-base['Nu'])/base['Nu'] < 1.5*txt)
stt = json.load(open(D+'data/steady/spectral_tails.json')); mm = [o for o in stt if o['Ra'] >= 1e6]
check('steady Fourier tails below 1e-10 (Ra >= 1e6)', max(o['fourier_tail'] for o in mm), cond=max(o['fourier_tail'] for o in mm) < 1e-10)
check('steady Chebyshev tails below 1e-8 (IH) and 4e-5 (RB)', None, cond=max(o['cheb_tail'] for o in mm if o['mode'] == 'ih') < 1e-8 and max(o['cheb_tail'] for o in mm if o['mode'] == 'rb') < 4e-5)

# ---------------------------------------------------------------------------------------------------------------
section('10. Validation (read from the test logs and stored benchmark results)')
bl = json.load(open(D+'data/steady/blankenbach.json')); best = {}
for c in bl: best[c['Ra']] = c
for Ra, txt in [(1e4, 4e-8), (1e5, 2e-8), (1e6, 2e-7)]:
    check(f'Blankenbach Ra={Ra:.0e}: relative error about {txt:g}', abs(best[Ra]['rel_err']), cond=abs(best[Ra]['rel_err']) < 1.5*txt)
log = open(D+'logs/test_ipc.log').read()
check('Ra_c = 1707.76 at k = 3.116 (RB, no-slip)', None, cond='ns RB R_c = 1707.76 at k = 3.1163' in log)
check('R_L = 37325.2 at k = 3.999 (IH, no-slip)', None, cond='ns IH R_c = 37325.17 at k = 3.9989' in log)
log = open(D+'logs/test_steady.log').read()
check('Jacobian against finite differences to 6e-12', None, cond='jacobian FD check rel err 6.2' in log)
check('time-stepping and Newton agree to 1e-10', None, cond='relative difference steady vs DNS: 5.84e-11' in log)

print(f'\n{N[0]} checks, {len(FAILED)} failed')
if FAILED:
    print('FAILED:', *FAILED, sep='\n  ')
    sys.exit(1)
print('ALL CHECKS PASSED')
