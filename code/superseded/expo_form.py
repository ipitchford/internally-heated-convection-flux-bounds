"""Exponential-polynomial form of the slab-correction symbols, for rigorous large-kappa tail bounds.
With X = e^kappa, Y = e^{kappa zeta}, S = e^{kappa s}, each symbol is N/Dn with N a Laurent polynomial in
X, Y, S whose coefficients are polynomials in kappa, zeta, s, and Dn a polynomial in kappa and X.
Reports, for each monomial X^a Y^b S^c, the exponent a + b zeta + c s relative to the leading denominator
exponent, maximised over the parameter rectangle, minus the claimed decay rate -(2 - zeta - s)."""
import sympy as sp, pickle, os, itertools, json
H = os.path.dirname(os.path.abspath(__file__))
k, z, s = sp.symbols('kappa z s', positive=True)
X, Y, S = sp.symbols('X Y S', positive=True)
a_sf, b_sf = pickle.load(open(os.path.join(H, 'slab_corr_sf.pkl'), 'rb'))
a_ns, b_ns = pickle.load(open(os.path.join(H, 'slab_corr.pkl'), 'rb'))
c_sf = a_sf*sp.sinh(k*z) + b_sf*k*z*sp.cosh(k*z)
c_ns = a_ns*(k*z*sp.cosh(k*z) - sp.sinh(k*z)) + b_ns*k*z*sp.sinh(k*z)
def to_XYS(e):
    e = e.rewrite(sp.exp)
    e = e.subs(sp.exp(k*z), Y)
    e = sp.expand_power_exp(sp.powsimp(sp.expand(e)))
    e = e.replace(lambda t: isinstance(t, sp.exp), lambda t: sp.exp(sp.expand(t.args[0])))
    def conv(t):
        ex = sp.expand(t.args[0]); out = 1; rest = 0
        for term in sp.Add.make_args(ex):
            cz = term.coeff(k*z) if term.has(z) else 0
            cs = term.coeff(k*s) if term.has(s) else 0
            ck = sp.simplify((term - cz*k*z - cs*k*s)/k)
            out *= Y**cz * S**cs * X**ck
        return out
    return e.replace(lambda t: isinstance(t, sp.exp), conv)
res = {}
for name, cexpr, dz in [('ns', c_ns, 2), ('sf', c_sf, 1)]:
    mu = k**2*sp.diff(cexpr, z, dz)
    if name == 'sf':
        mu = mu + s/2*k*sp.exp(-k)
    e = sp.together(to_XYS(mu))
    num, den = sp.fraction(e)
    num = sp.expand(num); den = sp.expand(den)
    P = sp.Poly(num, X, Y, S, 1/X, 1/Y, 1/S) if False else None
    terms = sp.Add.make_args(num)
    dterms = sp.Add.make_args(den)
    def expo(t):
        dX = sp.degree(t, X) if t.has(X) else 0
        # handle negative powers
        pw = sp.powsimp(t).as_powers_dict()
        return pw.get(X, 0), pw.get(Y, 0), pw.get(S, 0)
    dexp = max((expo(t)[0] for t in dterms))
    worst = -1e9; worst_t = None
    for t in terms:
        aX, bY, cS = expo(t)
        # exponent of the term divided by the leading denominator X^dexp, plus the decay rate (2 - z - s):
        # must be <= 0 on 0 <= z <= 1/2, 0 <= s <= 1 for a termwise bound to show decay like e^{-kappa(2-z-s)}
        val = max(float(aX - dexp + bY*zz + cS*ss + 2 - zz - ss) for zz in (0, 0.5) for ss in (0, 1))
        if val > worst: worst, worst_t = val, (aX, bY, cS)
    print(name, 'numerator terms:', len(terms), ' denominator leading X power:', dexp,
          ' max over terms of [exponent + (2 - zeta - s)]:', worst, 'at monomial X^a Y^b S^c =', worst_t, flush=True)
    res[name] = dict(n_terms=len(terms), den_lead=int(dexp), worst=worst)
