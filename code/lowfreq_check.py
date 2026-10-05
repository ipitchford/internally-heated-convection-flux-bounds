"""Low-frequency structure of the slab-correction symbols (Lemma 3.6).
mu_sf = kappa^2 d/dzeta c_sf(zeta,s), mu_ns = kappa^2 d^2/dzeta^2 c_ns(zeta,s), c = g_slab - g_halfspace.
Prints the Laurent/Taylor expansions in kappa at kappa = 0 with exact coefficients."""
import sympy as sp, pickle, os
H = os.path.dirname(os.path.abspath(__file__))
k, z, s = sp.symbols('kappa z s', positive=True)
a_sf, b_sf = pickle.load(open(os.path.join(H, 'slab_corr_sf.pkl'), 'rb'))
a_ns, b_ns = pickle.load(open(os.path.join(H, 'slab_corr.pkl'), 'rb'))
c_sf = a_sf*sp.sinh(k*z) + b_sf*k*z*sp.cosh(k*z)
c_ns = a_ns*(k*z*sp.cosh(k*z) - sp.sinh(k*z)) + b_ns*k*z*sp.sinh(k*z)
mu_sf = k**2*sp.diff(c_sf, z); mu_ns = k**2*sp.diff(c_ns, z, 2)
out = {}
for name, mu in [('sf', mu_sf), ('ns', mu_ns)]:
    ser = sp.series(mu, k, 0, 5).removeO()
    coeffs = [sp.factor(sp.simplify(ser.coeff(k, j))) for j in range(0, 5)]
    out[name] = coeffs
    print(name, 'mu = sum c_j kappa^j with'); [print('   c_%d =' % j, c) for j, c in enumerate(coeffs)]
r = sp.series(mu_sf + s/2*k*sp.exp(-k), k, 0, 4).removeO()
print('sf remainder r = mu + (s/2) kappa e^{-kappa}: kappa^0, kappa^1, kappa^2 coefficients:',
      [sp.factor(sp.simplify(r.coeff(k, j))) for j in range(0, 3)])
# endpoint identity zeta = 0, s = 1
print('mu_sf(kappa; 0, 1) + kappa e^{-kappa}/2 =', sp.simplify(mu_sf.subs({z: 0, s: 1}) + k*sp.exp(-k)/2))
