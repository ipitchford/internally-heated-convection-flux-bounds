import sympy as sp, mpmath as mp
k, z, s, c, y, r, u = sp.symbols('kappa z s c y r u', positive=True)
# half-space symbol m = kappa^2 g_hs for z<s
a, b = s-z, s+z
m = ((1+k*a)*sp.exp(-k*a) - (1+k*b+2*k**2*z*s)*sp.exp(-k*b))/(4*k)
d2 = sp.diff(m, z, 2)
claim = sp.Rational(1,4)*(k**2*a*sp.exp(-k*a) - k*sp.exp(-k*a) - sp.exp(-k*b)*(k**2*(b-4*s) + 2*s*k**3*z - k))
print('d2m - claim =', sp.simplify(d2 - claim))
# representation m = (kappa/4)[int_a^b c e^{-kc} dc - 2 z s e^{-kb}]
rep = k/4*(sp.integrate(c*sp.exp(-k*c),(c,a,b)) - 2*z*s*sp.exp(-k*b))
print('m - rep =', sp.simplify(m - rep))
# nu_j^(1) = || F^{-1}[|k|^j e^{-|k|}] ||_{L^1(R)}, kernel = (1/pi) Re[ j!/(1-iy)^{j+1} ].
# |kernel| has corners at the kernel's zeros, so the quadrature is split there (an unsplit quadrature loses
# about five digits for j = 3).
mp.mp.dps = 30
def _split(f, roots):
    pts = [0] + sorted(roots) + [mp.inf]
    return sum(abs(mp.quad(f, [pts[i], pts[i+1]])) for i in range(len(pts)-1))
def nu1(j):
    num = sp.expand(sp.re(sp.expand((1+sp.I*y)**(j+1))))
    roots = [float(x) for x in sp.Poly(num, y).nroots() if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0]
    f = sp.lambdify(y, sp.factorial(j)*num/(1+y**2)**(j+1)/sp.pi, 'mpmath'); return 2*_split(f, roots)
# nu_j^(2): kernel(r) = (1/2pi) int k^{j+1} e^{-k} J0(kr) dk = (1/2pi) (-d/dc)^{j+1} (c^2+r^2)^(-1/2) at c=1
def nu2(j):
    expr = sp.simplify(((-1)**(j+1)*sp.diff((c**2+r**2)**sp.Rational(-1, 2), c, j+1)).subs(c, 1))
    roots = [float(x) for x in sp.Poly(sp.expand(sp.numer(sp.together(expr))), r).nroots() if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0]
    f = sp.lambdify(r, expr*r, 'mpmath'); return _split(f, roots)   # (1/2pi)*2pi
for j in range(4):
    print('j=%d  nu^(1)=%.10f  nu^(2)=%.10f' % (j, nu1(j), nu2(j)))
print('checks: 2/pi=%.10f, 3sqrt3/(2pi)=%.10f, 6/pi=%.10f, 4/(3sqrt3)=%.10f' % (2/mp.pi, 3*mp.sqrt(3)/(2*mp.pi), 6/mp.pi, 4/(3*mp.sqrt(3))))
