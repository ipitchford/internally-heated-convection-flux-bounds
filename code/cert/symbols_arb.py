"""Ball-arithmetic (arb) evaluation of the slab-correction remainder symbols and their first two
kappa-derivatives, for certified bounds (Appendix B of the paper).

  ns:  r(kappa; zeta, s) = kappa^2 d^2/dzeta^2 c_ns,          c = g_slab - g_halfspace (clamped walls)
  sf:  r(kappa; zeta, s) = kappa^2 d/dzeta c_sf + (s/2) kappa e^{-kappa}   (simply supported walls)

Two representations, both free of removable singularities:
  * small/moderate kappa: fundamental system u_0..u_3 of (d^2 - kappa^2)^2 at zeta = 0 (entire power
    series in kappa^2 with factorial decay) and half-space data at zeta = 1 written as convergent series;
  * large kappa: the decaying basis e^{-kz}, kz e^{-kz}, e^{-k(1-z)}, k(1-z) e^{-k(1-z)} (4x4 solve).
Every function returns an arb_series of length 3 in kappa: (r, r', r''/2) enclosures valid for every
kappa in the kappa-ball and every (zeta, s) in the parameter balls."""
from flint import arb, arb_series, ctx
ctx.prec = 192
ctx.cap = 3

def ser(x):          # constant series
    return arb_series([x])

def kser(kball):     # kappa = kball + t
    return arb_series([kball, arb(1)])

def _fact_table(N):
    f = [arb(1)]
    for i in range(1, N+1): f.append(f[-1]*i)
    return f
FACT = _fact_table(200)

def _tail_ball(first_omitted_bound):
    """ball [-B, B] used to absorb a series tail whose absolute value is <= B"""
    return arb(0, first_omitted_bound)

def powser(coeff, x2, N, ratio_bound, deriv_bound_factor=1):
    """sum_{k=0}^{N} coeff(k) * x2^k for an arb_series x2 (= kappa^2 * const); tail bounded by
    2*|coeff(N+1) x2max^(N+1)| times deriv_bound_factor when the term ratio is < 1/2."""
    out = ser(0); p = ser(1)
    for k in range(N+1):
        out = out + p*coeff(k); p = p*x2
    return out

# ---------- fundamental system at zeta = 0, as series in kappa (kappa enters through x = kappa*zeta) ----------
def u_funcs(K, zeta, N=48):
    """u_j(zeta), u_j'(zeta), u_j''(zeta), j = 0..3, as arb_series in kappa (K = kappa series).
    u0 = sum (1-k) x^{2k}/(2k)!,                u1 = zeta sum (1-k) x^{2k}/(2k+1)!,
    u2 = (zeta^2/2) sum x^{2k}/(2k+1)!,          u3 = zeta^3 sum_{k>=1} k x^{2k-2}/(2k+1)!,   x = kappa zeta.
    Derivatives in zeta are taken termwise from the monomial form c_k kappa^{2k} zeta^{m_k}."""
    K2 = K*K
    res = {}
    # monomial coefficient tables: term k contributes a_k * kappa^{2k} * zeta^{e_k}
    specs = {0: (lambda k: arb(1-k)/FACT[2*k], lambda k: 2*k, 0),
             1: (lambda k: arb(1-k)/FACT[2*k+1], lambda k: 2*k+1, 0),
             2: (lambda k: arb(1)/(2*FACT[2*k+1]), lambda k: 2*k+2, 0),
             3: (lambda k: arb(k+1)/FACT[2*k+3], lambda k: 2*k+3, 0)}   # u3: reindexed k -> k+1
    for j, (a, e, _) in specs.items():
        vals = []
        for dz in range(3):
            tot = ser(0); Kp = ser(1)
            for k in range(N+1):
                ee = e(k)
                if ee - dz >= 0:
                    fall = 1
                    for t in range(dz): fall *= (ee - t)
                    tot = tot + Kp*(a(k)*fall*zeta**(ee-dz))
                Kp = Kp*K2
            vals.append(tot)
        res[j] = vals
    return res

def series_tail_check(kmax, zmax, N):
    """|term_k| <= kmax^{2k} zmax^{2k-2}/(2k-1)! up to polynomial factors; returns a crude bound on the
    omitted tail of every u-series (and its first two kappa- and zeta-derivatives)."""
    x = kmax*zmax if zmax > 0 else kmax
    b = arb(0); t = arb(1)
    for k in range(N+1, N+60):
        term = arb(kmax)**(2*k) * arb(max(zmax, 1))**(2*k+3) / FACT[2*k-1] * (2*k+3)**4
        b = b + term
    return b

# ---------- half-space data at zeta = 1 ----------
def hs_data_ns(K, s, N=80):
    """g_hs(1,s) and d_z g_hs(1,s) (clamped half-space), a = 1-s, b = 1+s:
       g  = e^{-kb}/4 * sum_{k>=3} [ (2s)^k/k! + a (2s)^{k-1}/(k-1)! ] kappa^{k-3}
       gz = e^{-kb}/4 * [ 2 s^2 - a sum_{k>=2} (2s)^k kappa^{k-1}/k! ]"""
    a, b = 1 - s, 1 + s
    eb = (-(K*b)).exp() if hasattr(K*b, 'exp') else None
    eb = arb_series.exp(-(K*b))
    S = ser(0); Kp = ser(1)
    for k in range(3, N+1):
        S = S + Kp*((2*s)**k/FACT[k] + a*(2*s)**(k-1)/FACT[k-1]); Kp = Kp*K
    T = ser(0); Kp = K
    for k in range(2, N+1):                      # sum_{k>=2} (2s)^k kappa^{k-1}/k!
        T = T + Kp*((2*s)**k/FACT[k]); Kp = Kp*K
    return eb*S/4, eb*(ser(2*s*s) - T*a)/4

def hs_data_sf(K, s, N=80):
    """stress-free half-space G0(1-s) - G0(1+s) = s e^{-kb}/(2 kappa) + g_reg, returned as
       (kappa^2 * pole part = s kappa e^{-kb}/2,  kappa^2 g_reg,  kappa^2 d_z^2 g_hs(1,s))."""
    a, b = 1 - s, 1 + s
    eb = arb_series.exp(-(K*b))
    S = ser(0); Kp = ser(1)
    for k in range(3, N+1):
        S = S + Kp*((2*s)**k/FACT[k] + a*(2*s)**(k-1)/FACT[k-1]); Kp = Kp*K
    T = ser(0); Kp = ser(1)
    for k in range(2, N+1):
        T = T + Kp*(a*(2*s)**(k-1)/FACT[k-1] - (2*s)**k/FACT[k]); Kp = Kp*K
    K2 = K*K
    return K*eb*(s/2), K2*eb*S/4, K2*eb*(ser(-4*s) + T*K)/4   # last: kappa^2 * g_zz
    
def solve2(m11, m12, m21, m22, r1, r2):
    det = m11*m22 - m12*m21
    return (r1*m22 - r2*m12)/det, (m11*r2 - m21*r1)/det

def r_small(kball, zeta, s, bc, Nu=48, Nh=80):
    """remainder symbol series for kappa in kball (small/moderate kappa), via the fundamental system"""
    K = kser(kball)
    U1 = u_funcs(K, arb(1), Nu); Uz = u_funcs(K, zeta, Nu)
    if bc == 'ns':
        g, gz = hs_data_ns(K, s, Nh)
        # c = p u2 + q u3 with c(1) = -g, c'(1) = -gz;  r = kappa^2 c''(zeta)
        p, q = solve2(U1[2][0], U1[3][0], U1[2][1], U1[3][1], -g, -gz)
        return K*K*(p*Uz[2][2] + q*Uz[3][2])
    # sf: c = p u1 + q u3 with c(1) = -g_hs(1,s), c''(1) = -g_zz(1,s); kappa^2 c = kappa^2 * (...)
    pole, greg, gzz = hs_data_sf(K, s, Nh)
    # kappa^2 c = -(pole) v0 - greg v0 - gzz v2, with v0, v2 the solutions with (v(1), v''(1)) = (1,0), (0,1)
    p0, q0 = solve2(U1[1][0], U1[3][0], U1[1][2], U1[3][2], ser(1), ser(0))
    p2, q2 = solve2(U1[1][0], U1[3][0], U1[1][2], U1[3][2], ser(0), ser(1))
    v0z = p0*Uz[1][1] + q0*Uz[3][1]; v2z = p2*Uz[1][1] + q2*Uz[3][1]
    mu = -(pole + greg)*v0z - gzz*v2z
    return mu + K*arb_series.exp(-K)*(s/2)

# ---------- large kappa: decaying basis ----------
def _basis(K, z, d):
    E1 = arb_series.exp(-(K*z)); E3 = arb_series.exp(-(K*(1-z))); y = 1 - z
    if d == 0: return [E1, K*z*E1, E3, K*y*E3]
    if d == 1: return [-K*E1, K*(ser(1)-K*z)*E1, K*E3, K*(K*y-1)*E3]
    if d == 2: return [K*K*E1, K*K*(K*z-2)*E1, K*K*E3, K*K*(K*y-2)*E3]

def _solve4(M, rhs):
    """Gaussian elimination with no pivoting (the matrix is diagonally dominant for kappa >= 3)"""
    A = [row[:] + [rhs[i]] for i, row in enumerate(M)]
    n = 4
    for c in range(n):
        for r_ in range(c+1, n):
            f = A[r_][c]/A[c][c]
            A[r_] = [A[r_][j] - f*A[c][j] for j in range(n+1)]
    x = [None]*n
    for r_ in range(n-1, -1, -1):
        acc = A[r_][n]
        for j in range(r_+1, n): acc = acc - A[r_][j]*x[j]
        x[r_] = acc/A[r_][r_]
    return x

def r_large(kball, zeta, s, bc):
    K = kser(kball)
    rows = (0, 1) if bc == 'ns' else (0, 2)
    z0, z1 = arb(0), arb(1)
    # order rows so that the matrix is close to block-triangular with nonzero pivots
    M = [_basis(K, z0, rows[0]), _basis(K, z0, rows[1]), _basis(K, z1, rows[0]), _basis(K, z1, rows[1])]
    a, b = 1 - s, 1 + s
    ea, eb = arb_series.exp(-(K*a)), arb_series.exp(-(K*b))
    if bc == 'ns':
        g1 = ((ser(1)+K*a)*ea - (ser(1)+K*b+K*K*(2*s))*eb)/(K*K*K*4)
        dg1 = (-(K*K*a)*ea + K*K*(ser(1-s)+K*(2*s))*eb)/(K*K*K*4)
        rhs = [ser(0), ser(0), -g1, -dg1]
        cs = _solve4(M, rhs)
        B = _basis(K, zeta, 2)
        return K*K*(cs[0]*B[0] + cs[1]*B[1] + cs[2]*B[2] + cs[3]*B[3])
    G = lambda x, e: ((ser(1)+K*x)*e/(K*K*K*4), (K*x - 1)*e/(K*4))   # G0 and G0'' at distance x
    Gm, Gp = G(a, ea), G(b, eb)
    rhs = [ser(0), ser(0), -(Gm[0]-Gp[0]), -(Gm[1]-Gp[1])]
    cs = _solve4(M, rhs)
    B = _basis(K, zeta, 1)
    mu = K*K*(cs[0]*B[0] + cs[1]*B[1] + cs[2]*B[2] + cs[3]*B[3])
    return mu + K*arb_series.exp(-K)*(s/2)
