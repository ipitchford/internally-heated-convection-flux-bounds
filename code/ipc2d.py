"""ipc2d: two-dimensional infinite-Prandtl-number Boussinesq convection in a horizontally periodic slab.

    -Lap u + grad p = R T e_z,  div u = 0,          (Stokes, velocity slaved to T)
    dT/dt + div(u T) = Lap T + Q                     (Q = 1: internal heating; Q = 0: Rayleigh-Benard)

z in [0,1] (Chebyshev-Gauss-Lobatto collocation), x in [0, Gamma) (Fourier, 3/2-rule dealiasing).
Velocity: for each horizontal wavenumber k the vertical velocity solves (D^2-k^2)^2 w = R k^2 T with
no-slip (w = Dw = 0) or stress-free (w = D^2 w = 0) walls; this linear map is precomputed with an
influence-matrix method built from two well-conditioned Dirichlet Helmholtz solves. u = (i/k) Dw.
Time stepping: SBDF2 (implicit diffusion, explicit advection), adaptive step via CFL with SBDF1 restart.
"""
import numpy as np, time, json

def cheb(N):
    """Trefethen's Chebyshev differentiation matrix on x_j = cos(pi j/N), j=0..N."""
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    Dm = np.outer(c, 1/c)/(dX + np.eye(N+1))
    Dm = Dm - np.diag(Dm.sum(1))
    return Dm, x

def clenshaw_curtis(N):
    """quadrature weights on x_j=cos(pi j/N) for int_{-1}^{1}"""
    theta = np.pi*np.arange(N+1)/N
    w = np.zeros(N+1); v = np.ones(N-1)
    if N % 2 == 0:
        w[0] = w[N] = 1/(N**2-1)
        for k in range(1, N//2): v -= 2*np.cos(2*k*theta[1:-1])/(4*k*k-1)
        v -= np.cos(N*theta[1:-1])/(N**2-1)
    else:
        w[0] = w[N] = 1/N**2
        for k in range(1, (N-1)//2+1): v -= 2*np.cos(2*k*theta[1:-1])/(4*k*k-1)
    w[1:-1] = 2*v/N
    return w

class IPC2D:
    def __init__(self, R, Gamma=4.0, Nx=256, Nz=64, bc='ns', Q=1.0, rb=False, seed=0):
        self.R, self.G, self.Nx, self.Nz, self.bc, self.Q, self.rb = R, Gamma, Nx, Nz, bc, Q, rb
        Dxi, xi = cheb(Nz)
        self.z = (1+xi)/2                      # z[0]=1 (top), z[Nz]=0 (bottom)
        self.D = 2*Dxi; self.D2 = self.D @ self.D
        self.wz = clenshaw_curtis(Nz)/2        # weights for int_0^1 dz
        self.k = 2*np.pi*np.arange(Nx//2+1)/Gamma
        self.Nk = len(self.k); self.M = 3*Nx//2
        self.x = np.arange(Nx)*Gamma/Nx
        self._build_velocity_maps()
        rng = np.random.default_rng(seed)
        Tphys = np.outer(self.base_profile(), np.ones(Nx)) + 1e-3*rng.standard_normal((Nz+1, Nx))*np.outer(np.sin(np.pi*self.z), np.ones(Nx))
        self.Th = np.fft.rfft(Tphys, axis=1)/Nx   # normalised coefficients
        self.Th[:, -1] = 0.0
        self.apply_bc(self.Th)
        self.t = 0.0; self.dt = None; self.prev = None

    def base_profile(self):
        return (1-self.z) if self.rb else self.z*(1-self.z)/2

    def apply_bc(self, Th):
        Th[0, :] = 0.0; Th[-1, :] = 0.0
        if self.rb: Th[-1, 0] = 1.0

    def _build_velocity_maps(self):
        N = self.Nz; I = np.eye(N+1)
        A = np.zeros((self.Nk, N+1, N+1)); AU = np.zeros((self.Nk, N+1, N+1))
        Z = I.copy(); Z[0, 0] = Z[N, N] = 0.0   # zero the boundary rows of the forcing
        for j, k in enumerate(self.k):
            if k == 0: continue
            L = self.D2 - k*k*I
            LD = L.copy(); LD[0, :] = 0; LD[0, 0] = 1; LD[N, :] = 0; LD[N, N] = 1
            LDi = np.linalg.inv(LD)
            P = LDi @ Z @ LDi @ Z                    # w_p from f (omega_p = 0 at walls, w_p = 0 at walls)
            if self.bc == 'ns':
                et = np.zeros(N+1); et[0] = 1; eb = np.zeros(N+1); eb[N] = 1
                wt = LDi @ Z @ (LDi @ et); wb = LDi @ Z @ (LDi @ eb)
                C = np.array([[self.D[0] @ wt, self.D[0] @ wb], [self.D[N] @ wt, self.D[N] @ wb]])
                B = np.vstack([self.D[0] @ P, self.D[N] @ P])
                coef = -np.linalg.solve(C, B)          # (2, N+1): alpha, beta as functionals of f
                Aw = P + np.outer(wt, coef[0]) + np.outer(wb, coef[1])
            else:
                Aw = P
            A[j] = Aw*(self.R*k*k)                  # w_hat = A T_hat
            AU[j] = (self.D @ Aw)*(self.R*k*k)/k    # u_hat = i * AU T_hat
        self.Aw, self.AU = A, AU

    @staticmethod
    def bmv(A, Th):
        """batched real-matrix x complex-vector: A (Nk,N,N), Th (N,Nk) complex -> (N,Nk)"""
        X = np.stack([Th.real.T, Th.imag.T], -1)
        Y = np.matmul(A, X)
        return (Y[..., 0] + 1j*Y[..., 1]).T

    def velocity_hat(self, Th):
        wh = self.bmv(self.Aw, Th)
        uh = 1j*self.bmv(self.AU, Th)
        return uh, wh

    def to_phys(self, fh):
        pad = np.zeros((fh.shape[0], self.M//2+1), complex); pad[:, :self.Nk] = fh
        return np.fft.irfft(pad, n=self.M, axis=1)*self.M

    def to_spec(self, f):
        out = np.fft.rfft(f, axis=1)[:, :self.Nk]/self.M
        out[:, -1] = 0.0
        return out

    def nonlinear(self, Th):
        uh, wh = self.velocity_hat(Th)
        T, u, w = self.to_phys(Th), self.to_phys(uh), self.to_phys(wh)
        uTh, wTh = self.to_spec(u*T), self.to_spec(w*T)
        Nh = -(1j*self.k[None, :]*uTh + self.D @ wTh)
        Nh[:, 0] += self.Q
        return Nh, (T, u, w)

    def _factor(self, sigma):
        N = self.Nz; I = np.eye(N+1)
        Hi = np.zeros((self.Nk, N+1, N+1))
        for j, k in enumerate(self.k):
            H = sigma*I - self.D2 + k*k*I
            H[0, :] = 0; H[0, 0] = 1; H[N, :] = 0; H[N, N] = 1
            Hi[j] = np.linalg.inv(H)
        return Hi

    def cfl(self, u, w):
        dx = self.G/self.M
        dz = np.abs(np.gradient(self.z))
        return np.max(np.abs(u))/dx + np.max(np.abs(w)/dz[:, None])

    def step(self, dt):
        Nh, (T, u, w) = self.nonlinear(self.Th)
        if self.dt != dt or self.prev is None:
            # (re)start with SBDF1
            if not hasattr(self, '_Hi1') or self._dt1 != dt:
                self._Hi1 = self._factor(1/dt); self._dt1 = dt
            rhs = self.Th/dt + Nh
            self.apply_bc(rhs)
            Tnew = self.bmv(self._Hi1, rhs)
            self.prev = (self.Th, Nh); self.dt = dt
            if not hasattr(self, '_Hi2') or self._dt2 != dt:
                self._Hi2 = self._factor(1.5/dt); self._dt2 = dt
        else:
            Tm, Nm = self.prev
            rhs = (2*self.Th - 0.5*Tm)/dt + 2*Nh - Nm
            self.apply_bc(rhs)
            Tnew = self.bmv(self._Hi2, rhs)
            self.prev = (self.Th, Nh)
        self.Th = Tnew; self.t += dt
        return T, u, w

    # ---------- diagnostics ----------
    def profiles(self, T, w):
        Tbar = T.mean(1); J = (w*T).mean(1); Wmax = np.abs(w).max(1)
        return Tbar, J, Wmax

    def fluxes(self, Th):
        T0 = Th[:, 0].real
        dT = self.D @ T0
        FB = dT[-1]; FT = -dT[0]   # F_B = d<T>/dz at z=0 (heat leaving downward), F_T = -d<T>/dz at z=1
        return FB, FT, np.sum(self.wz*T0)



    def run(self, t_end, dt0=1e-5, cfl_target=0.35, t_avg_start=None, diag_every=1, log_every=500,
            dt_max=2e-3, verbose=True, max_wall=None, snapshot_every=None, n_blocks=8):
        """Integrate to t_end.

        Diagnostics are evaluated at the time level t_n of the fields used in each step: the temperature
        coefficients, the velocity and the boundary fluxes all come from the same state. Time averages over
        [t_avg_start, t_end] are physical-time integrals (trapezoidal rule over every step) divided by the
        window length, so a variable time step does not bias them. The window is also split into n_blocks
        equal sub-windows whose integrals are returned for window-stability estimates."""
        dt = dt0; hist = []; n = 0; tw = time.time(); acc = None; last_change = 0; snaps = []
        prev_sample = None

        def sample(tn, Thn, T, w):
            FB, FT, Tm = self.fluxes(Thn)
            Tbar, J, Wmax = self.profiles(T, w)
            return dict(t=tn, FB=FB, FT=FT, Tm=Tm, wT=float(np.sum(self.wz*J)), Tbar=Tbar, J=J, Wmax=Wmax,
                        Tmax=float(T.max()), Tmin=float(T.min()), wmax=float(np.abs(w).max()))

        while self.t < t_end:
            tn, Thn = self.t, self.Th
            T, u, w = self.step(dt)          # fields at time tn, computed from the state Thn
            n += 1
            c = self.cfl(u, w)*dt
            averaging = t_avg_start is not None and tn >= t_avg_start
            if n % diag_every == 0 or averaging:
                q = sample(tn, Thn, T, w)
                if n % diag_every == 0:
                    hist.append((tn, q['FB'], q['FT'], q['Tm'], q['wT'], q['Tmax'], q['Tmin'], q['wmax'], dt))
                if averaging:
                    if acc is None:
                        L = (t_end - tn)/n_blocks
                        acc = dict(t0=tn, t1=tn, Tm0=q['Tm'], Tm1=q['Tm'], Wmax=q['Wmax'].copy(), Tmax=q['Tmax'],
                                   Tmin=q['Tmin'], nsteps=0, L=L,
                                   blocks=[dict(ta=tn+b*L, tb=tn+(b+1)*L, W=0.0, FB=0.0, FT=0.0, wT=0.0, Tm=0.0,
                                                Tbar=0*q['Tbar'], J=0*q['J']) for b in range(n_blocks)])
                    else:
                        h = tn - prev_sample['t']; tm = 0.5*(tn + prev_sample['t'])
                        b = min(int((tm - acc['t0'])/acc['L']), n_blocks-1)
                        B = acc['blocks'][b]; B['W'] += h
                        for key in ('FB', 'FT', 'wT', 'Tm'):
                            B[key] += 0.5*h*(q[key] + prev_sample[key])
                        B['Tbar'] += 0.5*h*(q['Tbar'] + prev_sample['Tbar'])
                        B['J'] += 0.5*h*(q['J'] + prev_sample['J'])
                    acc['Wmax'] = np.maximum(acc['Wmax'], q['Wmax']); acc['Tmax'] = max(acc['Tmax'], q['Tmax'])
                    acc['Tmin'] = min(acc['Tmin'], q['Tmin']); acc['t1'] = tn; acc['Tm1'] = q['Tm']
                    acc['nsteps'] += 1
                    prev_sample = q
                if verbose and n % log_every == 0:
                    print(f"t={tn:.5f} dt={dt:.2e} FB={q['FB']:.5f} FT={q['FT']:.5f} <wT>={q['wT']:.5f} "
                          f"<T>={q['Tm']:.5f} Tmax={q['Tmax']:.5f} Tmin={q['Tmin']:.2e} |w|max={q['wmax']:.2f} "
                          f"cfl={c:.2f} wall={time.time()-tw:.0f}s", flush=True)
                if not np.isfinite(q['FB']):
                    raise RuntimeError('blow-up')
            if snapshot_every and n % snapshot_every == 0:
                snaps.append((self.t, self.to_phys(self.Th), self.to_phys(self.velocity_hat(self.Th)[1])))
            if c > 1.3*cfl_target:                       # reduce immediately
                dt = dt*cfl_target/c; last_change = n
            elif c < 0.7*cfl_target and n - last_change > 50 and dt < dt_max:
                dt = min(dt*1.15, dt_max, dt*cfl_target/max(c, 1e-12)); last_change = n
            if max_wall is not None and time.time()-tw > max_wall:
                break
        return np.array(hist), acc, snaps
