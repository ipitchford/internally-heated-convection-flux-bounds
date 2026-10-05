"""Newton solver for steady 2D infinite-Prandtl convection (internally heated or Rayleigh-Benard).

Mirror-symmetric rolls: T(-x,z) = T(x,z), so the Fourier coefficients T_hat_j(z) are real.
Unknowns X[m, j] = T_hat_j(z_m), m = 0..Nz (Chebyshev-Lobatto, z_0 = 1 top), j = 0..Nk-2 (Nyquist dropped).
Residual  G = (D^2 - k_j^2) X + N(X)  in the interior rows; boundary rows impose the wall temperatures.
The Jacobian is assembled densely by applying the exact linearisation to batches of unit vectors.
"""
import numpy as np, scipy.linalg as sla, time
from ipc2d import IPC2D

class Steady:
    def __init__(self, R, Gamma, Nx, Nz, bc='ns', rb=False):
        self.s = IPC2D(R=R, Gamma=Gamma, Nx=Nx, Nz=Nz, bc=bc, Q=0.0 if rb else 1.0, rb=rb)
        s = self.s
        self.K = s.Nk - 1                     # modes kept
        self.N = s.Nz
        self.kk = s.k[:self.K]
        self.shape = (self.N+1, self.K)
        self.n = (self.N+1)*self.K

    # ----- batched transforms (last axis = modes) -----
    def phys(self, Fh):                     # Fh (..., N+1, K) complex -> (..., N+1, M)
        s = self.s; pad = np.zeros(Fh.shape[:-1] + (s.M//2+1,), complex); pad[..., :self.K] = Fh
        return np.fft.irfft(pad, n=s.M, axis=-1)*s.M

    def spec(self, f):
        return (np.fft.rfft(f, axis=-1)[..., :self.K]/self.s.M)

    def vel(self, Xh):                      # Xh (B, N+1, K) real/complex -> uh, wh
        s = self.s
        A = s.Aw[:self.K]; AU = s.AU[:self.K]
        Xt = np.moveaxis(Xh, 0, -1)          # (N+1, K, B)
        Xt = np.transpose(Xt, (1, 0, 2))     # (K, N+1, B)
        wh = np.matmul(A, Xt); uh = 1j*np.matmul(AU, Xt)
        back = lambda Y: np.transpose(Y, (2, 1, 0))  # (B, N+1, K)
        return back(uh), back(wh)

    def div_flux(self, uh1, wh1, Th2):
        """spectral coefficients of div(u1 T2) for batches"""
        u, w, T = self.phys(uh1), self.phys(wh1), self.phys(Th2)
        uT, wT = self.spec(u*T), self.spec(w*T)
        return 1j*self.kk*uT + np.einsum('ij,bjk->bik', self.s.D, wT)

    def residual(self, X):
        s = self.s; Xb = X[None].astype(complex)
        uh, wh = self.vel(Xb)
        G = (np.einsum('ij,bjk->bik', s.D2, Xb) - self.kk**2*Xb - self.div_flux(uh, wh, Xb))[0].real
        G[:, 0] += s.Q
        G[0, :] = X[0, :]; G[-1, :] = X[-1, :]
        if s.rb: G[-1, 0] -= 1.0
        return G

    def jac_apply(self, X, V):              # V (B, N+1, K) real
        s = self.s; Xb = X[None].astype(complex); Vc = V.astype(complex)
        uT, wT = self.vel(Xb); uV, wV = self.vel(Vc)
        lin = np.einsum('ij,bjk->bik', s.D2, Vc) - self.kk**2*Vc
        adv = self.div_flux(np.broadcast_to(uT, Vc.shape), np.broadcast_to(wT, Vc.shape), Vc) + self.div_flux(uV, wV, np.broadcast_to(Xb, Vc.shape))
        J = (lin - adv).real
        J[:, 0, :] = V[:, 0, :]; J[:, -1, :] = V[:, -1, :]
        return J

    def jacobian(self, X, batch=128):
        n = self.n; Jm = np.empty((n, n))
        for c0 in range(0, n, batch):
            c1 = min(n, c0+batch); B = c1-c0
            V = np.zeros((B, n)); V[np.arange(B), np.arange(c0, c1)] = 1.0
            Jm[:, c0:c1] = self.jac_apply(X, V.reshape(B, *self.shape)).reshape(B, n).T
        return Jm

    def newton(self, X, tol=1e-10, maxit=30, verbose=True, damping=True):
        for it in range(maxit):
            G = self.residual(X); r = np.linalg.norm(G)/np.sqrt(self.n)
            if verbose: print(f'  newton it {it}: |G| = {r:.3e}', flush=True)
            if r < tol: return X, True
            Jm = self.jacobian(X)
            dX = -sla.solve(Jm, G.ravel(), check_finite=False).reshape(self.shape)
            lam = 1.0
            if damping:
                while lam > 1e-3:
                    Xn = X + lam*dX
                    if np.linalg.norm(self.residual(Xn))/np.sqrt(self.n) < r*(1-0.3*lam): break
                    lam *= 0.5
            X = X + lam*dX
        G = self.residual(X); return X, np.linalg.norm(G)/np.sqrt(self.n) < tol*100

    # ----- diagnostics -----
    def diagnostics(self, X):
        s = self.s; T0 = X[:, 0]; dT = s.D @ T0
        FB, FT = dT[-1], -dT[0]
        Xb = X[None].astype(complex); uh, wh = self.vel(Xb)
        T, w = self.phys(Xb)[0], self.phys(wh)[0]
        J = (w*T).mean(1)
        return dict(FB=FB, FT=FT, Tmean=float(np.sum(s.wz*T0)), wT=float(np.sum(s.wz*J)),
                    Tmax=float(T.max()), Tmin=float(T.min()), wmax=float(np.abs(w).max()), Nu=-FB)

    @staticmethod
    def interp(Xold, Nz_old, Nz_new, K_new):
        """spectral interpolation of a state to a new Chebyshev grid and mode count"""
        from numpy.polynomial import chebyshev as C
        xo = np.cos(np.pi*np.arange(Nz_old+1)/Nz_old); xn = np.cos(np.pi*np.arange(Nz_new+1)/Nz_new)
        Xn = np.zeros((Nz_new+1, K_new)); K = min(K_new, Xold.shape[1])
        for j in range(K):
            coef = C.chebfit(xo, Xold[:, j], Nz_old)
            Xn[:, j] = C.chebval(xn, coef)
        return Xn
