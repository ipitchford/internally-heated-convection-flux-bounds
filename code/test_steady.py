import numpy as np, time
from steady2d import Steady
from ipc2d import IPC2D
np.random.seed(1)
st = Steady(R=3000, Gamma=2*np.pi/3.117, Nx=32, Nz=32, bc='ns', rb=True)
z = st.s.z
X = np.zeros(st.shape); X[:, 0] = 1-z; X[:, 1] = 0.1*np.sin(np.pi*z)
# Jacobian check by finite differences
V = np.random.randn(*st.shape); V[0] = 0; V[-1] = 0
eps = 1e-6
fd = (st.residual(X+eps*V) - st.residual(X-eps*V))/(2*eps)
ja = st.jac_apply(X, V[None])[0]
print('jacobian FD check rel err', np.linalg.norm(fd-ja)/np.linalg.norm(ja))
t0 = time.time(); X, ok = st.newton(X); print('converged', ok, '%.1fs'%(time.time()-t0))
d = st.diagnostics(X); print('steady Nu =', d['Nu'], 'FT=', d['FT'], 'Tmax', d['Tmax'])
# DNS check: time-integrate the same configuration to steady state
s = IPC2D(R=3000, Gamma=2*np.pi/3.117, Nx=32, Nz=32, bc='ns', Q=0.0, rb=True)
s.Th[:, 1] += 0.05*np.sin(np.pi*s.z)
h, acc, _ = s.run(t_end=6.0, dt0=1e-3, dt_max=5e-3, verbose=False)
print('DNS final Nu =', -h[-1,1], ' FT =', h[-1,2])
print('relative difference steady vs DNS: %.2e' % abs(-h[-1,1]/d['Nu']-1))
