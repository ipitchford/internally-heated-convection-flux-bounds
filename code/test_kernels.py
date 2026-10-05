from kernels import *
# consistency with earlier no-slip numbers, and symmetry K_{z,s} vs K_{s,z}
print('ns n=1', L1_norm(0.1, 0.4, 1, 'ns'), '(expect 7.3130e-03)')
print('sym ns', L1_norm(0.05, 0.2, 1, 'ns'), L1_norm(0.2, 0.05, 1, 'ns'))
print('sym sf', L1_norm(0.05, 0.2, 1, 'sf'), L1_norm(0.2, 0.05, 1, 'sf'))
# check sf kernel against sine series directly: K(y) = sum_j sin sin (1 - j pi |y|) e^{-j pi|y|}/(2 j pi)
z, s = 0.15, 0.45
Kc = kernel_from_symbol(lambda k: corr_symbol(k, z, s, 'sf', 0), 1, 2-z-s)
for y in [0.05, 0.3, 1.0, 2.5]:
    j = np.arange(1, 400001); al = j*np.pi
    ser = np.sum(np.sin(al*z)*np.sin(al*s)*(1-al*y)*np.exp(-al*y)/(2*al))
    print(y, K_hs(np.array([y]), z, s, 1, 'sf')[0] + Kc(y), ser)
# derivative symbols: finite-difference check
k = np.array([0.5, 2.0, 7.0]); h = 1e-4
for bc in ['ns','sf']:
    c0 = corr_symbol(k, 0.2, 0.6, bc, 0); cp = corr_symbol(k, 0.2+h, 0.6, bc, 0); cm = corr_symbol(k, 0.2-h, 0.6, bc, 0)
    print(bc, 'd1 fd', (cp-cm)/(2*h), corr_symbol(k, 0.2, 0.6, bc, 1))
    print(bc, 'd2 fd', (cp-2*c0+cm)/h**2, corr_symbol(k, 0.2, 0.6, bc, 2))
