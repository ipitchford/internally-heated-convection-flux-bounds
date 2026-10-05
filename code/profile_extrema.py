"""Finite-precision checks of a Chebyshev interpolant, not PDE certification."""
import numpy as np
from numpy.polynomial import Chebyshev, Polynomial


def inspect_profile(z, theta, fb, ft=None):
    z, theta = np.asarray(z), np.asarray(theta)
    if z.shape != theta.shape or not 0 < fb < 1:
        raise ValueError('Invalid profile or bottom flux')
    p = Chebyshev.fit(z, theta, len(z)-1, domain=[0, 1])
    slopes = (float(p.deriv()(0)), float(-p.deriv()(1)))
    cut = 1-fb
    points = []
    for lo, hi, line in [(0, cut, Polynomial([0, fb])),
                         (cut, 1, Polynomial([1-fb, -(1-fb)]))]:
        gap = p - line.convert(kind=Chebyshev, domain=[0, 1])
        roots = gap.deriv().roots()
        xs = [lo, hi] + [float(x.real) for x in roots
                         if abs(x.imag) < 1e-9 and lo < x.real < hi]
        points.extend((float(x), float(gap(x))) for x in xs)
    where, excess = max(points, key=lambda x: x[1])
    inner = (z > 1e-4) & (z < 1-1e-4)
    grid_ratio = float(np.max(theta[inner]/np.minimum(fb*z[inner], (1-fb)*(1-z[inner]))))
    return dict(grid_max_ratio=grid_ratio, grid_pass=grid_ratio <= 1+1e-9,
                bottom_slope=slopes[0], top_slope=slopes[1],
                top_limit_ratio=slopes[1]/(1-fb),
                measured_top_limit_ratio=None if ft is None else float(ft/(1-fb)),
                slope_compatible=slopes[0] <= fb+1e-9 and slopes[1] <= 1-fb+1e-9,
                interpolant_max_excess=excess, z_at_max_excess=where,
                interpolant_pass=excess <= 1e-11,
                assurance='floating-point spectral diagnostic; no continuum enclosure')
