from __future__ import division
import numpy as N
from pylab import plot
from NewCode.Integration import gauss_lobatto_coeffs
from NewCode.Utilities import gen_lagrange_polys
porder = 3
polys = gen_lagrange_polys(gauss_lobatto_coeffs(porder)[0]/2+0.5)
x = N.linspace(0,1,200)

for p in polys:
    plot(x, [p(xx)**2 for xx in x])
