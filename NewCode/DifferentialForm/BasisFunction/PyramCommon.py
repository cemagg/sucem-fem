from __future__ import division

import numpy as N
import scipy

from NewCode.Utilities import gen_lagrange_polys, partial, Struct
from NewCode.Integration import gauss_lobatto_coeffs


#
# See graglia99 for definitions. cov_coord_vecs correspond to grad(xi_i)
# and con_coord_vecs correspond to grad(xi_{i+1}) x grad(xi_{i+2})
#

l1,l2,l3 = v = N.array([[1,0,0], [0,1,0], [0,0,1]], N.float64)
cov_coord_vecs = N.array([l1, l2, -l1-l3, -l2-l3, l3], N.float64)
con_coord_vecs = {(0,1): l3,    (1,2): l3-l1,    (2,3): l3-l1-l2,    (3,0): l3-l2,
                  (1,4): l1,    (2,4): l2,       (3,4): -l1,         (0,4): -l2,
                  (1,0): -(l3), (2,1): -(l3-l1), (3,2): -(l3-l1-l2), (0,3): -(l3-l2),
                  (4,1): -(l1), (4,2): -(l2),    (4,3): l1,          (4,0): l2,
}



shapefuns = N.array([lambda c: c[2]*c[3]/(1-c[4]),
                     lambda c: c[1]*c[2]/(1-c[4]),
                     lambda c: c[0]*c[3]/(1-c[4]),
                     lambda c: c[0]*c[1]/(1-c[4]),
                     lambda c: c[4]])

cvcv = cov_coord_vecs
D_shapefuns = N.array([
    lambda c: (c[2]*c[3]*cvcv[4]-c[2]*cvcv[3]*(c[4]-1)-cvcv[2]*c[3]*(c[4]-1))/(c[4]-1)**2,
    lambda c: (c[1]*c[2]*cvcv[4]-c[1]*cvcv[2]*(c[4]-1)-cvcv[1]*c[2]*(c[4]-1))/(c[4]-1)**2,
    lambda c: (c[0]*c[3]*cvcv[4]-c[0]*cvcv[3]*(c[4]-1)-cvcv[0]*c[3]*(c[4]-1))/(c[4]-1)**2,
    lambda c: (c[0]*c[1]*cvcv[4]-c[0]*cvcv[1]*(c[4]-1)-cvcv[0]*c[1]*(c[4]-1))/(c[4]-1)**2,
    lambda c: cvcv[4]])
