from __future__ import division

import dolfin
from dolfin import cross, ds, dot
import numpy as N

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshGenerators import get_centred_cube

order = 2
lam = 1
domain_size = N.array([1*lam]*3)
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len)
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", order)
n1 = V.cell().n
n2 = dolfin.FacetNormal(mesh)
u_expr = dolfin.Expression(['1', '0', '0'])
x_hat = dolfin.Expression(['1', '0', '0'])
y_hat = dolfin.Expression(['0', '1', '0'])
z_hat = dolfin.Expression(['0', '0', '1'])
u = dolfin.interpolate(u_expr, V)
surface_domains = dolfin.MeshFunction('uint', mesh, 2)
surface_domains.set_all(0)
class Boundary1(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and dolfin.near(x[1], -lam/2)
    
boundary1 = Boundary1()
boundary1.mark(surface_domains, 1)

integrand_n1 = dot(cross(n1, u), z_hat)*ds(1)
integrand_n2 = dot(cross(n2, u), z_hat)*ds(1)

intg_n1 = dolfin.assemble(integrand_n1, exterior_facet_domains=surface_domains)
intg_n2 = dolfin.assemble(integrand_n2, exterior_facet_domains=surface_domains)

