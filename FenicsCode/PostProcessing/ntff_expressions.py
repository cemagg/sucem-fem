# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin
from dolfin import dot

get_r_hat = lambda : dolfin.Expression(
    ['sin(theta)*cos(phi)', 'sin(theta)*sin(phi)', 'cos(theta)'],
    theta=0, phi=0)

get_k0 = lambda : dolfin.Expression('k0', k0=0)

get_theta_hat = lambda : dolfin.Expression(
    ['cos(theta)*cos(phi)', 'cos(theta)*sin(phi)', '-sin(theta)'],
    theta=0, phi=0)

get_phi_hat = lambda : dolfin.Expression(['-sin(phi)', 'cos(phi)', '0.'], phi=0)

get_phase = lambda k0, rprime, r_hat : k0*dot(rprime, r_hat)

get_3d_vector = lambda : dolfin.Expression(['x1', 'x2', 'x3'],
                                           x1=0, x2=0, x3=0)
