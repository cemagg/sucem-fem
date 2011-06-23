import dolfin
from dolfin import dot

get_r_hat = lambda : dolfin.Expression(
    ['sin(theta)*cos(phi)', 'sin(theta)*sin(phi)', 'cos(theta)'])

get_k0 = lambda : dolfin.Expression('k0')

get_theta_hat = lambda : dolfin.Expression(
    ['cos(theta)*cos(phi)', 'cos(theta)*sin(phi)', '-sin(theta)'])

get_phi_hat = lambda : dolfin.Expression(['-sin(phi)', 'cos(phi)', '0.'])

get_phase = lambda k0, rprime, r_hat : k0*dot(rprime, r_hat)

get_3d_vector = lambda : dolfin.Expression(['x1', 'x2', 'x3'])
