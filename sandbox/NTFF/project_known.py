from __future__ import division

import dolfin
import math
import numpy
import numpy as N
import pickle

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshGenerators import get_centred_cube

class E_ff_Re(dolfin.Expression):
    part = 'real'
    freq = 1e9
    l = 1.
    I = 1.
    def eval(self, value, x):
        k0 = 2*math.pi*self.freq/c0
        r = max(math.sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2)), 1e-20)
        theta = math.acos(x[2]/r)
        phi = math.atan2(x[1], x[0])
        theta_hat_x = math.cos(theta)*math.cos(phi)
        theta_hat_y = math.cos(theta)*math.sin(phi)
        theta_hat_z = -math.sin(theta)
        H_phi = 1j*k0*self.l*self.I/(math.pi*4)*math.sin(theta)
        H_phi = H_phi*numpy.exp(-1j*k0*r)/r
        E_theta_complex = Z0*H_phi
        E_theta = getattr(E_theta_complex, self.part)
        value[0] = theta_hat_x*E_theta
        value[1] = theta_hat_y*E_theta
        value[2] = theta_hat_z*E_theta

    def value_shape(self):
        return (3,)

class E_ff_Im(dolfin.Expression):
    part = 'imag'
    freq = 1e9
    l = 1.
    I = 1.
    def eval(self, value, x):
        k0 = 2*math.pi*self.freq/c0
        r = max(math.sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2)), 1e-20)
        theta = math.acos(x[2]/r)
        phi = math.atan2(x[1], x[0])
        theta_hat_x = math.cos(theta)*math.cos(phi)
        theta_hat_y = math.cos(theta)*math.sin(phi)
        theta_hat_z = -math.sin(theta)
        H_phi = 1j*k0*self.l*self.I/(math.pi*4)*math.sin(theta)
        H_phi = H_phi*numpy.exp(-1j*k0*r)/r
        E_theta_complex = Z0*H_phi
        E_theta = getattr(E_theta_complex, self.part)
        value[0] = theta_hat_x*E_theta
        value[1] = theta_hat_y*E_theta
        value[2] = theta_hat_z*E_theta

    def value_shape(self):
        return (3,)



## Problem parameters
freq = 1e9
lam = c0/freq
I = 1.                                  # Dipole current
l = lam/1000                            # Dipole length
source_value = N.array([0,0,1.])*I*l
source_coord = N.array([0,0,0.]) 
## Discretisation settings
order = 3
domain_size = N.array([2*lam]*3)
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len, source_coord)
print 'No mesh elements: %d' % mesh.num_cells()
E_re = E_ff_Re()
E_re.I = I
E_re.freq = freq
E_re.l = l
E_im = E_ff_Re()
E_im.I = I
E_im.freq = freq
E_im.l = l

data = dict(order=order, freq=freq, domain_size=domain_size,
            max_edge_len=max_edge_len)

V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])
V2 = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(div)", data['order'])
x_re = dolfin.interpolate(E_re, V)
x_im = dolfin.interpolate(E_im, V)
x = N.array(x_re.vector()[:]) + 1j*N.array(x_im.vector()[:])
data['x'] = x

fname = 'interpdofs-%d-%s-%s.pickle' % (order, str(domain_size[0]), str(max_edge_len))
pickle.dump(data, open(fname, 'w'))

mesh.intersection_operator().clear()
fac = 1
mesh.coordinates()[:] *= fac
plot_pts = N.linspace(0, 0.95*lam, 101)[:, N.newaxis]*N.array([1,1,1])*fac
E_re_expr = N.zeros((len(plot_pts), 3))
E_im_expr = N.zeros((len(plot_pts), 3))
E_re_num = N.zeros((len(plot_pts), 3))
E_im_num = N.zeros((len(plot_pts), 3))
curl_E_re_num = N.zeros((len(plot_pts), 3))
curl_E_im_num = N.zeros((len(plot_pts), 3))
curl_x_re = dolfin.project(dolfin.curl(x_re), V2)
curl_x_im = dolfin.project(dolfin.curl(x_im), V2)

for i, pt in enumerate(plot_pts):
    E_re.eval(E_re_expr[i,:], pt)
    E_im.eval(E_im_expr[i,:], pt)
    E_re_num[i,:] = x_re(pt)
    E_im_num[i,:] = x_im(pt)
    curl_E_re_num[i,:] = curl_x_re(pt)
    curl_E_im_num[i,:] = curl_x_im(pt)


import pylab
r_pts = N.sqrt(N.sum(plot_pts**2, axis=1))/lam/fac
pylab.plot(r_pts, E_re_expr[:,0], label='re expr x')
pylab.plot(r_pts, E_re_expr[:,1], label='re expr y')
pylab.plot(r_pts, E_re_expr[:,2], label='re expr z')
pylab.plot(r_pts, E_re_num[:,0], label='re num x')
pylab.plot(r_pts, E_re_num[:,1], label='re num y')
pylab.plot(r_pts, E_re_num[:,2], label='re num z')
# pylab.plot(r_pts, curl_E_re_num[:,0], label='curl re num x')
# pylab.plot(r_pts, curl_E_re_num[:,1], label='curl re num y')
# pylab.plot(r_pts, curl_E_re_num[:,2], label='curl re num z')


pylab.legend(loc='best')
pylab.grid(1)
pylab.show()
