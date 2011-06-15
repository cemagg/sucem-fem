from __future__ import division

import pickle
import numpy as N
import dolfin

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from ntff import NTFF


#dolfin.parameters['form_compiler']['quadrature_degree'] = '8'
#dolfin.parameters['form_compiler']['quadrature_degree'] = '8'
fname = 'dofs-2-0.599584916-0.0499654096667.pickle'
#fname = 'dofs-2-1.199169832-0.0499654096667.pickle'
#fname = 'interpdofs-2-0.599584916-0.0499654096667.pickle'
theta_deg = N.linspace(0, 180, 91)
phi_deg = 0
data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh = get_centred_cube(data['domain_size'], data['max_edge_len'])
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])

ntff = NTFF(V)
ntff.set_dofs(data['x'])
ntff.set_frequency(data['freq'])
E_ff = N.array([ntff.calc_pt(th_deg, phi_deg) for th_deg in theta_deg])
E_theta = E_ff[:,0]
E_phi = E_ff[:,1]

import pylab 
pylab.plot(theta_deg, N.abs(E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(E_phi), label='|E_phi|')
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
pylab.plot(theta_deg, N.abs(E_theta_an), label='analytical')
pylab.legend()
pylab.grid(1)
pylab.show()

