from __future__ import division

import pickle
import numpy as N
import dolfin

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
from ntff import NTFF

fname = 'dofs_sph-2-sphere-r1m-6.pickle'
#fname = 'dofs_sph-1lam-3-sphere-r1m-6.pickle'
#fname = 'dofs_sph-2lam-2-sphere-r1m-15.pickle'
theta_deg = N.linspace(0, 180, 91)
phi_deg = 0
data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh_file = '../solvers/meshes/%s.femmesh' % data['mesh_id']
mesh = femmesh_2_dolfin_mesh(mesh_file)
mesh.coordinates()[:] *= lam*1
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

