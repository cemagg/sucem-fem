from __future__ import division

import pickle
import numpy as N
import dolfin

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
from ntff import NTFF
import ntff_variational
reload(ntff_variational)
from ntff_variational import VariationalNTFF

fname = 'dofs_sph-2-sphere-r1m-6.pickle'
#fname = 'dofs_sph-1lam-3-sphere-r1m-6.pickle'
#fname = 'dofs_sph-2lam-2-sphere-r1m-15.pickle'
theta_deg = N.linspace(0, 180, 71)
phi_deg = 0
testing_boost = 1 # Order by with to boost the testing function representation
data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh_file = '../solvers/meshes/%s.femmesh' % data['mesh_id']
mesh = femmesh_2_dolfin_mesh(mesh_file)
mesh.coordinates()[:] *= lam*1
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])
Vt = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order']+testing_boost)

varntff = VariationalNTFF(V, Vt)
varntff.set_k0(k0)
varntff.set_E_dofs(data['x'])
E_H_ff = N.array([varntff.calc_pt(th_deg, phi_deg) for th_deg in theta_deg])

ntff = NTFF(V)
ntff.set_dofs(data['x'])
ntff.set_frequency(data['freq'])
E_ff = N.array([ntff.calc_pt(th_deg, phi_deg) for th_deg in theta_deg])
# E_theta = E_ff[:,0]
# E_phi = E_ff[:,1]
_L = N.array(ntff._L)
_N = N.array(ntff._N)
_r_fac = N.array(ntff._r_fac)
E_theta = -_r_fac*(_L[:,1] + Z0*_N[:,0])
E_phi = _r_fac*(_L[:,0] - Z0*_N[:,1])


E_theta_new = _r_fac*(-_L[:,1] + E_H_ff[:,0])
E_phi_new = _r_fac*(_L[:,0] + E_H_ff[:,1])
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
err_old = N.sqrt((N.abs(E_phi)**2 + N.abs(-E_theta - E_theta_an)**2))
err_new = N.sqrt((N.abs(E_phi_new)**2 + N.abs(-E_theta_new - E_theta_an)**2))
err_old_th = N.sqrt(N.abs(-E_theta - E_theta_an)**2)
err_new_th = N.sqrt(N.abs(-E_theta_new - E_theta_an)**2)

relerr_old = 100*err_old/N.abs(E_theta_an)
relerr_new = 100*err_new/N.abs(E_theta_an)
relerr_old_th = 100*err_old_th/N.abs(E_theta_an)
relerr_new_th = 100*err_new_th/N.abs(E_theta_an)

pattern_energy_an = N.sqrt(N.sum(N.abs(E_theta_an)**2))
pattern_energy_old_th = N.sqrt(N.sum(N.abs(E_theta)**2))
pattern_energy_new_th = N.sqrt(N.sum(N.abs(E_theta_new)**2))
pattern_energy_old_ph = N.sqrt(N.sum(N.abs(E_phi)**2))    
pattern_energy_new_ph = N.sqrt(N.sum(N.abs(E_phi_new)**2))

pattern_err_old_th = 100*N.abs(pattern_energy_old_th - pattern_energy_an)/pattern_energy_an
pattern_err_new_th = 100*N.abs(pattern_energy_new_th - pattern_energy_an)/pattern_energy_an
pattern_err_old_ph = 100*N.abs(pattern_energy_old_ph)/pattern_energy_an
pattern_err_new_ph = 100*N.abs(pattern_energy_new_ph)/pattern_energy_an
pattern_err_old = pattern_err_old_th + pattern_err_old_ph
pattern_err_new = pattern_err_new_th + pattern_err_new_ph


beam_energy_an = N.sqrt(N.sum(N.abs(E_theta_an[18:-18])**2))     
beam_energy_old_th = N.sqrt(N.sum(N.abs(E_theta[18:-18])**2))    
beam_energy_new_th = N.sqrt(N.sum(N.abs(E_theta_new[18:-18])**2))
beam_energy_old_ph = N.sqrt(N.sum(N.abs(E_phi[18:-18])**2))    
beam_energy_new_ph = N.sqrt(N.sum(N.abs(E_phi_new[18:-18])**2))

beam_err_old_th = 100*N.abs(beam_energy_old_th - beam_energy_an)/beam_energy_an
beam_err_new_th = 100*N.abs(beam_energy_new_th - beam_energy_an)/beam_energy_an
beam_err_old_ph = 100*N.abs(beam_energy_old_ph)/beam_energy_an
beam_err_new_ph = 100*N.abs(beam_energy_new_ph)/beam_energy_an
beam_err_old = beam_err_old_th + beam_err_old_ph
beam_err_new = beam_err_new_th + beam_err_new_ph

#| Run | pat_surf | pat_vol | beam_surf | beam_vol | beam_th_surf | beam_th_vol |
print '| | %f | %f | %f | %f | %f | %f |' % (pattern_err_old, pattern_err_new,
                                             beam_err_old, beam_err_new,
                                             beam_err_old_th, beam_err_new_th)

import pylab 
pylab.plot(theta_deg, N.abs(E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(E_theta_new), label='|E_theta_new|')
pylab.plot(theta_deg, N.abs(E_phi), label='|E_phi|')
pylab.plot(theta_deg, N.abs(E_phi_new), label='|E_phi_new|')

pylab.plot(theta_deg, N.abs(E_theta_an), label='analytical')
pylab.legend()
pylab.grid(1)
pylab.show()

