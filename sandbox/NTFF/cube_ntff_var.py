from __future__ import division

import pickle
import numpy as N
import dolfin

import sys
sys.path.append('../../')
from FenicsCode.Consts import Z0, c0
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from surface_ntff import NTFF
import variational_ntff
reload(variational_ntff)
from variational_ntff import NTFF as var_NTFF

dolfin.parameters['optimize_form'] = True
dolfin.parameters['optimize'] = True
dolfin.parameters['optimize_use_dofmap_cache'] = True
dolfin.parameters['optimize_use_tensor_cache'] = True
dolfin.parameters['form_compiler']['optimize'] = True
dolfin.parameters['form_compiler']['cpp_optimize'] = True


#fname = 'reference_dofs-2-0.149896229-0.0499654096667.pickle'
fname = 'dofs-3-0.149896229-0.0499654096667.pickle'
#fname = 'dofs-3-0.599584916-0.0499654096667.pickle'
#fname = 'dofs-2-1.199169832-0.0499654096667.pickle'
#fname = 'interpdofs-2-0.599584916-0.0499654096667.pickle'
theta_deg = N.linspace(0, 180, 71)
#theta_deg = N.array([0, 13, 79, 90, 140, 179])
no_ff_pts = len(theta_deg)
phi_deg = N.zeros(no_ff_pts)
#phi_deg = N.linspace(0, 360, no_ff_pts)
testing_boost = 0 # Order by with to boost the testing function representation
data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh = get_centred_cube(data['domain_size'], data['max_edge_len'])
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])
Vt = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order']+testing_boost)

varntff = var_NTFF(V, Vt)
varntff.set_k0(k0)
varntff.set_dofs(data['x'])
E_H_ff = N.array([varntff.calc_pt_E_H(th_deg, ph_deg)
                  for th_deg, ph_deg in zip(theta_deg, phi_deg)])

ntff = NTFF(V)
ntff.set_dofs(data['x'])
ntff.set_frequency(data['freq'])
E_ff_old = N.array([ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
_L = N.array(ntff._L)
_N = N.array(ntff._N)
_r_fac = N.array(ntff._r_fac)
E_theta_old = -_r_fac*(_L[:,1] + Z0*_N[:,0])
E_phi_old = _r_fac*(_L[:,0] - Z0*_N[:,1])


E_theta = _r_fac*(-_L[:,1] + E_H_ff[:,0])
E_phi = _r_fac*(_L[:,0] + E_H_ff[:,1])
E_ff = N.vstack((E_theta, E_phi)).T
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
err_old = N.sqrt((N.abs(E_phi_old)**2 + N.abs(-E_theta_old - E_theta_an)**2))
err_new = N.sqrt((N.abs(E_phi)**2 + N.abs(-E_theta - E_theta_an)**2))
err_old_th = N.sqrt(N.abs(-E_theta_old - E_theta_an)**2)
err_new_th = N.sqrt(N.abs(-E_theta - E_theta_an)**2)

relerr_old = 100*err_old/N.abs(E_theta_an)
relerr_new = 100*err_new/N.abs(E_theta_an)
relerr_old_th = 100*err_old_th/N.abs(E_theta_an)
relerr_new_th = 100*err_new_th/N.abs(E_theta_an)

pattern_energy_an = N.sqrt(N.sum(N.abs(E_theta_an)**2))
pattern_energy_old_th = N.sqrt(N.sum(N.abs(E_theta_old)**2))
pattern_energy_new_th = N.sqrt(N.sum(N.abs(E_theta)**2))
pattern_energy_old_ph = N.sqrt(N.sum(N.abs(E_phi_old)**2))    
pattern_energy_new_ph = N.sqrt(N.sum(N.abs(E_phi)**2))

pattern_err_old_th = 100*N.abs(pattern_energy_old_th - pattern_energy_an)/pattern_energy_an
pattern_err_new_th = 100*N.abs(pattern_energy_new_th - pattern_energy_an)/pattern_energy_an
pattern_err_old_ph = 100*N.abs(pattern_energy_old_ph)/pattern_energy_an
pattern_err_new_ph = 100*N.abs(pattern_energy_new_ph)/pattern_energy_an
pattern_err_old = pattern_err_old_th + pattern_err_old_ph
pattern_err_new = pattern_err_new_th + pattern_err_new_ph


beam_energy_an = N.sqrt(N.sum(N.abs(E_theta_an[18:-18])**2))     
beam_energy_old_th = N.sqrt(N.sum(N.abs(E_theta_old[18:-18])**2))    
beam_energy_new_th = N.sqrt(N.sum(N.abs(E_theta[18:-18])**2))
beam_energy_old_ph = N.sqrt(N.sum(N.abs(E_phi_old[18:-18])**2))    
beam_energy_new_ph = N.sqrt(N.sum(N.abs(E_phi[18:-18])**2))

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
ff_result_data = dict(E_ff=E_ff, E_theta_analytical=E_theta_an,
                      theta_deg=theta_deg, phi_deg=phi_deg)
# pickle.dump(dict(ff_result_data=ff_result_data, nf_input_data=data),
#             open(fname.replace('dofs', 'variational_ntff'), 'w'))



import pylab 
pylab.plot(theta_deg, N.abs(E_theta_old), label='|E_theta_old|')
pylab.plot(theta_deg, N.abs(E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(E_phi_old), label='|E_phi_old|')
pylab.plot(theta_deg, N.abs(E_phi), label='|E_phi|')

pylab.plot(theta_deg, N.abs(E_theta_an), label='analytical')
pylab.legend()
pylab.grid(1)
pylab.show()

