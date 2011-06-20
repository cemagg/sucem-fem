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
fname = 'reference_dofs-2-0.149896229-0.0499654096667.pickle'
#fname = 'dofs-2-0.599584916-0.0499654096667.pickle'
#fname = 'dofs-2-1.199169832-0.0499654096667.pickle'
#fname = 'interpdofs-2-0.599584916-0.0499654096667.pickle'
theta_deg = N.array([0, 13, 79, 90, 140, 179])
no_ff_pts = len(theta_deg)
phi_deg = N.linspace(0, 360, no_ff_pts)
#phi_deg = N.zeros(91)
data = pickle.load(open(fname))
lam = c0/data['freq']
k0 = data['freq']*2*N.pi/c0
mesh = get_centred_cube(data['domain_size'], data['max_edge_len'])
V = dolfin.FunctionSpace(mesh, "Nedelec 1st kind H(curl)", data['order'])

ntff = NTFF(V)
ntff.set_dofs(data['x'])
ntff.set_frequency(data['freq'])
E_ff = N.array([ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
E_theta = E_ff[:,0]
E_phi = E_ff[:,1]
E_theta_an = Z0*N.exp(-1j*k0)*(1j*k0*1*lam/1000)/(4*N.pi)*N.sin(N.deg2rad(theta_deg))
ff_result_data = dict(E_ff=E_ff, E_theta_analytical=E_theta_an,
                      theta_deg=theta_deg, phi_deg=phi_deg)
pickle.dump(dict(ff_result_data=ff_result_data, nf_input_data=data),
            open(fname.replace('dofs', 'surface_ntff'), 'w'))
import pylab 
pylab.figure()
pylab.plot(theta_deg, N.abs(E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(E_phi), label='|E_phi|')
pylab.plot(theta_deg, N.abs(E_theta_an), label='analytical')
pylab.legend()
pylab.grid(1)
pylab.show()

