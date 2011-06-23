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

fname = 'reference_dofs-2-0.149896229-0.0499654096667.pickle'
#fname = 'dofs-3-0.599584916-0.0499654096667.pickle'
#fname = 'dofs-2-1.199169832-0.0499654096667.pickle'
#fname = 'interpdofs-2-0.599584916-0.0499654096667.pickle'
no_ff_pts = 40
theta_deg = N.linspace(0, 180, no_ff_pts)
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


def run_pts(pts):
    E_H_ff = N.array([varntff.calc_pt_E_H(th_deg, ph_deg)
                      for th_deg, ph_deg in zip(theta_deg, phi_deg)[0:pts]])

import time
run_pts(1)
start = time.time()
run_pts(no_ff_pts)
stop = time.time()
print 'runtime for %d far-field points: %f (s)' % (no_ff_pts, stop-start)

