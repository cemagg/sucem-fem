"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import random
import numpy as N
import scipy
#
# Local Imports
#
import sys
sys.path.append('../../../')
sys.path.append('../')
import NewCode
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.DiscretisedSystem import CurlCurlNewmark
from NewCode.Meshes import CalculateConnectivity
from NewCode.Meshes import Conversions
from NewCode.Meshes.MeshIO import Femmesh

import test_data
from FenicsCode.Consts import eps0, mu0, c0, Z0

parameters = dict(test_data.problem_data)
freq = parameters['f'];
lam = c0/freq
k_0 = 2*N.pi*freq/c0
source_coord = N.array([0,0,0.]) + 1e-5
source_value = N.array([0,0,1.])*parameters['I']*parameters['l']

import dolfin as dol
source_point = dol.Point(*source_coord)
dol_mesh = dol.UnitTetrahedron()
dol_mesh.coordinates()[:] *= lam
listmesh = Conversions.dolfin_mesh_2_listmesh(dol_mesh)
listmesh = CalculateConnectivity.get_all_connectivities(listmesh)
## End dolfin mesh setup

mesh = Mesh.Mesh(listmesh)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 6
CurlCurlNewmark.triIntegrationOrder = 6

print 'Mesh elements: ', len(mesh.elements)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
order = 1
newmarkSystem = CurlCurlNewmark(mesh, order=order, BC=free, useQ=True)
totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
weights, elPerm = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: source_value, r0=source_coord)
drive_dofnos = elPerm[1]
drive_el, drive_coord = PostProc.LocatePoints(newmarkSystem.disc.mesh, [source_coord])
print "Getting test points"
test_pts = test_data.r_1[0:10]
test_elnos, test_el_coords = PostProc.LocatePoints(
    newmarkSystem.disc.mesh, test_pts)
print "Done"

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform.D             # RHS proportional to d/dt of current

M = newmarkSystem.massMatrix()
S = newmarkSystem.stiffnessMatrix()
S_0 = newmarkSystem.boundarySurfaceMatrix()

no_dofs = M.shape[0]
print 'no_dofs: ', no_dofs
b = N.zeros(no_dofs, dtype=N.complex128)
rhs_contrib = 1j*k_0*Z0*weights
b[drive_dofnos] = rhs_contrib
print "Building system matrix"
A = S - k_0**2*M + 1j*k_0*S_0

import scipy.sparse.linalg
print "Solving system matrix"
A_lu = scipy.sparse.linalg.factorized(A)
x = A_lu(b)

dofs = newmarkSystem.disc.newDOFs()
reconparms = list(dofs.calc_reconparms(test_elnos, test_el_coords))
E_1 = dofs.recon_fromparms(reconparms, x)
