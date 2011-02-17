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
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmark

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 2

print 'Mesh elements: ', len(mesh.elements)

dt = 1/30.
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
newmarkSystem = CurlCurlNewmark(mesh, order=1, BC=cb)
totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
d=1/40.
weights1, elPerm1 = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,d/2.]))
weights2, elPerm2 = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,-d/2.]))
dof_weights = dict(zip(elPerm1[1], weights1))
for dof, weight in zip(elPerm2[1], weights2):
    if dof in dof_weights: dof_weights[dof] += weight
    else: dof_weights[dof] = weight

drive_dofnos = N.array(dof_weights.keys(), N.int)
weights = N.array(dof_weights.values(), N.float64)
print "Done"

print "Getting test points"
direction = N.array([0, 1., 0])
test_pts = PostProc.MakeLine(1/10.*direction, 0.9*direction, 121)
test_elnos, test_el_coords = PostProc.LocatePoints(
    newmarkSystem.disc.mesh, test_pts)
print "Done"
newmarkSystem.addReconstructedLogger(test_elnos, test_el_coords)

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform.D             # RHS proportional to d/dt of current
newmarkSystem.setDriveDOFs(drive_dofnos, weights, drv_fun)
newmarkSystem.useLU = True
#newmarkSystem.useLU = False
print "Setting timestep"
newmarkSystem.setTimestep(dt)
print "Done"
newmarkSystem.step(180)
import pickle
pickle.dump({'loggedReconstructed': newmarkSystem.loggedReconstructed,
             'weights': weights, 'drive_dofnos': drive_dofnos, 'drive_dofnos': drive_dofnos},
            file('+sphere-r1-two-pt-drive-intg-runstuff.pickle', 'w'))



