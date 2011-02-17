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
newmarkSystem = CurlCurlNewmark(mesh, order=1, BC=cb, useQ=False)
totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
weights, elPerm = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,0.]))
drive_dofnos = elPerm[1]
print "Done"
print "Getting test points"
direction = N.array([0, 1., 0])
#test_pts = PostProc.MakeLine(1/10.*direction, 9*direction, 121)
test_pts = N.array([0, 0.5875, 0])
test_elnos, test_el_coords = PostProc.LocatePoints(
    newmarkSystem.disc.mesh, test_pts)
print "Done"
newmarkSystem.addReconstructedLogger(test_elnos, test_el_coords)

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform.D             # RHS proportional to d/dt of current
newmarkSystem.setDriveDOFs(drive_dofnos, weights, drv_fun)
#newmarkSystem.useLU = True
newmarkSystem.useLU = False
print "Setting timestep"
newmarkSystem.setTimestep(dt)
print "Done"
newmarkSystem.step(180)
import pickle
pickle.dump({'loggedReconstructed': newmarkSystem.loggedReconstructed,
             'weights': weights, 'elPerm': elPerm, 'drive_dofnos': drive_dofnos},
            file('+sphere-r2-PEC-runstuff.pickle', 'w'))

rs = newmarkSystem.loggedReconstructed[0]
rsvals = rs['vals']
ts65 = N.array([ts[65] for ts in rsvals])
plot(N.arange(len(ts65))*dt, ts65[:,0], label='E_x')
plot(N.arange(len(ts65))*dt, ts65[:,1], label='E_y')
plot(N.arange(len(ts65))*dt, ts65[:,2], label='E_z')
legend()

