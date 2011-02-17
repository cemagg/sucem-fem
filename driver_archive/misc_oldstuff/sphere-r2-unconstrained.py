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
newmarkSystem = CurlCurlNewmark(mesh, order=1, BC=free, useQ=False)
totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
weights, elPerm = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,0.]))
drive_dofnos = elPerm[1]
print "Done"
print "Getting test points"
direction = N.array([0, 1., 0])
test_pts = PostProc.MakeLine(1/10.*direction, 1.9*direction, 241)
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
            file('+sphere-r2-unconstrained-runstuff.pickle', 'w'))


"""
Plotty stuff:

import analyse_fftstuff
E_dofs = N.array(newmarkSystem.loggedDOFs.values()[0]['vals']).transpose()
analyse_this = analyse_fftstuff.AnalyseThingy(E_dofs, dt, drv_fun)
analyse_this.plotFreqbars(an_freqs)
analyse_this.plotRange(20*log10(N.sum(N.abs(analyse_this.dfts_series), axis=0)/N.abs(analyse_this.fs)))

"""


"""
Results:

Using CT/LN newmark beta=0.3, dt = 0.02
                freq

mesh edge len   0.83333333,1.11803399,1.20185043x2,1.30170828x2,1.41421356,1.42400062,1.56347192x2

0.5/6             1/546       3/733      4/788          4/853      5/927       11/933     6/1024
0.5/12            

"""
