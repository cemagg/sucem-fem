"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import random
import numpy as N
import scipy
from scipy import integrate

#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet

import coax_test_setup

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)


Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 4

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

coax_test_setup.setupMatcher(mesh)
wgm = coax_test_setup.wgm
onPort = coax_test_setup.onPort
dt = coax_test_setup.dt

direchFree = lambda edge: edge.index in wgm.work_edgeno_set
freeE = wgm.superFree

# All lowest order

newmarkSystem = CurlCurlNewmarkDirechlet(mesh, order=1, BC=freeE, useQ=False)
newmarkSystem.setDirechBCs(direchFree)
newmarkSystem.direchSys.dofs.dofArray[:] = wgm.dofs

drv_fun = coax_test_setup.drv_fun
newmarkSystem.drive_fun = drv_fun
newmarkSystem.useLU = True
print "Setting timestep"
newmarkSystem.setTimestep(dt)
print "Done"


direction = coax_test_setup.direction
xypos = coax_test_setup.xypos
test_pts = coax_test_setup.test_pts
test_elnos, test_el_coords = coax_test_setup.test_elnos, coax_test_setup.test_el_coords
newmarkSystem.addReconstructedLogger(test_elnos, test_el_coords)

newmarkSystem.step(1400)

rs = newmarkSystem.loggedReconstructed[0]
rsvals = N.array(rs['vals'])
# plot(test_pts[:,2],rsvals[1000][:,1])

import analyse_fftstuff
class AnalyseThingyBoxcar(analyse_fftstuff.AnalyseThingy):
    windowFun = staticmethod(scipy.signal.boxcar)
    
at = AnalyseThingyBoxcar([rsvals[:,200,0]], dt, drv_fun)
at.xmin = 10
at.xmax = 150*4

TF = at.dfts_series[0]/at.fs

plot(N.arange(at.xmin, at.xmax)*at.df, N.unwrap(N.angle(TF[at.xmin:at.xmax])))



