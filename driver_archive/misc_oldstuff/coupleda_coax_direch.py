"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
from itertools import izip
import os, sys
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
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, SystemMatrix, PostProc
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.DifferentialForm import DiscretiserDOFs
from NewCode.DiscretisedSystem import CoupledFirstOrderSystemHardSource
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets

import coax_test_setup

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 4

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
constrained = lambda x: False

coax_test_setup.setupMatcher(mesh)
wgm = coax_test_setup.wgm
onPort = coax_test_setup.onPort
dt = coax_test_setup.dt

fs = wgm.fs

freeE = wgm.superFree
freeB = lambda ent: freeE(ent) or wgm.onPort(ent)
freeH = free
freeD = free 
direchFree = lambda edge: edge.index in wgm.work_edgeno_set
# All lowest order
#orders = {'E':(1,True), 'B':(1,True)}
coupledSystem = CoupledFirstOrderSystemHardSource(
    mesh, dt=dt, BCs=Struct(E=freeE, B=freeB, H=freeH, D=freeD))
coupledSystem.setDirechBCs(Struct(
    E=direchFree, B=constrained, H=constrained, D=constrained))
T, omega = 0.5, 5.523599
# def drv_fun(dt, n):
#     t = dt*n
#     return (1. - N.exp(-(t/2/T)**2))*N.sin(omega*t)

drv_fun = coax_test_setup.drv_fun
coupledSystem.setSourceDOFs(wgm.dofs, drv_fun)


direction = coax_test_setup.direction
xypos = coax_test_setup.xypos
test_pts = coax_test_setup.test_pts
test_elnos, test_el_coords = coax_test_setup.test_elnos, coax_test_setup.test_el_coords
coupledSystem.addReconstructedLogger('E', test_elnos, test_el_coords)


DiscretiserDOFs.PformDiscretiserBase.useLU = True
coupledSystem.useLU = True
coupledSystem.step(1400)

rs = coupledSystem.loggedReconstructed.E[0]
rsvals = N.array(rs['vals'])

import analyse_fftstuff
class AnalyseThingyBoxcar(analyse_fftstuff.AnalyseThingy):
    windowFun = staticmethod(scipy.signal.boxcar)
at = AnalyseThingyBoxcar([rsvals[:,200,0]], dt, drv_fun)
at.xmin = 300
at.xmax = 150*8

TF = at.dfts_series[0]/at.fs

plot(N.arange(at.xmin, at.xmax)*at.df, N.unwrap(N.angle(TF[at.xmin:at.xmax]))/5)


