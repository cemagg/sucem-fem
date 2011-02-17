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
from NewCode.Meshes import BrickMesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, SystemMatrix, PostProc
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DifferentialForm import DiscretiserDOFs
from NewCode.DiscretisedSystem import BrickCoupledFirstOrderSystemBDirechlet

from wg_dof_matching_stuff import RoughWGMatcher
import brick_cavity

class BrickRoughWGMatcher(RoughWGMatcher):
    DiscretiserModule = BrickDiscretiser

a,b,l = [getattr(BrickRoughWGMatcher, x) for x in ("a", "b", "l")]

mesh = BrickMesh.Mesh(brick_cavity.make_rect_cavity_brick_listmesh(
    a,b,l, 1/10.))

print 'Mesh elements: ', len(mesh.elements)

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

wgm = BrickRoughWGMatcher(mesh)

fs = wgm.fs
onPort = wgm.onPort

freeE = wgm.superFree
freeB = lambda ent: freeE(ent) or wgm.onPort(ent)
direchFree = lambda edge: edge.index in wgm.work_edgeno_set
dt = 0.01
# All lowest order
orders = {'E':(1,True), 'B':(1,True)}
coupledSystem = BrickCoupledFirstOrderSystemBDirechlet(
    mesh, dt=dt, BCs=Struct(E=freeE, B=freeB), disc_orders=orders)
coupledSystem.setDirechBCs(Struct(E=direchFree, B=lambda x: False))
coupledSystem.direchSys.dofs.E.dofArray[:] = wgm.dofs

T, omega = 0.5, 5.523599
def drv_fun(dt, n):
    t = dt*n
    return (1. - N.exp(-(t/2/T)**2))*N.sin(omega*t)

coupledSystem.drv_fun = drv_fun


direction = N.array([0, 0, 1.], N.float64)
xypos = [wgm.a/2, wgm.b/2, 0]
test_pts = PostProc.MakeLine(0.25*direction, 9.75*direction, 400) + xypos
test_elnos, test_el_coords = PostProc.LocatePoints(coupledSystem.mesh, test_pts)
coupledSystem.addReconstructedLogger('E', test_elnos, test_el_coords)


DiscretiserDOFs.PformDiscretiserBase.useLU = True
coupledSystem.useLU = True
coupledSystem.step(2048)

rs = coupledSystem.loggedReconstructed.E[0]
rsvals = rs['vals']
plot(test_pts[:,2],rsvals[1000][:,1])



