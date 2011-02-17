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
from NewCode.Meshes import  BrickMesh
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet

from wg_dof_matching_stuff import RoughWGMatcher
import brick_cavity

class BrickRoughWGMatcher(RoughWGMatcher):
    DiscretiserModule = BrickDiscretiser

a,b,l = [getattr(BrickRoughWGMatcher, x) for x in ("a", "b", "l")]

mesh = BrickMesh.Mesh(brick_cavity.make_rect_cavity_brick_listmesh(
    a,b,l, 1/10.))

print 'Mesh elements: ', len(mesh.elements)

class BrickCurlCurlNewmarkDirechlet(CurlCurlNewmarkDirechlet):
    DiscretiserModule = BrickDiscretiser

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

wgm = BrickRoughWGMatcher(mesh)
fs = wgm.fs
onPort = wgm.onPort

direchFree = lambda edge: edge.index in wgm.work_edgeno_set
freeE = wgm.superFree

dt = 0.01
# All lowest order

newmarkSystem = BrickCurlCurlNewmarkDirechlet(mesh, order=1, BC=freeE, useQ=False)
newmarkSystem.setDirechBCs(direchFree)
newmarkSystem.direchSys.dofs.dofArray[:] = wgm.dofs

T, omega = 0.5, 5.523599
def drv_fun(dt, n):
    t = dt*n
    return (1. - N.exp(-(t/2/T)**2))*N.sin(omega*t)

newmarkSystem.drive_fun = drv_fun
newmarkSystem.useLU = True
print "Setting timestep"
newmarkSystem.setTimestep(dt)
print "Done"


direction = N.array([0, 0, 1.], N.float64)
xypos = [wgm.a/2, wgm.b/2, 0]
test_pts = PostProc.MakeLine(0.25*direction, 9.75*direction, 400) + xypos
test_elnos, test_el_coords = PostProc.LocatePoints(newmarkSystem.disc.mesh, test_pts)
newmarkSystem.addReconstructedLogger(test_elnos, test_el_coords)

cd = N.array([1,0,0], N.float64)
test_crossection = PostProc.MakeLine(0.1*cd, cd*wgm.a*0.9, 51) + [0,0,2.2]
newmarkSystem.addReconstructedLogger(
    *PostProc.LocatePoints(newmarkSystem.disc.mesh, test_crossection))

newmarkSystem.step(2048)

rs = newmarkSystem.loggedReconstructed[0]
rsc = newmarkSystem.loggedReconstructed[1]
rsvals = rs['vals']
rsvalsc = rsc['vals']
plot(test_pts[:,2],rsvals[1000][:,1])

