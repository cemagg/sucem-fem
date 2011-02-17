"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import numpy as N
import scipy
import random
from scipy import signal
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
from NewCode.DiscretisedSystem import CoupledFirstOrderSystem

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 4

dt = 1/150.
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
orders = {'E':(2,True), 'H':(2,True), 'D':(1,False), 'B':(1,False)}
coupledSystem = CoupledFirstOrderSystem(
    mesh, dt=dt, BCs=Struct(E=cb, H=free, D=free, B=cb),
    disc_orders=orders)

for discName, disc in coupledSystem.discs.items():
    print discName, disc.basisSet.info, disc.totalDOFs

totalDOFs_D = coupledSystem.discs.D.totalDOFs
workDOFs_D = min(int(N.ceil(totalDOFs_D/100.)), 1000)
totalDOFs_E = coupledSystem.discs.E.totalDOFs
workDOFs_E = min(int(N.ceil(totalDOFs_E/100.)), 1000)

print "Getting source DOFs"
# weights, elPerm = coupledSystem.dofs.D.calcProjPointfunRHS_with_elPerm(
#     matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,0.]))
# drive_dofnos = elPerm[1]
coupledSystem.dofs.D.matchPointFunction(matchfun=lambda r: N.array([0,0,1.]),
                                        r0=N.array([0,0,0.]))
drive_dofnos = N.argwhere(N.abs(coupledSystem.dofs.D.dofArray) > 0).flatten()
weights = coupledSystem.dofs.D.dofArray[drive_dofnos]
coupledSystem.dofs.D.zero()
print "Done"
print "Getting test points"
direction = N.array([0, 1., 0])
test_pts = PostProc.MakeLine(1/10.*direction, .9*direction, 241)
test_elnos, test_el_coords = PostProc.LocatePoints(
    coupledSystem.mesh, test_pts)
print "Done"
coupledSystem.addReconstructedLogger('E', test_elnos, test_el_coords)

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform             # RHS proportional to the current

coupledSystem.setDriveDOFs_J(drive_dofnos, weights, drv_fun)
dofs = coupledSystem.dofs
from NewCode.MatrixUtils import add_diagonal_preconditioner
add_diagonal_preconditioner(dofs.E.matrix.mass())
add_diagonal_preconditioner(dofs.H.matrix.mass())
coupledSystem.step(100)

from pylab import *
rs = coupledSystem.loggedReconstructed.E[0]
rsvals = rs['vals']
ts120 = N.array([ts[120] for ts in rsvals])
plot(N.arange(len(ts120))*dt, ts120[:,0], label='E_x')
plot(N.arange(len(ts120))*dt, ts120[:,1], label='E_y')
plot(N.arange(len(ts120))*dt, ts120[:,2], label='E_z')
xlabel('time (s)')
ylabel('E field magnitude at r = [0, 0.5, 0]')
legend()
