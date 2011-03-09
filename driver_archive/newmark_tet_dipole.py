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
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode import Mesh
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.DiscretisedSystem import CurlCurlNewmark

order = 1
eMAGUSImport.init('workspace')
mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())


print 'Mesh elements: ', len(mesh.elements)

dt = 1/30.
runtime = 6.
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
#newmarkSystem = CurlCurlNewmark(mesh, order=order, BC=free, useQ=True)
newmarkSystem = CurlCurlNewmark(mesh, order=order, BC=cb, useQ=False)

totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
weights, elPerm = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,0.]))
drive_dofnos = elPerm[1]
print "Done"
print "Getting test points"
direction = N.array([0, 1., 0])
#test_pts = PostProc.MakeLine(1/10.*direction, 9*direction, 121)
test_pts = N.array([[0, 0.5875, 0]])
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
newmarkSystem.step(int(N.ceil(runtime/dt)))
import pickle
# pickle.dump({'loggedReconstructed': newmarkSystem.loggedReconstructed,
#              'weights': weights, 'elPerm': elPerm, 'drive_dofnos': drive_dofnos},
#             file('+newmark_tet_dipole_tmp.pickle', 'w'))

from pylab import plot, xlabel, ylabel, legend, show

rs = newmarkSystem.loggedReconstructed[0]
rsvals = rs['vals']
ts = N.array([ts[0] for ts in rsvals])
plot(N.arange(len(ts))*dt, ts[:,0], label='E_x')
plot(N.arange(len(ts))*dt, ts[:,1], label='E_y')
plot(N.arange(len(ts))*dt, ts[:,2], label='E_z')
xlabel('time (s)')
ylabel('E field magnitude at r = [%f, %f, %f]' %
       (test_pts[0][0], test_pts[0][1], test_pts[0][2]) )
legend()