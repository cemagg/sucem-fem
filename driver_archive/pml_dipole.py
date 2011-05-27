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
sys.path.append('../')
import NewCode
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.DiscretisedSystem import PMLSystem

h = 1/7.5
geom_size = 2.
no_PML_cells = 5
PML_m = 3
sigma_factor = 1
runtime = 6.

a,b,c = [geom_size + h*no_PML_cells*2]*3
order = 1

mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(
    a,b,c, [h, h, h], grid_offset=[-a/2,-b/2,-c/2]))

real_h = mesh.gridStepSize[0]

PML_pos = a/2 - real_h*no_PML_cells
PML_len = real_h*no_PML_cells
wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/real_h/wave_imp*sigma_factor

sigma_u = lambda u: N.where(N.abs(u) >= PML_pos,
                            PML_sigma_max*((N.abs(u)-PML_pos)/PML_len)**PML_m,
                            0.)
    
sigma_x_fn = lambda r: sigma_u(r.T[0])
sigma_y_fn = lambda r: sigma_u(r.T[1])
sigma_z_fn = lambda r: sigma_u(r.T[2])
sigma_fns = Struct(x=sigma_x_fn, y=sigma_y_fn, z=sigma_z_fn)

print 'Mesh elements: ', len(mesh.elements)

dt = 1/75.
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
E_order = (order, True)
B_order = (order, True)
freeE = freeB = cb
system = PMLSystem(mesh, disc_orders=dict(E=E_order, B=B_order),
                   BCs=Struct(E=freeE, B=freeB))
totalDOFs = system.discs.E.totalDOFs

print "Getting source DOFs"
weights, elPerm = system.dofs.E.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=N.array([0,0,0.]))
drive_dofnos = elPerm[1]
print "Done"
print "Getting test points"
direction = N.array([0, 1., 0])
test_pts = N.array([[0, 0.5875, 0]])
test_elnos, test_el_coords = PostProc.LocatePoints(
    system.mesh, test_pts)
print "Done"
system.addReconstructedLogger('E', test_elnos, test_el_coords)
system.set_sigmas(sigma_fns)

current_waveform = Waveforms.get_d_gaussian(fc=1, tpr=-60)
drv_fun = current_waveform             
system.setDriveDOFs_J(drive_dofnos, weights, drv_fun)
print "Setting timestep"
system.setTimestep(dt)
print "Done"
system.step(int(N.ceil(runtime/dt)))
import pickle

from pylab import plot, xlabel, ylabel, legend, show

rs = system.loggedReconstructed.E[0]
rsvals = rs['vals']
ts = N.array([ts[0] for ts in rsvals])
plot(N.arange(len(ts))*dt, ts[:,0], label='E_x')
plot(N.arange(len(ts))*dt, ts[:,1], label='E_y')
plot(N.arange(len(ts))*dt, ts[:,2], label='E_z')
xlabel('time (s)')
ylabel('E field magnitude at r = [0, 0.5875, 0]')
legend()
