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
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.DiscretisedSystem import PMLSystem

h = 1.000000000001/7.5
order = 1
no_PML_cells = 5
geom_size = N.array([6,3,3.])*h
dipole_pt = N.array([-1.5,0,0])*h
dipole_dir = N.array([0,0,1.])


#geom_size = N.array([2]*3)
g_a, g_b, g_c = geom_size
geom_max = geom_size/2
geom_min = -geom_max

meas_max = geom_size/2*(1-1e-10)
meas_min = -meas_max


meshinfo_no_PML = BrickMeshGen.make_rect_cavity_brick_listmesh(
    g_a, g_b, g_c, [h, h, h], grid_offset=[-g_a/2,-g_b/2,-g_c/2])
real_h = meshinfo_no_PML['GridStepSize'] 


PML_m = 3
sigma_factor = 1
runtime = 6.
measure_x, measure_y, measure_z = N.vstack((meas_min, [1e-10]*3, meas_max)).T
#measure_x *= .95
#measure_y *= 0.95 ; measure_z *= .95
measure_pts = N.array(
    [(x,y,z) for x in measure_x for y in measure_y for z in measure_z])

PML_len = real_h*no_PML_cells
a,b,c = geom_size + PML_len*2


mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(
    a, b, c, real_h, grid_offset=[-a/2,-b/2,-c/2]))

PML_pos = N.array([a,b,c])/2 - real_h*no_PML_cells

wave_imp = 1.
PML_sigma_max = 0.8*(PML_m+1)/real_h/wave_imp*sigma_factor

def sigma_u(coord_no, r) :
    u = r.T[coord_no] ; sigma_max = PML_sigma_max[coord_no]
    pos = PML_pos[coord_no] ; plen = PML_len[coord_no]
    return N.where(N.abs(u) >= pos,
                   sigma_max*((N.abs(u)-pos)/plen)**PML_m, 0.)
    
sigma_x_fn = partial(sigma_u,0)
sigma_y_fn = partial(sigma_u,1)
sigma_z_fn = partial(sigma_u,2)
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
    matchfun=lambda r: dipole_dir, r0=dipole_pt)
drive_dofnos = elPerm[1]
print "Done"
print "Getting test points"
test_elnos, test_el_coords = PostProc.LocatePoints(
    system.mesh, measure_pts)
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
#pickle.dump(Struct(loggedReconstructed=system.loggedReconstructed),
#            file('+pml_%dcell_dipole-refltest_3rd.pickle' % no_PML_cells, 'w'))


#from pylab import *
rs = system.loggedReconstructed.E[0]
rsvals = rs['vals']
look_pts = (4,-5)                       # [+-xmax, 0, 0]
tsa, tsb = N.array([[ts[d] for ts in rsvals] for d in look_pts])
figure(1)
plot(N.arange(len(tsa))*dt, tsa[:,0], label='E_xa')
plot(N.arange(len(tsa))*dt, tsa[:,1], label='E_ya')
plot(N.arange(len(tsa))*dt, tsa[:,2], label='E_za')
xlabel('time (s)')
ylabel('E field magnitude at r = %s' % str(measure_pts[look_pts[0]]))
legend()
title('order %d using %d PML cells' %(order, no_PML_cells))
figure(2)
plot(N.arange(len(tsb))*dt, tsb[:,0], label='E_xb')
plot(N.arange(len(tsb))*dt, tsb[:,1], label='E_yb')
plot(N.arange(len(tsb))*dt, tsb[:,2], label='E_zb')
xlabel('time (s)')
ylabel('E field magnitude at r = %s' % str(measure_pts[look_pts[1]]))
legend()
title('order %d using %d PML cells' %(order, no_PML_cells))

