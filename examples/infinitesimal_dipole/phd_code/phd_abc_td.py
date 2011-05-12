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
sys.path.append('../../../')
import NewCode
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.DiscretisedSystem import CurlCurlNewmark
from NewCode.Meshes import CalculateConnectivity
from NewCode.Meshes import Conversions
from NewCode.Meshes.MeshIO import Femmesh

#meshfile = '../../../workspace/sphere-r1m-6.femmesh'
# listmesh = CalculateConnectivity.get_all_connectivities(
#     Femmesh.get_femmesh_as_listmesh(meshfile))

source_coord = N.array([0,0,0.])

## Set up dolfin mesh
import dolfin as dol
source_point = dol.Point(*source_coord)
domain_size = N.array([4.]*3)
max_edge_len = 1/6.
domain_subdivisions = N.array(N.ceil(domain_size/max_edge_len), N.uint)
dol_mesh = dol.UnitCube(*domain_subdivisions)
# Transform dol_mesh to correct dimensions
dol_mesh.coordinates()[:] *= domain_size
dol_mesh.coordinates()[:] -= domain_size/2
## Translate dol_mesh slightly so that source coordinate lies at centroid of an element
source_elnos = dol_mesh.all_intersected_entities(source_point)
closest_elno = source_elnos[(N.argmin([source_point.distance(dol.Cell(dol_mesh, i).midpoint())
                                  for i in source_elnos]))]
centre_pt = dol.Cell(dol_mesh, closest_elno).midpoint()
centre_coord = N.array([centre_pt.x(), centre_pt.y(), centre_pt.z()])
# There seems to be an issue with the intersect operator if the
# dol_mesh coordinates are changed after calling it for the first
# time. Since we called it to find the centroid, we should init a
# new dol_mesh
dol_mesh_coords = dol_mesh.coordinates().copy()
dol_mesh = dol.UnitCube(*domain_subdivisions)
dol_mesh.coordinates()[:] = dol_mesh_coords
dol_mesh.coordinates()[:] -= centre_coord
listmesh = Conversions.dolfin_mesh_2_listmesh(dol_mesh)
listmesh = CalculateConnectivity.get_all_connectivities(listmesh)
## End dolfin mesh setup

mesh = Mesh.Mesh(listmesh)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 6
CurlCurlNewmark.triIntegrationOrder = 6

print 'Mesh elements: ', len(mesh.elements)

dt = 1/30.
runtime = 6
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
order = 1
newmarkSystem = CurlCurlNewmark(mesh, order=order, BC=free, useQ=True)
totalDOFs = newmarkSystem.disc.totalDOFs

print "Getting source DOFs"
weights, elPerm = newmarkSystem.dofs.calcProjPointfunRHS_with_elPerm(
    matchfun=lambda r: N.array([0,0,1.]), r0=source_coord)
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
#             file('+sphere-r1-LTQN-ABC-runstuff.pickle', 'w'))

from pylab import plot, xlabel, ylabel, legend, show
rs = newmarkSystem.loggedReconstructed[0]
rsvals = rs['vals']
ts = N.array([ts[0] for ts in rsvals])
plot(N.arange(len(ts))*dt, ts[:,0], label='E_x')
plot(N.arange(len(ts))*dt, ts[:,1], label='E_y')
plot(N.arange(len(ts))*dt, ts[:,2], label='E_z')
xlabel('time (s)')
ylabel('E field magnitude at r = [0, 0.5875, 0]')
legend()
