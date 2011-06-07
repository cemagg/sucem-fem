from __future__ import division

import sys
sys.path.append('../../')
import numpy as N
import os
import dolfin
import FenicsCode
import FenicsCode.BoundaryConditions.ABC
from FenicsCode.Consts import eps0, mu0, c0
from FenicsCode.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from FenicsCode.Sources import point_source
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system
from FenicsCode.Utilities.MeshIO import femmesh_2_dolfin_mesh
from FenicsCode.PostProcessing import Reconstruct


## Problem parameters
freq = 1e9
lam = c0/freq
I = 1.                                  # Dipole current
l = lam/1000                            # Dipole length
source_value = N.array([0,0,1.])*I*l
source_coord = N.array([0,0,0.]) 
## Discretisation settings
order = 2
mesh_id_list = [
    'sphere-r1m-4',
    'sphere-r1m-5',
    'sphere-r1m-6',
    'sphere-r1m-7',
    'sphere-r1m-8',
    'sphere-r1m-9',
    'sphere-r1m-10',
    'sphere-r1m-15',
    ]
mesh_id = mesh_id_list[-1]
mesh_file = '../solvers/meshes/%s.femmesh' % mesh_id
mesh = femmesh_2_dolfin_mesh(mesh_file)
mesh.coordinates()[:] *= 2*lam
# Request information:
field_pts = N.array([lam,0,0])*(N.arange(88)/100+1/10)[:, N.newaxis]

## Implementation
material_mesh_func = dolfin.MeshFunction('uint', mesh, 3)
material_mesh_func.set_all(0)
materials = {0:dict(eps_r=1, mu_r=1),}
abc = FenicsCode.BoundaryConditions.ABC.ABCBoundaryCondition()
abc.set_region_number(1)
bcs = FenicsCode.BoundaryConditions.container.BoundaryConditions()
bcs.add_boundary_condition(abc)
dp = DrivenProblemABC()
dp.set_mesh(mesh)
dp.set_basis_order(order)
dp.set_material_regions(materials)
dp.set_region_meshfunction(material_mesh_func)
dp.set_boundary_conditions(bcs)
current_sources = point_source.CurrentSources()
dipole_source = point_source.PointCurrentSource()
dipole_source.set_position(source_coord)
dipole_source.set_value(source_value)
current_sources.add_source(dipole_source)
dp.set_sources(current_sources)
dp.init_problem()
dp.set_frequency(freq)

A = dp.get_LHS_matrix()
b = dp.get_RHS()
print 'solve using scipy bicgstab'
x = solve_sparse_system ( A, b, preconditioner_type='diagonal' )

recon = Reconstruct(dp.function_space)
recon.set_dof_values(x)
E_field = recon.reconstruct_points(field_pts)

import pickle
fname = 'dofs_sph-2lam-%d-%s.pickle' % (order, mesh_id)
pickle.dump(
    dict(x=x, order=order, mesh_id=mesh_id, freq=freq),
    open(fname, 'w'))

from pylab import *
r1 = field_pts[:]/lam
x1 = r1[:,0]
sys.path.append('../../examples/infinitesimal_dipole')
from blog_example import analytical_pts, analytical_result
E_ana = N.abs(analytical_result)
E_num = E_field
figure()
plot(x1, N.abs(E_num[:,0]), '-g', label='x_num')
plot(x1, N.abs(E_num[:,1]), '-b', label='y_num')
plot(x1, N.abs(E_num[:,2]), '-r', label='z_num')
plot(analytical_pts, E_ana, '--r', label='z_ana')
ylabel('E-field Magnitude')
xlabel('Distance (wavelengths)')
legend(loc='best')
grid(True)
show()
