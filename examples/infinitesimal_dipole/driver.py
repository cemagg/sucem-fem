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
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from FenicsCode.PostProcessing import Reconstruct

sys.path.insert(0, '../../')
from test_data import problem_data
del sys.path[0]

## Problem parameters
I = problem_data['I']                   # Dipole current
l = problem_data['l']                   # Dipole length
source_value = N.array([0,0,1.])*I*l
freq = problem_data['f']
lam = c0/freq
source_coord = N.array([0,0,0.]) 
## Discretisation settings
order = 2
domain_size = N.array([2*lam]*3)
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len, source_coord)
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

from pylab import *
r1 = field_pts[:]/lam
x1 = r1[:,0]
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
