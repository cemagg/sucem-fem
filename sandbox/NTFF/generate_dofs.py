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
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system, UMFPACKSolver
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
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
domain_size = N.array([lam]*3)*0.5
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len, source_coord)
# Request information:
field_pts = N.array(
    [N.max(domain_size)/2-max_edge_len/2,0,0]
    )*(N.arange(88)/100+1/10)[:, N.newaxis]

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
#print 'solve using scipy bicgstab'
#x = solve_sparse_system ( A, b, preconditioner_type='diagonal' )
print 'solve using scipy UMFPACK'
x = UMFPACKSolver(A).solve(b)

recon = Reconstruct(dp.function_space)
recon.set_dof_values(x)
E_field = recon.reconstruct_points(field_pts)

import pickle
fname = 'dofs-%d-%s-%s.pickle' % (order, str(domain_size[0]), str(max_edge_len))
pickle.dump(
    dict(x=x, order=order, domain_size=domain_size,
         max_edge_len=max_edge_len, freq=freq),
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
