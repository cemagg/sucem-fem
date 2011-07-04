# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import sys
import numpy as N
import os
import dolfin
sys.path.insert(0, '../../')
import FenicsCode.Sources.current_source
import FenicsCode.BoundaryConditions.ABC
from FenicsCode.Consts import eps0, mu0, c0
from FenicsCode.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from FenicsCode.PostProcessing import Reconstruct
import FenicsCode.Utilities.Optimization
del sys.path[0]
from FillamentSource import FillamentCurrentSource


FenicsCode.Utilities.Optimization.set_dolfin_optimisation(True)
## Problem parameters
freq =  1.0e+9                          # Frequency
lam = c0/freq
l = lam/4                               # Dipole length
I = 1.0                                 # Dipole current
source_direction = N.array([0,0,1.])    # Source orientation
source_centre = N.array([0,0,0.])        # Position of the source
source_endpoints =  N.array(
    [-source_direction*l/2, source_direction*l/2]) + source_centre

## Discretisation settings
order = 1
domain_size = N.array([lam]*3)
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len)
# Request information:
field_pts = N.array([1,0,0])*(N.linspace(0, domain_size[0]/2.1))[:, N.newaxis]

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
current_sources = FenicsCode.Sources.current_source.CurrentSources()
dipole_source = FillamentCurrentSource()
dipole_source.set_source_endpoints(source_endpoints)
dipole_source.set_value(I)
current_sources.add_source(dipole_source)
dp.set_sources(current_sources)
dp.init_problem()
dp.set_frequency(freq)

A = dp.get_LHS_matrix()
b = dp.get_RHS()
print 'solve using scipy bicgstab'
x = solve_sparse_system ( A, b)

recon = Reconstruct(dp.function_space)
recon.set_dof_values(x)
E_field = recon.reconstruct_points(field_pts)

from pylab import *

r1 = field_pts[:]/lam
x1 = r1[:,0]
#E_ana = N.abs(analytical_result)
E_num = E_field
figure()
plot(x1, N.abs(E_num[:,0]), '-g', label='x_num')
plot(x1, N.abs(E_num[:,1]), '-b', label='y_num')
plot(x1, N.abs(E_num[:,2]), '-r', label='z_num')
#plot(analytical_pts, E_ana, '--r', label='z_ana')
ylabel('E-field Magnitude')
xlabel('Distance (wavelengths)')
legend(loc='best')
grid(True)
show()
