## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import sys
import numpy as N
import os
import dolfin
sys.path.insert(0, '../../')
import sucemfem.Sources.current_source
from sucemfem.BoundaryConditions import ABCBoundaryCondition, BoundaryConditions
from sucemfem.Consts import eps0, mu0, c0
from sucemfem.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from sucemfem.Sources import point_source
from sucemfem.Utilities.LinalgSolvers import solve_sparse_system
from sucemfem.Utilities.MeshGenerators import get_centred_cube
from sucemfem.PostProcessing import Reconstruct
from test_data import problem_data
del sys.path[0]

from pylab import *
from blog_example import analytical_pts, analytical_result

## Problem parameters
I = problem_data['I']                   # Dipole current
l = problem_data['l']                   # Dipole length
source_value = N.array([0,0,1.])*I*l
freq = problem_data['f']
lam = c0/freq
source_coord = N.array([0,0,0.]) 
## Discretisation settings
order = 1
domain_size = N.array([2*lam]*3)
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len, source_coord)
# Request information:
field_pts = N.array([lam,0,0])*(N.arange(88)/100+1/10)[:, N.newaxis]

## Implementation
material_mesh_func = dolfin.MeshFunction('uint', mesh, 3)
material_mesh_func.set_all(0)
materials = {0:dict(eps_r=1, mu_r=1),}
abc = ABCBoundaryCondition()
abc.set_region_number(1)
bcs = BoundaryConditions()
bcs.add_boundary_condition(abc)
dp = DrivenProblemABC()
dp.set_mesh(mesh)
dp.set_basis_order(order)
dp.set_material_regions(materials)
dp.set_region_meshfunction(material_mesh_func)
dp.set_boundary_conditions(bcs)
current_sources = sucemfem.Sources.current_source.CurrentSources()
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

r1 = field_pts[:]/lam
x1 = r1[:,0]
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
