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
import setup_mesh
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
mesh = setup_mesh.get_centred_cube(domain_size, max_edge_len, source_coord)

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

dp.init_problem()
