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
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 

import dolfin
import numpy as np
import sys
modpath = '../../'
sys.path.append(modpath)

import sucemfem
from sucemfem.BoundaryConditions import ABCBoundaryCondition, BoundaryConditions
from sucemfem.BoundaryConditions import PECWallsBoundaryCondition 
from sucemfem.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from sucemfem.Sources.fillament_current_source import FillamentCurrentSource
import sucemfem.Utilities.LinalgSolvers
from sucemfem.Sources.PostProcess import ComplexVoltageAlongLine

sys.path.remove(modpath)

## Problem specifications
meshname = 'patch'
freqs = np.linspace(2.75, 3.2, 17)*1e9
er = 2.2
materials = {0:dict(eps_r=1),
             20:dict(eps_r=er)}
PEC_label = 10
If = 8.9e-3                             # feed inset
### Consider reading these from the .geo file using the suggestions in
### http://article.gmane.org/gmane.comp.cad.gmsh.general/3381
Lx = 31.18e-3;                          # Patch x-length from patch.geo
Hs = 2.87e-3;                           # Substrate height from patch.geo
Of = Lx/2 - If                         # Feed offeset from centre
feed_pts = np.array([[-Of, 0, 0], [-Of, 0, Hs]])
I = 1.                                  # Source current
##
## Solution settings
order = 1                               # Discretisation order
solver_types = ['umfpack', 'bicgstab', 'gmres']
solver_type = solver_types[1]
##
## Calculations
mesh = dolfin.Mesh(meshname+'.xml')
materials_mesh_function = dolfin.MeshFunction('uint', mesh, meshname + '_physical_region.xml')
pec_mesh_function = dolfin.MeshFunction('uint', mesh, meshname + '_facet_region.xml')
mesh.init(3,1)
print '%d elements with %d edges' % (mesh.num_cells(), mesh.num_edges())

abc = ABCBoundaryCondition()
abc.set_region_number(1)
bcs = BoundaryConditions()
bcs.add_boundary_condition(abc)
pec_walls = PECWallsBoundaryCondition()
pec_walls.init_with_meshfunction (pec_mesh_function, PEC_label)
bcs.add_boundary_condition(pec_walls)

## Set up high level problem class
dp = DrivenProblemABC()
dp.set_mesh(mesh)
dp.set_basis_order(order)
dp.set_material_regions(materials)
dp.set_region_meshfunction(materials_mesh_function)
dp.set_boundary_conditions(bcs)
## Set up current fillament source
current_sources = sucemfem.Sources.current_source.CurrentSources()
fillament_source = FillamentCurrentSource()
fillament_source.set_source_endpoints(feed_pts)
fillament_source.set_value(I)
current_sources.add_source(fillament_source)
## Set source in problem container
dp.set_sources(current_sources)
print 'init problem'
dp.init_problem()
Zs = []

def get_solver(type, A):
    ## Solve. Choose sparse solver if UMFPack runsout of memory
    if type == 'umfpack':
        print 'solve using UMFPack'
        solver = sucemfem.Utilities.LinalgSolvers.UMFPACKSolver(A)
    elif type == 'bicgstab':
        print 'solve using scipy bicgstab'
        solver = sucemfem.Utilities.LinalgSolvers.BiCGStabSolver(
            A, preconditioner_type='diagonal')
    elif type == 'gmres':
        print 'solve using scipy gmres'
        solver = sucemfem.Utilities.LinalgSolvers.GMRESSolver(
            A, preconditioner_type='diagonal')
    return solver
        


for freq in freqs:
    print 'set frequency %f' %freq
    dp.set_frequency(freq)
    print 'gettix LHS matrix'
    ## Get sytem LHS matrix and RHS Vector
    A = dp.get_LHS_matrix()
    print 'gettix RHS'
    b = dp.get_RHS()
    solver = get_solver(solver_type, A)
    def cb(solverobj, res):
        if solverobj._callback_count % 1001 == 1:
            print 'Frequency %f residual at iteration %d: %f' % (
                freq, solverobj._callback_count, res)
    solver.add_user_callback(cb)
    x = solver.solve(b)
    complex_voltage = ComplexVoltageAlongLine(dp.function_space)
    complex_voltage.set_dofs(x)

    volts = -complex_voltage.calculate_voltage(*feed_pts)
    Z = volts/I
    print 'impedance at %f Hz: ', Z
    Zs.append(Z)

import pickle
pickle.dump(dict(Zs=Zs, freqs=freqs), open('patch_o-%d_%s_Z.pickle' % (order, solver_type), 'w'))
