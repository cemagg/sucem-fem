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
import FenicsCode.Utilities.LinalgSolvers
import FenicsCode.Utilities.Optimization
from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from FenicsCode.Consts import eps0, mu0, c0
from FenicsCode.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from FenicsCode.Utilities.LinalgSolvers import solve_sparse_system
from FenicsCode.Sources.fillament_source import FillamentCurrentSource
from FenicsCode.PostProcessing import surface_ntff
from FenicsCode.Testing.ErrorMeasures import normalised_RMS
import pylab
from FenicsCode.Testing.Analytical import current_fillament_farfield
del sys.path[0]


FenicsCode.Utilities.Optimization.set_dolfin_optimisation(True)
### Postprocessing requests
theta_deg = N.linspace(0, 180, 181)
no_ff_pts = len(theta_deg)
phi_deg = N.zeros(no_ff_pts)


### Problem parameters
freq =  1.0e+9                          # Frequency
lam = c0/freq
l = lam/4                               # Dipole length
I = 1.0                                 # Dipole current
source_direction = N.array([0,0,1.])    # Source orientation
source_centre = N.array([0,0,0.])        # Position of the source
source_endpoints =  N.array(
    [-source_direction*l/2, source_direction*l/2]) + source_centre

### Discretisation settings
order = 2
domain_size = N.array([lam]*3)*1
max_edge_len = lam/6
mesh = get_centred_cube(domain_size, max_edge_len)

### Implementation
#
## Set up materials function with all free-space
material_mesh_func = dolfin.MeshFunction('uint', mesh, 3)
material_mesh_func.set_all(0)
materials = {0:dict(eps_r=1, mu_r=1),}
## Set up 1st-order analytical ABC
abc = FenicsCode.BoundaryConditions.ABC.ABCBoundaryCondition()
abc.set_region_number(1)
bcs = FenicsCode.BoundaryConditions.container.BoundaryConditions()
bcs.add_boundary_condition(abc)
## Set up high level problem class
dp = DrivenProblemABC()
dp.set_mesh(mesh)
dp.set_basis_order(order)
dp.set_material_regions(materials)
dp.set_region_meshfunction(material_mesh_func)
dp.set_boundary_conditions(bcs)
## Set up current fillament source
current_sources = FenicsCode.Sources.current_source.CurrentSources()
fillament_source = FillamentCurrentSource()
fillament_source.set_source_endpoints(source_endpoints)
fillament_source.set_value(I)
current_sources.add_source(fillament_source)
## Set source in problem container
dp.set_sources(current_sources)
dp.init_problem()
dp.set_frequency(freq)

## Get sytem LHS matrix and RHS Vector
A = dp.get_LHS_matrix()
b = dp.get_RHS()
## Solve. Choose spare solver if UMFPack runsout of memory
# print 'solve using scipy bicgstab'
# x = solve_sparse_system ( A, b, preconditioner_type='diagonal')
print 'solve using UMFPack'
umf_solver = FenicsCode.Utilities.LinalgSolvers.UMFPACKSolver(A)
x = umf_solver.solve(b)

## Post-process solution to obtain far-field
print 'calculating far field'
surf_ntff = surface_ntff.NTFF(dp.function_space)
surf_ntff.set_dofs(x)
surf_ntff.set_frequency(freq)
surf_E_ff = N.array([surf_ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
surf_E_theta = surf_E_ff[:,0]
surf_E_phi = surf_E_ff[:,1]

## Calculate some errors relative to the analytical solution
an_E_theta = [current_fillament_farfield.eval_E_theta(freq, l, I, th)
              for th in N.deg2rad(theta_deg)]
start=10 ; stop=-10                     # Don't include the very ends
err = normalised_RMS(
    surf_E_theta[start:stop], an_E_theta[start:stop], surf_E_phi[start:stop])
err_theta = normalised_RMS(surf_E_theta[start:stop], an_E_theta[start:stop])
err_abs_theta = normalised_RMS(N.abs(surf_E_theta[start:stop]),
                               N.abs(an_E_theta[start:stop]))
print 'Far-field RMS error: ', err

print 'plotting'
pylab.figure()
pylab.plot(theta_deg, N.abs(surf_E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(surf_E_phi), label='|E_phi|')
an_E_theta = [FenicsCode.Testing.Analytical.current_fillament_farfield.eval_E_theta(freq, l, I, th) for th in N.deg2rad(theta_deg)]
pylab.plot(theta_deg, N.abs(an_E_theta), label='analytical')
pylab.legend()



