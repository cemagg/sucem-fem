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
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import sys
import numpy as N
import os
import dolfin
sys.path.insert(0, '../../')
import sucemfem.Sources.current_source
import sucemfem.BoundaryConditions.ABC
import sucemfem.Utilities.LinalgSolvers
import sucemfem.Utilities.Optimization
from sucemfem.Utilities.MeshGenerators import get_centred_cube
from sucemfem.Consts import eps0, mu0, c0
from sucemfem.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from sucemfem.Utilities.LinalgSolvers import solve_sparse_system
from sucemfem.Sources.fillament_current_source import FillamentCurrentSource
from sucemfem.PostProcessing import surface_ntff
from sucemfem.Testing.ErrorMeasures import normalised_RMS
import pylab
from sucemfem.Testing.Analytical import current_fillament_farfield
del sys.path[0]


sucemfem.Utilities.Optimization.set_dolfin_optimisation(True)
### Postprocessing requests
theta_deg = N.linspace(0, 180, 181)
no_ff_pts = len(theta_deg)
phi_deg = N.zeros(no_ff_pts)


### Problem parameters
freq =  1.0e+9                          # Frequency
lam = c0/freq
l = lam/10                            # Dipole length
I = 1.0                                 # Dipole current
source_direction = N.array([0,0,1.])    # Source orientation
source_centre = N.array([0,0,0.])        # Position of the source
source_endpoints =  N.array(
    [-source_direction*l/2, source_direction*l/2]) + source_centre

### Discretisation settings
order = 2
domain_size = N.array([lam]*3)*0.25
max_edge_len = lam/36
mesh = get_centred_cube(domain_size, max_edge_len)

fname = 'data/f-%f_o-%d_s-%f_l-%f_h-%f' % (freq, order,
                                           domain_size[0], l/lam,
                                           max_edge_len/lam)
meshfilename = fname+'_mesh.xml'
materialsfilename = fname+'_materials.xml'
print fname
### Implementation
#
## Set up materials function with all free-space
material_mesh_func = dolfin.MeshFunction('uint', mesh, 3)
material_mesh_func.set_all(0)
materials = {0:dict(eps_r=1, mu_r=1),}
## Set up 1st-order analytical ABC
abc = sucemfem.BoundaryConditions.ABC.ABCBoundaryCondition()
abc.set_region_number(1)
bcs = sucemfem.BoundaryConditions.container.BoundaryConditions()
bcs.add_boundary_condition(abc)
## Set up high level problem class
dp = DrivenProblemABC()
dp.set_mesh(mesh)
dp.set_basis_order(order)
dp.set_material_regions(materials)
dp.set_region_meshfunction(material_mesh_func)
dp.set_boundary_conditions(bcs)
## Set up current fillament source
current_sources = sucemfem.Sources.current_source.CurrentSources()
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
umf_solver = sucemfem.Utilities.LinalgSolvers.UMFPACKSolver(A)
x = umf_solver.solve(b)

import pickle
with open(fname+'.pickle', 'w') as f:
    pickle.dump(dict(x=x, meshfile=meshfilename,
                     materialsfile=materialsfilename,
                     material_properties=materials,
                     freq=freq, order=order, 
                     source_endpoints=source_endpoints,
                     I=I),f)
dolfin.File(meshfilename)  << mesh
dolfin.File(materialsfilename)  << material_mesh_func


# ## Post-process solution to obtain far-field
# print 'calculating far field'
# surf_ntff = surface_ntff.NTFF(dp.function_space)
# surf_ntff.set_dofs(x)
# surf_ntff.set_frequency(freq)
# surf_E_ff = N.array([surf_ntff.calc_pt(th_deg, ph_deg)
#                 for th_deg, ph_deg in zip(theta_deg, phi_deg)])
# surf_E_theta = surf_E_ff[:,0]
# surf_E_phi = surf_E_ff[:,1]

# ## Calculate some errors relative to the analytical solution
# an_E_theta = [current_fillament_farfield.eval_E_theta(freq, l, I, th)
#               for th in N.deg2rad(theta_deg)]
# start=10 ; stop=-10                     # Don't include the very ends
# err = normalised_RMS(
#     surf_E_theta[start:stop], an_E_theta[start:stop], surf_E_phi[start:stop])
# err_theta = normalised_RMS(surf_E_theta[start:stop], an_E_theta[start:stop])
# err_abs_theta = normalised_RMS(N.abs(surf_E_theta[start:stop]),
#                                N.abs(an_E_theta[start:stop]))
# print 'Far-field RMS error: ', err

# print 'plotting'
# pylab.figure()
# pylab.plot(theta_deg, N.abs(surf_E_theta), label='|E_theta|')
# pylab.plot(theta_deg, N.abs(surf_E_phi), label='|E_phi|')
# an_E_theta = [sucemfem.Testing.Analytical.current_fillament_farfield.eval_E_theta(freq, l, I, th) for th in N.deg2rad(theta_deg)]
# pylab.plot(theta_deg, N.abs(an_E_theta), label='analytical')
# pylab.legend()



