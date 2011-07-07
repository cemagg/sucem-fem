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
import FenicsCode.Utilities.LinalgSolvers

from FenicsCode.Utilities.MeshGenerators import get_centred_cube
from FenicsCode.PostProcessing import Reconstruct
import FenicsCode.Utilities.Optimization
del sys.path[0]
from FillamentSource import FillamentCurrentSource


FenicsCode.Utilities.Optimization.set_dolfin_optimisation(True)
## Postprocessing requests
theta_deg = N.linspace(0, 180, 181)
no_ff_pts = len(theta_deg)
phi_deg = N.zeros(no_ff_pts)


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
order = 3
domain_size = N.array([lam]*3)*1
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
dipole_source.no_integration_points = 1000
dipole_source.set_source_endpoints(source_endpoints)
dipole_source.set_value(I)
current_sources.add_source(dipole_source)
dp.set_sources(current_sources)
dp.init_problem()
dp.set_frequency(freq)

A = dp.get_LHS_matrix()
b = dp.get_RHS()
print 'solve using scipy bicgstab'
x = solve_sparse_system ( A, b, preconditioner_type='diagonal')
# print 'solve using UMFPack'
# umf_solver = FenicsCode.Utilities.LinalgSolvers.UMFPACKSolver(A)
# x = umf_solver.solve(b)

import pickle
fname = 'data/dofs-%d-%s-%s-%s-%s.pickle' % (
    order, str(domain_size[0]), str(max_edge_len), str(l),
    str(dipole_source.no_integration_points))
pickle.dump(
    dict(x=x, order=order, domain_size=domain_size,
         max_edge_len=max_edge_len, freq=freq, l=l, I=I,
         source_integration_points=dipole_source.no_integration_points),
    open(fname, 'w'))



print 'calculating far field'
from FenicsCode.PostProcessing import surface_ntff
surf_ntff = surface_ntff.NTFF(dp.function_space)
surf_ntff.set_dofs(x)
surf_ntff.set_frequency(freq)
surf_E_ff = N.array([surf_ntff.calc_pt(th_deg, ph_deg)
                for th_deg, ph_deg in zip(theta_deg, phi_deg)])
surf_E_theta = surf_E_ff[:,0]
surf_E_phi = surf_E_ff[:,1]



print 'plotting'
import pylab
pylab.figure()
pylab.plot(theta_deg, N.abs(surf_E_theta), label='|E_theta|')
pylab.plot(theta_deg, N.abs(surf_E_phi), label='|E_phi|')
import analytical
an_E_theta = [analytical.eval_E_theta(freq, l, I, th) for th in N.deg2rad(theta_deg)]
pylab.plot(theta_deg, N.abs(an_E_theta), label='analytical')
pylab.legend()
start=10 ; stop=-10
from FenicsCode.Testing.ErrorMeasures import normalised_RMS

err = normalised_RMS(
    surf_E_theta[start:stop], an_E_theta[start:stop], surf_E_phi[start:stop])
err_theta = normalised_RMS(surf_E_theta[start:stop], an_E_theta[start:stop])
err_abs_theta = normalised_RMS(N.abs(surf_E_theta[start:stop]),
                               N.abs(an_E_theta[start:stop]))
print fname, err, err_theta

# recon = Reconstruct(dp.function_space)
# recon.set_dof_values(x)
# E_field = recon.reconstruct_points(field_pts)

# from pylab import *

# r1 = field_pts[:]/lam
# x1 = r1[:,0]
# #E_ana = N.abs(analytical_result)
# E_num = E_field
# figure()
# plot(x1, N.abs(E_num[:,0]), '-g', label='x_num')
# plot(x1, N.abs(E_num[:,1]), '-b', label='y_num')
# plot(x1, N.abs(E_num[:,2]), '-r', label='z_num')
# #plot(analytical_pts, E_ana, '--r', label='z_ana')
# ylabel('E-field Magnitude')
# xlabel('Distance (wavelengths)')
# legend(loc='best')
# grid(True)
# show()
