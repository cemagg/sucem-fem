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
"""
This is a functional test for the far-field solution of a constant-current fillament
"""

import numpy as N
import dolfin
import unittest
import sucemfem.Sources.current_source
import sucemfem.Utilities.LinalgSolvers
import sucemfem.Utilities.Optimization
from sucemfem.BoundaryConditions import ABCBoundaryCondition, BoundaryConditions
from sucemfem.Utilities.MeshGenerators import get_centred_cube
from sucemfem.Consts import eps0, mu0, c0
from sucemfem.ProblemConfigurations.EMDrivenProblem import DrivenProblemABC
from sucemfem.Sources.fillament_current_source import FillamentCurrentSource
from sucemfem.PostProcessing import surface_ntff
from sucemfem.Testing.ErrorMeasures import normalised_RMS
from sucemfem.Testing.Analytical import current_fillament_farfield

class test_current_fillament(unittest.TestCase):
    def test_ff_error(self):
        sucemfem.Utilities.Optimization.set_dolfin_optimisation(True)
        ### Postprocessing requests
        theta_deg = N.linspace(10, 170, 161)
        no_ff_pts = len(theta_deg)
        phi_deg = N.zeros(no_ff_pts)
        ### Problem parameters
        freq =  1.0e+9                          # Frequency
        lam = c0/freq
        l = lam/4                               # Dipole length
        I = 1.0                                 # Dipole current
        source_direction_z = N.array([0,0,1.])    # Source orientation
        source_direction_x = N.array([1.,0,0])    # Source orientation
        source_direction_y = N.array([0,1,0.])    # Source orientation
        source_centre = N.array([0,0,0.])        # Position of the source
        source_endpoints_z =  N.array(
            [-source_direction_z*l/2, source_direction_z*l/2]) + source_centre
        source_endpoints_x =  N.array(
            [-source_direction_x*l/2, source_direction_x*l/2]) + source_centre
        source_endpoints_y =  N.array(
            [-source_direction_y*l/2, source_direction_y*l/2]) + source_centre

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
        abc = ABCBoundaryCondition()
        abc.set_region_number(1)
        bcs = BoundaryConditions()
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
        fillament_source.no_integration_points = 1000
        fillament_source.set_source_endpoints(source_endpoints_z)
        fillament_source.set_value(I)
        current_sources.add_source(fillament_source)
        ## Set source in problem container
        dp.set_sources(current_sources)
        dp.init_problem()
        dp.set_frequency(freq)

        ## Get sytem LHS matrix and RHS Vector
        A = dp.get_LHS_matrix()
        b_z = dp.get_RHS()
        fillament_source.set_source_endpoints(source_endpoints_x)
        b_x = dp.get_RHS()
        fillament_source.set_source_endpoints(source_endpoints_y)
        b_y = dp.get_RHS()

        #import pdb ; pdb.set_trace()
        A
        print 'solve using UMFPack'
        umf_solver = sucemfem.Utilities.LinalgSolvers.UMFPACKSolver(A)
        x_z = umf_solver.solve(b_z)
        x_x = umf_solver.solve(b_x)
        x_y = umf_solver.solve(b_y)

        ## Post-process solution to obtain far-field
        print 'calculating far field'
        surf_ntff = surface_ntff.NTFF(dp.function_space)
        surf_ntff.set_dofs(x_z)
        surf_ntff.set_frequency(freq)
        surf_E_ff_z = N.array([surf_ntff.calc_pt(th_deg, ph_deg)
                        for th_deg, ph_deg in zip(theta_deg, phi_deg)])
        surf_E_theta_z = surf_E_ff_z[:,0]
        surf_E_phi_z = surf_E_ff_z[:,1]
        surf_ntff.set_dofs(x_x)
        surf_E_ff_x = N.array([surf_ntff.calc_pt(th_deg+90, ph_deg)
                        for th_deg, ph_deg in zip(theta_deg, phi_deg)])
        surf_E_theta_x = surf_E_ff_x[:,0]
        surf_E_phi_x = surf_E_ff_x[:,1]
        surf_ntff.set_dofs(x_y)
        surf_E_ff_y = N.array([surf_ntff.calc_pt(th_deg+90, ph_deg)
                        for th_deg, ph_deg in zip(theta_deg, phi_deg)])
        surf_E_theta_y = surf_E_ff_y[:,0]
        surf_E_phi_y = surf_E_ff_y[:,1]

        ## Calculate some errors relative to the analytical solution
        an_E_theta = [current_fillament_farfield.eval_E_theta(freq, l, I, th)
                      for th in N.deg2rad(theta_deg)]
        err_z = normalised_RMS(
            surf_E_theta_z, an_E_theta, surf_E_phi_z)
        err_theta_z = normalised_RMS(surf_E_theta_z, an_E_theta)
        err_abs_theta_z = normalised_RMS(N.abs(surf_E_theta_z),
                                       N.abs(an_E_theta))
        err_x = normalised_RMS(
            surf_E_theta_x, an_E_theta, surf_E_phi_x)
        err_theta_x = normalised_RMS(surf_E_theta_x, an_E_theta)
        err_abs_theta_x = normalised_RMS(N.abs(surf_E_theta_x),
                                       N.abs(an_E_theta))
        err_y = normalised_RMS(
            surf_E_theta_y, an_E_theta, surf_E_phi_y)
        err_theta_y = normalised_RMS(surf_E_theta_y, an_E_theta)
        err_abs_theta_y = normalised_RMS(N.abs(surf_E_theta_y),
                                       N.abs(an_E_theta))
        print 'Far-field RMS error: ', err_z, err_x, err_y
        # Expected error for lam/6 mesh, 2nd order discretisation,
        # lam/4 current fillament source is ~4.685%
        self.assertTrue(err_z < 4.7)      
        self.assertTrue(err_x < 4.7)      
        self.assertTrue(err_z < 4.7)      
    


