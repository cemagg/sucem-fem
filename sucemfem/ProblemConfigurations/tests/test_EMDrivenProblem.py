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

import unittest
import pickle
import numpy as N
import dolfin

from sucemfem.Testing import Meshes
from sucemfem.Testing import Paths
from sucemfem.Sources import point_source, current_source
from sucemfem.BoundaryConditions import ABCBoundaryCondition, BoundaryConditions
from sucemfem.ProblemConfigurations import EMDrivenProblem

# Module under test:

class test_DrivenProblemABC(unittest.TestCase):
    """Integration test for DrivenProblemABC class"""
    def setUp(self):
        self.source_coord = N.array([0,0,0.])
        self.source_value = N.array([2,-1,3.])
        self.frequency = 1e8
        testmesh = Meshes.InscribedTet()
        self.mesh = testmesh.get_dolfin_mesh()
        self.material_mesh_func = dolfin.MeshFunction('uint', self.mesh, 3)
        self.material_mesh_func.set_all(0)
        self.materials = {0:dict(eps_r=1, mu_r=1),}
        self.abc = ABCBoundaryCondition()
        self.abc.set_region_number(1)
        self.bcs = BoundaryConditions()
        self.bcs.add_boundary_condition(self.abc)
        self.current_sources = current_source.CurrentSources()
        self.dipole_source = point_source.PointCurrentSource()
        self.dipole_source.set_position(self.source_coord)
        self.dipole_source.set_value(self.source_value)
        self.current_sources.add_source(self.dipole_source)
        self.DUT = EMDrivenProblem.DrivenProblemABC()
        self.DUT.set_mesh(self.mesh)
        self.DUT.set_basis_order(1)
        self.DUT.set_material_regions(self.materials)
        self.DUT.set_region_meshfunction(self.material_mesh_func)
        self.DUT.set_boundary_conditions(self.bcs)
        self.DUT.set_sources(self.current_sources)
        
    def test_get_LHS_matrix(self):
        self.DUT.set_frequency(self.frequency)
        self.DUT.init_problem()
        # test data generated using a suspected-working version on 31 May 2011
        actual_LHSmat = self.DUT.get_LHS_matrix().todense()
        desired_file = Paths.get_module_path_filename('LHS_matrix.npy', __file__)
        desired_LHSmat = N.load(desired_file)
        self.assertTrue(N.allclose(
            actual_LHSmat, desired_LHSmat, rtol=1e-10, atol=3e-15))

    def test_get_RHS(self):
        self.DUT.set_frequency(self.frequency)
        self.DUT.init_problem()
        actual_RHS = self.DUT.get_RHS()
        # test data generated using a suspected-working version on 31 May 2011
        desired_file = Paths.get_module_path_file('RHS_vector.pickle', __file__)
        desired_RHS    = pickle.load(desired_file)
        self.assertTrue(N.allclose(
            actual_RHS, desired_RHS, rtol=1e-12, atol=1e-16))
        
