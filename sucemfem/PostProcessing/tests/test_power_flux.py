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

import unittest
import os
import pickle
import numpy as np
import dolfin
from sucemfem.Testing import Paths
from sucemfem.Consts import c0

# Module under test
from sucemfem.PostProcessing import power_flux


class test_SurfaceFlux(unittest.TestCase):
    data_dir = 'data'
    test_data_file = 'power_flux_reference_dofs_f-1000000000.000000_o-2_s-0.074948_l-0.100000_h-0.166667.pickle'
    
    def setUp(self):
        data_dir = Paths.get_module_path_filename(self.data_dir,  __file__)
        data_file = open(os.path.join(data_dir, self.test_data_file))
        data = self.data = pickle.load(data_file)
        self.mesh = dolfin.Mesh(os.path.join(data_dir, data['meshfile']))
        self.discretisation_order = data['order']
        self.function_space = dolfin.FunctionSpace(
            self.mesh, "Nedelec 1st kind H(curl)", self.discretisation_order)
        self.k0 = 2*np.pi*data['freq']/c0
        self.DUT = power_flux.SurfaceFlux(self.function_space)

    def test_calc_flux(self):
        self.DUT.set_k0(self.k0)
        self.DUT.set_dofs(self.data['x'])
        # Desired result as calculated with gitrev 953c7063b02547f5233a29ced884ba5af0fd0fe3
        self.assertAlmostEqual(self.DUT.calc_flux(), 26.3586862457)



class test_SurfaceFlux(unittest.TestCase):
    data_dir = 'data'
    test_data_file = 'power_flux_reference_dofs_f-1000000000.000000_o-2_s-0.074948_l-0.100000_h-0.166667.pickle'
    
    def setUp(self):
        data_dir = Paths.get_module_path_filename(self.data_dir,  __file__)
        data_file = open(os.path.join(data_dir, self.test_data_file))
        data = self.data = pickle.load(data_file)
        self.mesh = dolfin.Mesh(os.path.join(data_dir, data['meshfile']))
        self.discretisation_order = data['order']
        self.function_space = dolfin.FunctionSpace(
            self.mesh, "Nedelec 1st kind H(curl)", self.discretisation_order)
        self.k0 = 2*np.pi*data['freq']/c0
        self.DUT = power_flux.VariationalSurfaceFlux(self.function_space)

    def test_calc_flux(self):
        self.DUT.set_k0(self.k0)
        self.DUT.set_dofs(self.data['x'])
        # Desired result as calculated with gitrev 953c7063b02547f5233a29ced884ba5af0fd0fe3
        self.assertAlmostEqual(self.DUT.calc_flux(), (26.9015957282-7.85747240242e-08j))
        
