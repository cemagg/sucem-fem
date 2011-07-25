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
import numpy as np
import dolfin
from FenicsCode.Testing import Paths
from FenicsCode.Utilities.MeshGenerators import get_centred_cube

# Module under test:
from FenicsCode.Sources.fillament_current_source import FillamentCurrentSource
from FenicsCode.Sources.fillament_source import FillamentSource
        
class FillamentEnvironment(object):
    def __init__(self, datafile):
        test_data = pickle.load(datafile)
        self.desired_dofnos = test_data['dofnos']
        self.desired_rhs_contribs = test_data['rhs_contribs']
        self.discretisation_order = test_data['order']
        self.mesh = get_centred_cube(
            test_data['domain_size'], test_data['max_edge_len'])
        self.discretisation_space = dolfin.FunctionSpace(
            self.mesh, "Nedelec 1st kind H(curl)", self.discretisation_order)
        self.no_integration_points = test_data['no_integration_points']
        self.I = test_data['I']
        self.source_endpoints = test_data['source_endpoints']
        
class test_fillament_current_source(unittest.TestCase):
    test_data_file = 'fillament_source_test_data.pickle'
    rtol=1e-10
    atol=1e-12
    
    def setUp(self):
        desired_file = Paths.get_module_path_file(self.test_data_file, __file__)
        env = self.environment = FillamentEnvironment(desired_file)
        self.DUT = FillamentCurrentSource()
        self.DUT.set_function_space(env.discretisation_space)
        self.DUT.set_no_integration_points(env.no_integration_points)
        self.DUT.set_source_endpoints(env.source_endpoints)
        self.DUT.set_value(env.I)

    def test_source(self):
        actual_dofnos, actual_rhs = self.DUT.get_contribution()
        sortind = np.argsort(actual_dofnos)
        actual_dofnos = actual_dofnos[sortind].copy()
        actual_rhs = actual_rhs[sortind].copy()
        self.assertTrue(np.all(actual_dofnos == self.environment.desired_dofnos))
        self.assertTrue(np.allclose(actual_rhs, self.environment.desired_rhs_contribs,
                                    rtol=self.rtol, atol=self.atol))

    def test_negative_direction(self):
        dofnos, rhs = self.DUT.get_contribution()
        sortind = np.argsort(dofnos)
        dofnos = dofnos[sortind]
        rhs = rhs[sortind]
        self.DUT.set_source_endpoints(-self.DUT.source_endpoints)
        dofnos_neg, rhs_neg = self.DUT.get_contribution()
        sortind = np.argsort(dofnos_neg)
        dofnos_neg = dofnos_neg[sortind]
        rhs_neg = rhs_neg[sortind]
        self.assertTrue(np.all(dofnos == dofnos_neg))
        self.assertTrue(np.allclose(rhs, -rhs_neg))
        
class test_fillament_source(unittest.TestCase):
    def setUp(self):
        self.mesh = dolfin.UnitCube(1,1,1)
        self.V = dolfin.FunctionSpace(self.mesh, "Nedelec 1st kind H(curl)", 1)
        self.DUT = FillamentSource(self.V)
        self.fillament_current = 2
        self.fillament_endpoints = np.array([[0,0,0.7], [0,0,0.2]])
        self.source_parameters = dict(
            I=self.fillament_current, endpoints=self.fillament_endpoints)

    def test_set_source_parameters(self):
        self.DUT.set_source_parameters(self.source_parameters)
        self.assertTrue(np.all(self.DUT.direction == [0,0,-1]))
        self.assertAlmostEqual(self.DUT.length, 0.5)

    def test_get_current_source(self):
        self.DUT.set_source_parameters(self.source_parameters)
        cs = self.DUT.get_current_source()
        self.assertTrue(isinstance(cs, FillamentCurrentSource))
        self.assertEqual(cs.value, self.fillament_current)
        self.assertTrue(np.all(cs.source_endpoints == self.fillament_endpoints))
