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
        
class test_fillament_source(unittest.TestCase):
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
        

