from __future__ import division

import unittest
import pickle
import numpy as N
import dolfin
from FenicsCode.Testing import Paths
from FenicsCode.Utilities.MeshGenerators import get_centred_cube

# Module under test:
import ntff


class NTFFEnvironment(object):
    def __init__(self, datafile):
        test_data = pickle.load(datafile)
        ff_result_data = test_data['ff_result_data']
        nf_input_data = test_data['nf_input_data']
        self.desired_E_ff = ff_result_data['E_ff']
        self.desired_E_theta = self.desired_E_ff[:,0]
        self.desired_E_phi = self.desired_E_ff[:,1]
        self.theta_coords = ff_result_data['theta_deg']
        self.phi_coords = ff_result_data['phi_deg']
        self.discretisation_order = nf_input_data['order']
        self.mesh = get_centred_cube(
            nf_input_data['domain_size'], nf_input_data['max_edge_len'])
        self.discretisation_space = dolfin.FunctionSpace(
            self.mesh, "Nedelec 1st kind H(curl)", self.discretisation_order)
        self.discretisation_dofs = nf_input_data['x']
        self.frequency = nf_input_data['freq']
        
class test_surface_ntff(unittest.TestCase):
    test_data_file = 'reference_surface_ntff-2-0.149896229-0.0499654096667.pickle'
    rtol=1e-12
    atol=1e-16
    
    def setUp(self):
        desired_file = Paths.get_module_path_file(self.test_data_file, __file__)
        self.environment = NTFFEnvironment(desired_file)
        self.DUT = ntff.NTFF(self.environment.discretisation_space)

    def test_ff(self):
        env = self.environment
        self.DUT.set_frequency(env.frequency)
        self.DUT.set_dofs(env.discretisation_dofs)
        self.DUT.init_calc()
        actual_E_ff = [self.DUT.calc_pt(th_deg, ph_deg)
                       for th_deg, ph_deg in zip(env.theta_coords, env.phi_coords)]
        self.assertTrue(N.allclose(actual_E_ff, env.desired_E_ff,
                                   rtol=self.rtol, atol=self.atol))
        
