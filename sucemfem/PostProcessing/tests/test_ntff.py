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
import pickle
import numpy as N
import dolfin
from sucemfem.Testing import Paths
from sucemfem.Utilities.MeshGenerators import get_centred_cube

# Module under test:
from sucemfem.PostProcessing import surface_ntff
from sucemfem.PostProcessing import variational_ntff

class test_interpolant(unittest.TestCase):
    test_data_file = 'data/interpolant_test_data.pickle'
    rtol=1e-10
    atol=1e-7

    def setUp(self):
        desired_file = Paths.get_module_path_file(self.test_data_file, __file__)
        self.desired_data = pickle.load(desired_file)
            
class test_interpolant_expression(test_interpolant):
    def setUp(self):
        super(test_interpolant_expression, self).setUp()
        self.DUT = variational_ntff.TransformTestingExpression()
        
    def test_expression_re(self):
        dd = self.desired_data
        k0, ahats, rhat = dd['k0'], dd['ahats'], dd['rhat']
        ahats, coords = dd['ahats'], dd['coords']
        for i, ahat in enumerate(ahats):
            desired_vals = dd['vals'][i].real
            self.DUT.set_parms(rhat, ahat, k0)
            expr_r = self.DUT.get_expression()[0]
            actual_vals = N.zeros((len(coords),3), N.float64)
            for j, coord in enumerate(coords):
                expr_r.eval(actual_vals[j,:], coord)
            self.assertTrue(N.allclose(actual_vals, desired_vals,
                                       rtol=self.rtol, atol=self.atol))
        
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
    test_data_file = 'data/reference_surface_ntff-2-0.149896229-0.0499654096667.pickle'
    rtol=1e-12
    atol=1e-7
    
    def setUp(self):
        desired_file = Paths.get_module_path_file(self.test_data_file, __file__)
        self.environment = NTFFEnvironment(desired_file)
        self.DUT = surface_ntff.NTFF(self.environment.discretisation_space)

    def test_ff(self):
        env = self.environment
        self.DUT.set_frequency(env.frequency)
        self.DUT.set_dofs(env.discretisation_dofs)
        actual_E_ff = [self.DUT.calc_pt(th_deg, ph_deg)
                       for th_deg, ph_deg in zip(env.theta_coords, env.phi_coords)]
        self.assertTrue(N.allclose(actual_E_ff, env.desired_E_ff,
                                   rtol=self.rtol, atol=self.atol))
        

class test_variational_ntff(test_surface_ntff):
    test_data_file = 'data/reference_variational_ntff-2-0.149896229-0.0499654096667.pickle'
    def setUp(self):
        desired_file = Paths.get_module_path_file(self.test_data_file, __file__)
        self.environment = NTFFEnvironment(desired_file)
        self.DUT = variational_ntff.NTFF(self.environment.discretisation_space)
    
