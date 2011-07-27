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

# Module under test
from sucemfem.PostProcessing import CalcEMFunctional

class test_CalcEMFunctional(unittest.TestCase):
    def setUp(self):
        self.mesh = dolfin.UnitCube(2,2,2)
        self.function_space = dolfin.FunctionSpace(
            self.mesh, "Nedelec 1st kind H(curl)", 1)
        nodofs = self.function_space.dofmap().global_dimension()
        self.E_dofs = np.random.random(nodofs) + 1j*np.random.random(nodofs)
        self.g_dofs = np.random.random(nodofs) + 1j*np.random.random(nodofs)
        self.k0 = np.random.rand()*10
        self.DUT = CalcEMFunctional(self.function_space)

    def test_symmetry(self):
        self.DUT.set_k0(self.k0)
        self.DUT.set_g_dofs(self.g_dofs)
        self.DUT.set_E_dofs(self.E_dofs)
        val1 = self.DUT.calc_functional()
        self.DUT.set_g_dofs(self.E_dofs)
        self.DUT.set_E_dofs(self.g_dofs)
        val2 = self.DUT.calc_functional()        
        self.assertEqual(val1, val2)
