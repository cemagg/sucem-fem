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
import numpy as np
import dolfin

# Module under test:
from sucemfem.Sources.PostProcess import VoltageAlongLine
from sucemfem.Sources.PostProcess import ComplexVoltageAlongLine

class test_VoltageAlongLine(unittest.TestCase):
    def setUp(self):
        self.mesh = dolfin.UnitCube(3,3,3)
        self.V = dolfin.FunctionSpace(self.mesh, "Nedelec 1st kind H(curl)", 2)
        self.u = dolfin.interpolate(
            dolfin.Expression(('0','0', '2*x[2]')), self.V)
        self.DUT = VoltageAlongLine(self.u)

    def test_calculate_voltage(self):
        # Should result in 1v
        p1v_pts = (np.array([[0.5,0.35,0], [0.5,0.5,1]]),
                   np.array([[0,0,0], [1,1,1]]))
        # Should result in -1v
        m1v_pts = (np.array([[0.25,0.5,1-3e-16], [0.5,0.75,3e-16]]),
                   np.array([[1,1,1], [0,0,0]]))
        # Should result in 0v
        p0v_pts = np.array([[0,0,1], [1,1,1]])
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p1v_pts[0]), 1.)
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p1v_pts[1]), 1.)
        self.assertAlmostEqual(self.DUT.calculate_voltage(*m1v_pts[0]), -1.)
        self.assertAlmostEqual(self.DUT.calculate_voltage(*m1v_pts[1]), -1.)
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p0v_pts), 0.)


class test_ComplexVoltageAlongLine(unittest.TestCase):
    def setUp(self):
        self.mesh = dolfin.UnitCube(3,3,3)
        self.V = dolfin.FunctionSpace(self.mesh, "Nedelec 1st kind H(curl)", 3)
        self.u_r = dolfin.interpolate(
            dolfin.Expression(('0','0', '2*x[2]')), self.V)
        self.u_i = dolfin.interpolate(
            dolfin.Expression(('0','0', '-x[2]*x[2]')), self.V)
        self.x = self.u_r.vector().array() + 1j*self.u_i.vector().array()
        self.DUT = ComplexVoltageAlongLine(self.V)
        self.DUT.set_dofs(self.x)

    def test_calculate_voltage(self):
        # Should result in (1. - 1j/3) V
        p1v_pts = (np.array([[0.5,0.35,0], [0.5,0.5,1]]),
                   np.array([[0,0,0], [1,1,1]]))
        # Should result in (-1. + 1j/3) v
        m1v_pts = (np.array([[0.25,0.5,1-3e-16], [0.5,0.75,3e-16]]),
                   np.array([[1,1,1], [0,0,0]]))
        # Should result in 0v
        p0v_pts = np.array([[0,0,1], [1,1,1]])
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p1v_pts[0]), (1 - 1j/3))
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p1v_pts[1]), (1 - 1j/3))
        self.assertAlmostEqual(self.DUT.calculate_voltage(*m1v_pts[0]), -(1 - 1j/3))
        self.assertAlmostEqual(self.DUT.calculate_voltage(*m1v_pts[1]), -(1 - 1j/3))
        self.assertAlmostEqual(self.DUT.calculate_voltage(*p0v_pts), 0.)


