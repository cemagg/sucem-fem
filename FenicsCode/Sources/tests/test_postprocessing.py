# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import unittest
import numpy as np
import dolfin

# Module under test:
from FenicsCode.Sources.PostProcess import VoltageAlongLine

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

    
