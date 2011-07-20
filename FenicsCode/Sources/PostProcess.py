from __future__ import division

import dolfin
import numpy as np
from scipy.integrate import romberg
from FenicsCode.Utilities.Geometry import unit_vector, vector_length
from FenicsCode.Utilities.Converters import as_dolfin_vector

class VoltageAlongLine(object):
    """Measure voltage along a straight line between two points"""
    def __init__(self, field_function):
        self.field_function = field_function

    def calculate_voltage(self, start_pt, end_pt):
        fn = self.field_function
        delta = end_pt - start_pt
        l_hat = unit_vector(delta)
        # Evaluate E . l_hat where E is the electric field vector and
        # l_hat is a unit vector along the integration path.
        eval_fn = lambda l: np.dot(l_hat, fn(*(start_pt + delta*l)))
        # Integrate over unit length
        intg = romberg(eval_fn, 0, 1)
        # Multiply by inteval length to de-normalise
        interval_len = vector_length(delta)
        intg = intg*interval_len
        return intg

class ComplexVoltageAlongLine(object):
    """Measure complex voltage along a straight line between two points"""
    def __init__(self, function_space):
        self.function_space = function_space

    def set_dofs(self, dofs):
        self.dofs = dofs
        self.x_r = as_dolfin_vector(self.dofs.real)
        self.x_i = as_dolfin_vector(self.dofs.imag)
        self.E_r = dolfin.Function(self.function_space, self.x_r)
        self.E_i = dolfin.Function(self.function_space, self.x_i)
        self.real_voltage = VoltageAlongLine(self.E_r)
        self.imag_voltage = VoltageAlongLine(self.E_i)

    def calculate_voltage(self, start_pt, end_pt):
        return (self.real_voltage.calculate_voltage(start_pt, end_pt) +
                1j*self.imag_voltage.calculate_voltage(start_pt, end_pt))
    

        
        
