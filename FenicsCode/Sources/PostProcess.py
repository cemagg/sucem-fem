from __future__ import division

import numpy as np
from scipy.integrate import romberg
from FenicsCode.Utilities.Geometry import unit_vector, vector_length

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


        
