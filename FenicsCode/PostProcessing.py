import dolfin
import numpy as N

class Reconstruct(object):
    """Reconstruct field values, dealing with complex numbers as required"""

    def __init__(self, function_space):
        """Initialise reconstruction for dolfin FunctionSpace object function_space"""
        self.function_space = function_space

    def set_dof_values(self, x):
        """Set dof values x as numpy.array"""
        self.dof_values = x

    def reconstruct_points(self, points):
        """Reconstruct function at points

        Where points is an n x d array for a d-dimensional function
        """
        u = dolfin.Function(self.function_space)
        # N.require() to make sure that the array is continuous.
        u.vector()[:] = N.require(N.real(self.dof_values), requirements='C')
        finite_element = self.function_space.dolfin_element()
        bf_value_dimension = finite_element.value_dimension(0)
        values_re = N.zeros((len(points), bf_value_dimension), dtype=N.float64)
        values = values_re
        for i, fp in enumerate(points):
            values_re[i,:] = u(fp) 
        # Add imaginary part, if any
        if N.any(N.iscomplex(self.dof_values)):
            values_im = N.zeros_like(values_re)
            for i, fp in enumerate(points):
                u.vector()[:] = N.require(N.imag(self.dof_values), requirements='C')
                values_im[i,:] = u(fp)
            values = values_re + 1j*values_im
                
        return values
