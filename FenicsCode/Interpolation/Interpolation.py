# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as np
import dolfin

class SurfaceInterpolant(object):
    """Calculate a complex surface interpolant on a 3D domain"""
    def __init__(self, function_space):
        V = self.function_space = function_space
        self.u_r = dolfin.Function(V)
        self.u_i = dolfin.Function(V)
        
    def set_interpolant(self, interpolant):
        """Set interpolant f([x,y,z]) -> [f_x, f_y, f_z]"""
        self.interpolant = interpolant

        class int_r(dolfin.Expression):
            def __init__(self, interpolant):
                self.interpolant = interpolant
                
            def value_shape(self):
                return (3,)

            def eval(self, value, x):
                value[:] = self.interpolant(x).real

        class int_i(dolfin.Expression):
            def __init__(self, interpolant):
                self.interpolant = interpolant

            def value_shape(self):
                return (3,)

            def eval(self, value, x):
                value[:] = self.interpolant(x).imag

        self.interpolant_expression_Re = int_r(interpolant)
        self.interpolant_expression_Im = int_i(interpolant)
        
    def set_interpolant_expression(self, expr_r, expr_i):
        """Set interpolant as dolfin expression

        Can be used for faster evaluation if the interpolant is
        available as a dolfin expression. Two interpolants have to be
        passed, expr_r for the real part and expr_i for the imaginary
        part.
        """
        self.interpolant_expression_Re = expr_r
        self.interpolant_expression_Im = expr_i

    def calculate_interpolation(self):
        def boundary(x, on_boundary):
            return on_boundary
        V = self.function_space
        bc_Re = dolfin.DirichletBC(
            V, self.interpolant_expression_Re, boundary)
        bc_Im = dolfin.DirichletBC(
            V, self.interpolant_expression_Im, boundary)
        bc_Re.apply(self.u_r.vector())
        bc_Im.apply(self.u_i.vector())
        x = np.zeros(self.u_r.vector().size(), np.complex128)
        x[:] = self.u_r.vector().array() + 1j*self.u_i.vector().array()
        return x
