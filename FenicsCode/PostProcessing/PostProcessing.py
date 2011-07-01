# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin
import numpy as N
from dolfin import dx, dot, curl

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

class CalcEMFunctional(object):
    """Evaluate EM functional, assuming freespace

    Evaluate <curl(E), curl(g)> - k0^2 <E, g>
    where E is the E-field solution and g is an H(curl) testing function
    
    """

    def __init__(self, function_space, testing_space=None):
        V = self.function_space = function_space
        if testing_space is None:
            self.testing_space = V
        else:
            self.testing_space = testing_space
        Vt = self.testing_space
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)
        self.g_r = dolfin.Function(Vt)
        self.g_i = dolfin.Function(Vt)
        self.dx = dx
        self.cell_domains = None
        
    def set_k0(self, k0):
        self.k0 = k0
        self.form_r, self.form_i = self._get_forms()

    def set_cell_domains(self, cell_domains, cell_mark_value):
        """Set cell domain mesh function

        Needed if the functional is to be calculated over only part of the domain
        """
        self.cell_domains = cell_domains
        self.cell_mark_value = cell_mark_value
        self.dx = dx(self.cell_mark_value)
        
    def _get_forms(self):
        E_r, E_i, g_r, g_i = self.E_r, self.E_i, self.g_r, self.g_i
        k0 = self.k0
        form_r = (dot(curl(E_r), curl(g_r)) - dot(curl(E_i), curl(g_i)) \
                  + k0**2*dot(E_r, g_r) - dot(E_i, g_i))*self.dx
        form_i = (dot(curl(E_r), curl(g_i)) + dot(curl(E_i), curl(g_r)) \
                  + k0**2*dot(E_r, g_i) - dot(E_i, g_r))*self.dx
        return form_r, form_i

    def set_E_dofs(self, E_dofs):
        x_r = N.real(E_dofs).copy()
        x_i = N.imag(E_dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i
        
    def set_g_dofs(self, g_dofs):
        x_r = N.real(g_dofs).copy()
        x_i = N.imag(g_dofs).copy()
        self.g_r.vector()[:] = x_r
        self.g_i.vector()[:] = x_i

    def calc_functional(self):
        """Calculate functional using given E, g and k0 values"""
        I_r = dolfin.assemble(self.form_r, cell_domains=self.cell_domains)
        I_i = dolfin.assemble(self.form_i, cell_domains=self.cell_domains)

        return I_r + 1j*I_i
