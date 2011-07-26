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
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin
import numpy as N
from dolfin import dx, dot, curl

__all__ = ['Reconstruct', 'CalcEMFunctional']

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
        self.epsr_function = None
        self.mur_function = None
        self._form_compiler_parameters = None
        self.dirty = True
        
    def set_k0(self, k0):
        self.k0 = k0
        self.dirty = True


    def set_cell_domains(self, cell_domains, cell_mark_value):
        """Set cell domain mesh function

        Needed if the functional is to be calculated over only part of the domain
        """
        self.cell_domains = cell_domains
        self.cell_mark_value = cell_mark_value
        self.dx = dx(self.cell_mark_value)
        self.dirty = True
        
    def set_epsr_function(self, permittivity_function):
        self.epsr_function = epsr_function
        self.dirty = True
        
    def set_mur_function(self, mur_function):
        self.mur_function = mur_function
        self.dirty = True

    def set_quadrature_degree(self, quadrature_degree):
        """Optionally set quadrature degre, otherwise dolfin auto is used"""
        self.quadrature_degree = quadrature_degree
        if self._form_compiler_parameters is None:
            self._form_compiler_parameters = {}
        self._form_compiler_parameters.update(dict(
            quadrature_degree=quadrature_degree))

    def _get_epsr_function(self):
        if self.epsr_function is not None:
            epsr_func = self.epsr_function
        else:
            epsr_func = dolfin.Constant(1)

        return epsr_func
   
    def _get_mur_function(self):
        if self.epsr_function is not None:
            mur_func = self.mur_function
        else:
            mur_func = dolfin.Constant(1)

        return mur_func

    def _get_forms(self):
        if self.dirty:
            E_r, E_i, g_r, g_i = self.E_r, self.E_i, self.g_r, self.g_i
            k0 = self.k0
            eps_r = self._get_epsr_function()
            mu_r = self._get_mur_function()
            form_r = (dot(curl(E_r)/mu_r, curl(g_r)) - dot(curl(E_i)/mu_r, curl(g_i)) \
                      - k0**2*dot(eps_r*E_r, g_r) - dot(eps_r*E_i, g_i))*self.dx
            form_i = (dot(curl(E_r)/mu_r, curl(g_i)) + dot(curl(E_i)/mu_r, curl(g_r)) \
                  - k0**2*dot(eps_r*E_r, g_i) + dot(eps_r*E_i, g_r))*self.dx
            self.form_r, self.form_i = form_r, form_i
            self.dirty = False

        return self.form_r, self.form_i

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
        form_r, form_i = self._get_forms()
        I_r = dolfin.assemble(
            self.form_r, cell_domains=self.cell_domains,
            form_compiler_parameters=self._form_compiler_parameters)
        I_i = dolfin.assemble(
            self.form_i, cell_domains=self.cell_domains,
            form_compiler_parameters=self._form_compiler_parameters)

        return I_r + 1j*I_i

