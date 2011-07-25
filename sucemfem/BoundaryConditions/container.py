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
from __future__ import division
from sucemfem import Forms


class BoundaryConditions(object):
    def __init__(self):
        self.boundary_conditions = {}
        
    def add_boundary_condition(self, boundary_condition):
        bc_num = boundary_condition.region_number
        # boundary condition numbers have to be unique
        assert(bc_num not in self.boundary_conditions)
        self.boundary_conditions[bc_num] = boundary_condition

    def apply_essential(self, A, b=None):
        """
        Apply essential boundary conditions to system matrix A and optional RHS b
        """
        for bc_num, bc in self.boundary_conditions.items():
            apply_fn = bc.get_essential_application_func()
            if b: apply_fn(A, b)
            else: apply_fn(A)

    def get_linear_form(self):
        """Get boundary conditions contribution to RHS linear form
        """
        lin_form = Forms.NullForm()
        for bc_num, bc in self.boundary_conditions.items():
            lin_form = lin_form + bc.get_linear_form()

        return lin_form

    def get_bilinear_form(self):
        """Get boundary conditions contribution to bilinear form
        """

        bilin_form = Forms.NullForm()
        for bc_num, bc in self.boundary_conditions.items():
            bilin_form = bilin_form + bc.get_bilinear_form()

        return bilin_form

    def set_function_space(self, function_space):
        """Set all the boundary conditions in the collection's function_space
        """
        self.function_space = function_space
        for bc in self.boundary_conditions.values():
            bc.set_function_space(function_space)

    
