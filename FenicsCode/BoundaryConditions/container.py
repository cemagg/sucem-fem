# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division
from FenicsCode import Forms


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

    
