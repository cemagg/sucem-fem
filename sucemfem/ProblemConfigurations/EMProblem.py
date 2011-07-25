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
# Evan Lezar <mail@evanlezar.com>

import dolfin

from sucemfem import Forms 
from sucemfem import Materials 
from sucemfem import SystemMatrices

from sucemfem.BoundaryConditions import BoundaryConditions

class EMProblem(object):
    """
    A base class for solving electromagnetic problems
    """
    def __init__ (self):
        self.element_type = "Nedelec 1st kind H(curl)"
        self.mesh = None
        self.order = None
        self.function_space = None
        self.material_regions = None
        self.region_meshfunction = None
        self.boundary_conditions = BoundaryConditions()
    
    def get_global_dimension(self):
        """Return total number of system dofs, including Dirichlet constrained dofs
        """
        return self.function_space.dofmap().global_dimension()
    
    def set_mesh(self, mesh):
        """Set the mesh associated with the problem.
        
        @param mesh: A dolfin mesh that defines the geometric discretisation of the problem.
        """
        self.mesh = mesh;
    
    def set_basis_order(self, order):
        """Set the order of the basis functions used in the discretisation of the problem.
        
        @param order: The polynomial order of the basis functions.
        """
        self.basis_order = order
        if self.mesh is not None:
            self._init_function_space()
    
    def set_boundary_conditions(self, bcs):
        """
        Set the boundary conditions for the problem based on the keyword arguments passed
        or with the boundary condition object provided
        
        @param bcs: A BoundaryConditions object containing a collection of boundary conditions, 
            or a single BoundaryCondition.
        """
        if type(bcs) == BoundaryConditions: 
            self.boundary_conditions = bcs
        else:
            self.boundary_conditions.add_boundary_condition ( bcs )
    
    def set_material_regions(self, material_regions):
        """Set material region properties
        
        See documentation of L{Materials.MaterialPropertiesFactory} for input format
        """
        self.material_regions = material_regions
    
    def set_region_meshfunction(self, region_meshfunction):
        self.region_meshfunction = region_meshfunction
        
    def _init_boundary_conditions(self):
        """Initialise the boundary conditions associated with the problem.
        """
        self.boundary_conditions.set_function_space(self.function_space)
    
    def _init_combined_forms (self):
        """Initialise the Dolfin forms for the problem.
        """
        self.combined_forms = self.FormCombiner()
        self.combined_forms.set_interior_forms(self.interior_forms)
        self.combined_forms.set_boundary_conditions(self.boundary_conditions)
        
    def _init_function_space (self):
        """If required, initialise a dolfin function space from the stored mesh, element_type, and basis function order information.
        """
        if self.function_space is None:
            self.function_space = dolfin.FunctionSpace(
                self.mesh, self.element_type, self.basis_order)
    
    def _init_interior_forms(self):
        """Initialise the Galerkin interior forms for the problem.
        """
        self.interior_forms = Forms.EMGalerkinInteriorForms()
        self.interior_forms.set_material_functions(self.material_functions)
        self.interior_forms.set_function_space(self.function_space)
            
    def _init_material_properties (self):
        mat_props_fac = Materials.MaterialPropertiesFactory ( self.material_regions )
        mat_func_fac = Materials.MaterialFunctionFactory(
            mat_props_fac.get_material_properties(), 
            self.region_meshfunction, 
            self.mesh )
        self.material_functions = mat_func_fac.get_material_functions ( 'eps_r', 'mu_r' )
    
    def _init_system_matrices (self, matrix_class=None):
        """Initialise the system matrices associated with the problem.
        
        @keyword matrix_class: An optional dolfin class to use for matrix storage.
            (default: None).
        """
        bilin_forms = self.combined_forms.get_forms()
        sysmats = SystemMatrices.SystemMatrices()
        if matrix_class is not None:
            sysmats.set_matrix_class ( matrix_class )
        sysmats.set_matrix_forms(bilin_forms)
        sysmats.set_boundary_conditions(self.boundary_conditions)
        self.system_matrices = sysmats.calc_system_matrices()

    def init_problem(self):
        """Perform the final initialisation of the problem components.
        """
        self._init_function_space()
        self._init_material_properties()
        self._init_interior_forms()
        self._init_boundary_conditions()
        self._init_combined_forms()
        self._init_system_matrices()
        
