# Authors:
# Neilen Marais <nmarais@gmail.com>
# Evan Lezar <mail@evanlezar.com>

import dolfin

from FenicsCode import Forms 
from FenicsCode import Materials 
from FenicsCode import SystemMatrices

import FenicsCode.BoundaryConditions.container

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
        self.boundary_conditions = FenicsCode.BoundaryConditions.container.BoundaryConditions()
    
    def get_global_dimension(self):
        """Return total number of system dofs, including Dirichlet constrained dofs
        """
        return self.function_space.dofmap().global_dimension()
    
    def set_mesh(self, mesh):
        self.mesh = mesh;
    
    def set_basis_order(self, order):
        self.basis_order = order
        if self.mesh is not None:
            self._init_function_space()
    
    def set_boundary_conditions(self, bcs):
        """
        Set the boundary conditions for the problem based on the keyword arguments passed
        or with the boundary condition object provided
        
        @param bcs: A BoundaryConditions object that stores the boundary conditions to be
            applied to the problem
        """
        self.boundary_conditions = bcs
    
    def set_material_regions(self, material_regions):
        """Set material region properties

        See documentation of Materials.MaterialPropertiesFactory for input format
        """
        self.material_regions = material_regions
    
    def set_region_meshfunction(self, region_meshfunction):
        self.region_meshfunction = region_meshfunction
        
    def _init_boundary_conditions(self):
        self.boundary_conditions.set_function_space(self.function_space)
    
    def _init_combined_forms (self):
        self.combined_forms = self.FormCombiner()
        self.combined_forms.set_interior_forms(self.interior_forms)
        self.combined_forms.set_boundary_conditions(self.boundary_conditions)
        
    def _init_function_space (self):
        if self.function_space is None:
            self.function_space = dolfin.FunctionSpace(
                self.mesh, self.element_type, self.basis_order)
    
    def _init_interior_forms(self):
        self.interior_forms = Forms.EMGalerkinInteriorForms()
        self.interior_forms.set_material_functions(self.material_functions)
        self.interior_forms.set_function_space(self.function_space)
            
    def _init_material_properties (self):
        mat_props_fac = Materials.MaterialPropertiesFactory(self.material_regions)
        mat_func_fac = Materials.MaterialFunctionFactory(
            mat_props_fac.get_material_properties(), 
            self.region_meshfunction, 
            self.mesh )
        self.material_functions = mat_func_fac.get_material_functions ( 'eps_r', 'mu_r' )
    
    def _init_system_matrices (self, matrix_class=None):
        bilin_forms = self.combined_forms.get_forms()
        sysmats = SystemMatrices.SystemMatrices()
        if matrix_class is not None:
            sysmats.set_matrix_class ( matrix_class )
        sysmats.set_matrix_forms(bilin_forms)
        sysmats.set_boundary_conditions(self.boundary_conditions)
        self.system_matrices = sysmats.calc_system_matrices()

    def init_problem(self):
        self._init_function_space()
        self._init_material_properties()
        self._init_interior_forms()
        self._init_boundary_conditions()
        self._init_combined_forms()
        self._init_system_matrices()
        
