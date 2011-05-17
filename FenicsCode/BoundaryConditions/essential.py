from __future__ import division
import dolfin
from FenicsCode.BoundaryConditions import BoundaryCondition


class EssentialBoundaryCondition(BoundaryCondition):
    """Essential boundary condition class

    Sets up Dirichlet type boundary conditions for prescribed boundary
    values. The boundary conditions are applied to the user-specified
    function space self.function_space. If a non-zero BC is required,
    set_boundary_value_expression() needs to be called.

    The boundary region can be specified using a pre-defined mesh
    function (init_with_meshfunction()) or using a SubDomain object
    (init_with_subdomain())

    Also see documentation of the parent class
    """
    
    # Default boundary condition is zero, i.e. PEC
    boundary_value_expression = dolfin.Expression(("0.0", "0.0", "0.0"), degree=1)

    def set_boundary_value_expression(self, boundary_value_expression):
        """Set value expression for boundary condition

        The EssentialBoundaryCondition class defaults to zero
        (i.e. PEC) boundary condition if
        set_boundary_value_expression() is not called

        Parameters
        ----------

        boundary_value_expression -- dolfin expression object

        """
        self.boundary_value_expression = boundary_value_expression

    def init_with_meshfunction(self, mesh_function, region_number):
        """Initialise using a mesh function to indicate boundary region

        Parameters
        ----------

        mesh_function -- A 'uint' mesh function defined on the facets
            (i.e. of dimension one lower than mesh geometry) of the
            computational mesh

        region_number -- Region number (as defined in mesh_function)
            that this boundary condition should be applied to
        """
        self.region_number = region_number
        self.mesh_function = mesh_function
        self.subdomain = None

    def init_with_subdomain(self, subdomain, region_number):
        """Initialise using a mesh function to indicate boundary region

        Parameters
        ----------

        subdomain -- A SubDomain subclass defining the region of
            application for the boundary condition. The user must
            ensure that it does not overlap any other boundary
            condition region.

        region_number -- A unique region number for this boundary
            condition. The user must ensure that no other boundary
            condition uses the same region number.
            
        """
        self.region_number = region_number
        self.subdomain = subdomain
        self.mesh_function = None

    def get_essential_application_func(self, function_space=None):
        """Return an essential boundary condition application function.

        See parent class documentation for more details
        """
        if function_space: V = function_space
        else: V = self.function_space
        u_bdry = self.boundary_value_expression
        
        if self.mesh_function:
            self._dirichletBC = dolfin.DirichletBC(
                V, u_bdry, self.mesh_function, self.region_number)
        elif self.subdomain:
            self._dirichletBC = dolfin.DirichletBC(V, u_bdry, self.subdomain)
        else: raise Exception("Neither mesh function nor subdomain specified!")

        return self._dirichletBC.apply
