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
#    boundary_value_expression = dolfin.Expression(("0.0", "0.0", "0.0"), degree=1)

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

    def __init_boundary_value_expression ( self ):
        
        expr = ()
        for i in range(self.function_space.mesh().geometry().dim()):
            expr += ('0.0',)
        
        self.boundary_value_expression = dolfin.Expression(expr, degree=1)

    def get_essential_application_func(self, function_space=None):
        """Return an essential boundary condition application function.

        See parent class documentation for more details
        """
        if function_space is not None: self.set_function_space ( function_space )
                
        if self.boundary_value_expression is None:
            self.__init_boundary_value_expression()
            
        self._dirichletBC = dolfin.DirichletBC(self.function_space, 
                                               self.boundary_value_expression, 
                                               self.mesh_function, 
                                               self.region_number)
        
        return self._dirichletBC.apply
