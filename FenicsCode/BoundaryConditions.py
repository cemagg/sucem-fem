from __future__ import division

import dolfin
from FenicsCode import Forms

class BoundaryCondition(object):
    """Boundary condition base class

    A boundary condition can contribute to a system in three ways,
    viz. as an essential constraint on the trial space, a bilinear
    form contribution to the LHS and as a linear form contribution to
    the RHS.

    This base class describes a standard interface for obtaining the
    contributions of a basis function to those three aspects. It also
    provides nil-contribution default implementation.
    """

    def get_essential_application_func(self, function_space=None):
        """Return an essential boundary condition application function.

        Parameters
        ----------

        function_space -- Optional dolfin function space to use for
            constructing the essential boundary condition. If None is
            specified, the function space stored in self is used.

        Return value
        ------------

        apply(A, [b]) -- A function that applies the essential
            component of the boundary condition to the system matrix
            equation with the matrix A on LHS and the optional vector
            b on the RHS.
        """
        return lambda x: None

    def get_linear_form(self, test_function=None):
        """Return boundary condition's  linear form contribution as a dolfin form

        Parameters
        ----------

        test_function -- Optional dolfin testing function to use in
            the form construction. If None is specified, the testing
            function stored in self is used

        Return value
        ------------
        
        linear_form -- A dolfin linear form. The contribution of the
            linear form should be added to the form used to calculate
            the RHS of the system matrix equation that is eventually
            solved.
        """
        return Forms.NullForm()

    def get_bilinear_form(self, test_function=None, trial_function=None):
        """Return boundary condition's  bilinear form contribution as a dolfin form

        Parameters
        ----------

        test_function -- Optional dolfin testing function to use in
            the form construction. If None is specified, the testing
            function stored in self is used

        trial_function -- Optional dolfin trial function to use in
            the form construction. If None is specified, the trial
            function stored in self is used

        Return value
        ------------
        
        bilinear_form -- A dolfin linear form. This linear form should
            be added to the form used to calculate the RHS of the
            system matrix equation that is eventually solved.
        """
        return Forms.NullForm()

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

    def set_function_space(self, function_space):
        """Set function space on which the essential boundary condition is to be applied.

        Parameters
        ----------

        function_space -- dolfin function space object 

        """
        self.function_space = function_space

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
        
