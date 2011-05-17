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
