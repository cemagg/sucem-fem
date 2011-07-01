# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin

from dolfin import inner, dot, dx, curl

class NullForm(object):
    """Null form that can be added to a dolfin form while making a zero contribution
    """
    _integrals = []

    def __add__(self, other):
        return other


class GalerkinInteriorForms(object):
    """Basic Galerkin System forms for the interior part of the problem

    The essential feature of a Galerkin system is that it uses the
    same test and trial spaces, barring the effects of boundary
    conditions.

    The interior forms are the basic system form without taking into
    account boundary conditions or excitations.  Only the standard
    'mass' (identity) and 'stiffness' (exterior derivative . exterior
    derivative) forms should be definied.
    """

    def set_material_functions(self, material_functions):
        self.material_functions = material_functions

    def set_function_space(self, function_space):
        """Set function space used to solve the problem"""
        V = self.function_space = function_space
        
        self.test_function = dolfin.TestFunction(V)
        self.trial_function = dolfin.TrialFunction(V)

    def get_mass_form(self):
        """Get 'mass'/identity matrix form"""
        raise NotImplementedError()
        
    def get_stiffness_form(self):
        """Get 'stiffness' / curl . curl matrix form"""
        raise NotImplementedError()

class EMGalerkinInteriorForms(GalerkinInteriorForms):
    """Weak Galerkin interior forms for EM vector wave equation

    The standard 'mass' (v.u) and 'stiffness' (curl(v).curl(u)) forms
    are defined.
    """
    def get_mass_form(self):
        """Get 'mass'/identity matrix form"""
        u = self.trial_function
        v = self.test_function
        eps_r = self.material_functions['eps_r']
        m = eps_r*inner(v, u)*dx
        return m
    
    def get_stiffness_form(self):
        """Get 'stiffness' / curl . curl matrix form"""
        u = self.trial_function
        v = self.test_function
        mu_r = self.material_functions['mu_r']
        s = dot(curl(v), curl(u))/mu_r*dx
        return s


class CombineGalerkinForms(object):
    """Base class for problem-specific form combination logic"""

    def set_interior_forms(self, interior_forms):
        """Set interior_forms with an instance of GalerkinInteriorForms"""
        self.interior_forms = interior_forms

    def set_boundary_conditions(self, boundary_conditions):
        """Set boundary_conditions with an instance of BoundaryConditions"""
        self.boundary_conditions = boundary_conditions

    def get_combined_forms(self):
        """Get combined forms

        Return value
        ------------

        combined_forms[matrix_name] -- a dict with forms as required
        to calculate system matrix 'matrix_name'. E.g.

        {'M': mass_form}

        where 'M' is the name of the matrix and 'mass_form' is the
        form required to calculate the matrix. 

        This function will raise a NotImplementedError. Concrete
        classes should sub-class CombineForms to provide an
        implementation. The concrete class should use a combination of
        self.interior_forms and self.boundary_conditions to calculate
        the final forms.
        """
        raise NotImplementedError("Use a concrete class")

