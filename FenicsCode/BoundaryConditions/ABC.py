# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin
from dolfin import inner, cross, ds

from FenicsCode.BoundaryConditions import BoundaryCondition

class ABCBoundaryCondition(BoundaryCondition):
    """
    Implements a first-order ABC, i.e. assuming 
    n x curl(E) + j*k_0*(n x (n x E)) = 0.

    Currently the ABC is applied to the whole exterior. Should be made
    to work with mesh function region numbers in the future.

    """
    def get_bilinear_form(self, test_function=None, trial_function=None):
        """Get bilinear form for implementing the ABC
        """
        V = self.function_space 
        if test_function is None:
            v = dolfin.TestFunction(self.function_space)
        else:
            v = test_function
        if trial_function is None:
            u = dolfin.TrialFunction(self.function_space)
        else:
            u = trial_function

        n = V.cell().n
        s_0 = inner(cross(n, v), cross(n, u))*ds
        return s_0
