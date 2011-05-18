import dolfin
from dolfin import inner, cross, ds

from FenicsCode.BoundaryConditions import BoundaryCondition

class ABCBoundaryCondition(BoundaryCondition):
    def get_bilinear_form(self, test_function=None, trial_function=None):
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
