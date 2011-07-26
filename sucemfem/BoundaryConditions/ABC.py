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
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
import dolfin
from dolfin import inner, cross, ds

from sucemfem.BoundaryConditions import BoundaryCondition

class ABCBoundaryCondition(BoundaryCondition):
    """
    Implements a first-order ABC, i.e. assuming 
    
    M{n S{times} curl(E) + j*k_0*(n S{times} (n S{times} E)) = 0.}
    
    @todo: Currently the ABC is applied to the whole exterior. Should be made
    to work with mesh function region numbers in the future.
    """
    def get_bilinear_form(self, test_function=None, trial_function=None):
        """Calculate and return the bilinear form associated with the boundary condition.
        
        The form is calculated as: M{<n S{times} v, n S{times} u>}, where M{v} and M{u} are the 
        test and trial functions, respectively.
        
        @keyword test_function: override the stored test function space. 
            (default: None.)
        @keyword trial_function: override the stored trial function space.
            (default: None.)
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
