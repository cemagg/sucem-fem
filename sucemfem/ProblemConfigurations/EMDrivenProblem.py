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

from __future__ import division

import numpy as N
import dolfin

from sucemfem import Forms 
from sucemfem import SystemMatrices
from sucemfem.Consts import c0, Z0
from sucemfem.Utilities.Converters import dolfin_ublassparse_to_scipy_csr
from EMProblem import EMProblem

class CombineForms(Forms.CombineGalerkinForms):
    def get_forms(self):
        m = self.interior_forms.get_mass_form()
        s = self.interior_forms.get_stiffness_form()
        ABC_form = self.boundary_conditions.get_bilinear_form()
        return dict(M=m, S=s, S_0=ABC_form)

class DrivenProblemABC(EMProblem):
    """Set up driven problem, potentially terminated by an ABC.
    
    Assumes lossless, frequency independent materials, and that the
    boundary bilinear form is:
        
    dot(cross(n, u), cross(n, v))
    
    where n is a face normal and u,v are the trial and testing
    functions. All forms are assumed to be real valued. They real
    forms will be combined into a complex system matrix of the form

    A = S - k0**2*M + 1j*k_0*S_0

    where S, M are stiffness and mass matrices and k0 is the
    freespace wave-number
    """

    FormCombiner = CombineForms

    def set_sources(self, sources):
        self.sources = sources

    def set_frequency(self, frequency):
        """Set simulation frequency in Hz"""
        self.frequency = frequency

    def get_LHS_matrix(self):
        k0 = 2*N.pi*self.frequency/c0
        M = dolfin_ublassparse_to_scipy_csr(self.system_matrices['M'])
        S = dolfin_ublassparse_to_scipy_csr(self.system_matrices['S'])
        S_0 = dolfin_ublassparse_to_scipy_csr(self.system_matrices['S_0'])
        return S - k0**2*M + 1j*k0*S_0

    def get_RHS(self):
        RHS = N.zeros(self.get_global_dimension(), N.complex128)
        dofnos, contribs = self._get_RHS_contributions()
        k0 = 2*N.pi*self.frequency/c0
        contribs = -1j*k0*Z0*contribs
        RHS[dofnos] += contribs
        return RHS

    def _init_system_matrices (self):
        """Initialise the system matrices associated with the problem. 
        Matrices are stored in dolfin.uBLASSparseMatrix format.
        """ 
        EMProblem._init_system_matrices(
            self, matrix_class=dolfin.uBLASSparseMatrix )

    def _get_RHS_contributions(self):
        self.sources.set_function_space(self.function_space)
        self.sources.init_sources()
        return self.sources.get_source_contributions()
