
from __future__ import division

import numpy as N
import dolfin

from FenicsCode.ProblemConfigurations import EMVectorWaveEigen
from FenicsCode import Forms 
from FenicsCode import SystemMatrices
from FenicsCode.Consts import c0, Z0
from FenicsCode.Utilities.Converters import dolfin_ublassparse_to_scipy_csr

class CombineForms(Forms.CombineGalerkinForms):
    def get_forms(self):
        m = self.interior_forms.get_mass_form()
        s = self.interior_forms.get_stiffness_form()
        ABC_form = self.boundary_conditions.get_bilinear_form()
        return dict(M=m, S=s, S_0=ABC_form)

class DrivenProblemABC(EMVectorWaveEigen.EigenProblem):
    """Set up driven problem, potentially terminated by an ABC

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

    CombineForms = CombineForms
    
    def set_boundary_conditions(self, bcs):
        """Set boundary conditions with a BoundaryConditions object"""
        self.boundary_conditions = bcs

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

    def get_global_dimension(self):
        """Return total number of system dofs, including Dirichlet constrained dofs
        """
        return self.function_space.dofmap().global_dimension()

    def get_RHS(self):
        RHS = N.zeros(self.get_global_dimension(), N.complex128)
        dofnos, contribs = self._get_RHS_contributions()
        k0 = 2*N.pi*self.frequency/c0
        contribs = 1j*k0*Z0*contribs
        RHS[dofnos] += contribs
        return RHS

    def _get_boundary_conditions(self):
        try:
            bcs = self.boundary_conditions
        except AttributeError:
            raise AttributeError('set_boundary_conditions() method must be called first')
        bcs.set_function_space(self.function_space)
        return bcs

    def _get_system_matrices(self):
        bilin_forms = self.combined_forms.get_forms()
        sysmats = SystemMatrices.SystemMatrices()
        sysmats.set_matrix_class(dolfin.uBLASSparseMatrix)
        sysmats.set_matrix_forms(bilin_forms)
        sysmats.set_boundary_conditions(self.boundary_conditions)
        return sysmats.calc_system_matrices()

    def _get_RHS_contributions(self):
        self.sources.set_function_space(self.function_space)
        self.sources.init_sources()
        return self.sources.get_source_contributions()
