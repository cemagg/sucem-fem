import dolfin
from FenicsCode import Forms 
from FenicsCode import Materials 
from FenicsCode import BoundaryConditions 
from FenicsCode import SystemMatrices

from scipy.sparse.linalg.eigen.arpack import speigs
import FenicsCode.BoundaryConditions.essential
import FenicsCode.BoundaryConditions.container
from EMProblem import EMProblem

class CombineForms(Forms.CombineGalerkinForms):
    def get_forms(self):
        m = self.interior_forms.get_mass_form()
        s = self.interior_forms.get_stiffness_form()
        # I suspect that using e.g. the bilinear part of a first order
        # ABC in an eigen solution could potentially be used for eigen
        # solution, so I'm including the bilinear boundary condition
        # forms here
        BC_bilin_forms = self.boundary_conditions.get_bilinear_form()
        return dict(M=m, S=s, BC=BC_bilin_forms)

class EigenProblem(EMProblem):
    FormCombiner = CombineForms        

class DefaultEigenSolver(object):
    def set_eigenproblem(self, eigenproblem):
        """Sets initialised instance of EigenProblem to solve"""
        self.eigenproblem = eigenproblem

    def set_sigma(self, sigma):
        """Spectrum shift (sigma) to apply to k^2"""
        self.sigma = sigma

    def solve_problem(self, nev, ncv=None):
        """Solve problem for eigenvalues and eigenvectors
        Input Values
        -------------
        @param nev: Number of eigenpairs to compute
        @param ncv: Number of Arnoldi basisvectors to use. If None, default to 2*nev+1
        Return Values
        -------------
        (eigs_w, eigs_v) with

        eigs_w -- array of n_eig eigen values corresponding to k^2 for the
            mode's resonant wavenumber

        eigs_v -- 2D array ofshape n_eig x problem_dim eigen vectors
            corresponding to the modal distributions of the eigenvalues.
            eigs_v[i] is the eigenvector corresponding to eigs_w[i]
        
        """
        M = self.eigenproblem.system_matrices['M']
        S = self.eigenproblem.system_matrices['S']
        solve_mat = S - self.sigma*M
        lu = dolfin.LUSolver(solve_mat)
        lu.parameters["reuse_factorization"] = True
        lu.parameters["report"] = False
        bb = dolfin.Vector(M.size(0))
        xx = dolfin.Vector(M.size(0))
        def sigma_solve(b):
            bb[:] = b
            lu.solve(xx, bb)
            return xx[:]
        M_matvec = lambda x: M*x
        eigs_w, eigs_v = speigs.ARPACK_gen_eigs(
            M_matvec, sigma_solve, M.size(0), self.sigma, nev, ncv=ncv)

        return eigs_w, eigs_v.T
