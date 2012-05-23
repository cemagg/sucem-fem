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
from sucemfem import Forms 
from sucemfem import Materials 
from sucemfem import BoundaryConditions 
from sucemfem import SystemMatrices

#from scipy.sparse.linalg.eigen.arpack import speigs
#from scipy.sparse.linalg.eigen.arpack import eigs

from scipy.sparse.linalg import arpack
#from scipy.sparse.linalg import eigs


from EMProblem import EMProblem

import numpy as np

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
        """Solve problem and return the eigenvalues and eigenvectors
        
        @param nev: Number of eigenpairs to compute
        @keyword ncv: Number of Arnoldi basisvectors to use. 
            (default: 2*nev+1).
        
        @rtype: (C{numpy.array}, C{numpy.array})
        @return: (eigs_w, eigs_v) -- A tupple consisting of the eigenvalues and eigenvectors of the eigensystem.
            The eigenvalues are returned as an array of n_eig values corresponding to k^2 for the mode's resonant wavenumber.
            The eigenvectors are returned as a 2D array of shape (n_eig, problem_dim), with row i corresponding to the 
            modal distributions associated with the i-th eigenvalue. 
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
              
        #speigs in ARPACK has been removed in scipy 0.9/0.10
        #eigs is now in scipy.sparse.linalg.arpack
        #therefore force M_matvec to have "shape" and "dtype" attributes to comply with "eigs" requirements
        #also perform the sigma_solve function for spectrum shift to sigma
        class RM:
            def __init__(self, M):
                self.M = M
                self.shape = (M.size(0), M.size(1))
                self.dtype = np.dtype('d')                              
            def matvec(self, x):                                
                return sigma_solve(self.M*x)
                               
#        eigs_w, eigs_v = speigs.ARPACK_gen_eigs(
#            M_matvec, sigma_solve, M.size(0), self.sigma, nev, ncv=ncv)
                                            
        eigs_w, eigs_v= arpack.eigs(RM(M), k=nev, sigma=self.sigma, which='LM', ncv=ncv)
      
        return eigs_w, eigs_v.T