"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
from numpy import *
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Utilities import RMS
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode import DifferentialForm
from NewCode.DifferentialForm import Discretiser
from NewCode import Integration, SystemMatrix, PostProc
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DifferentialForm.BasisFunction import Oneform

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

print 'Mesh elements: ', len(mesh.elements)
#eigsys = NewCode.DiscretisedSystem.AudioEigen(mesh, 1, mixed=False)
eigsys = NewCode.DiscretisedSystem.VectorWaveEigen(mesh, 4, mixed=True)
#eigsys = NewCode.DiscretisedSystem.VectorWaveEigen(mesh, 4, mixed=False, BC=lambda x:True)
print 'Getting S'
S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()
#w = scipy.linalg.eig(S, M, right=False, overwrite_a=True, overwrite_b=True)

from scipy.sandbox.arpack import speigs

#sigma = 0.030 ; which = 'LM' ; nev = 1 # For finding modes close to TE101
sigma = 0.030 ; which = 'LM' ; nev = 11 # For finding first 11 modes
#sigma = 0.0005 ; which = 'LR' ; nev = 11 # For finding modes > 0
def sigma_solve_gen(A):
    from NewCode.MatrixUtils import MatrixSolver
    ms = MatrixSolver(5e-9)
    cnt = [0]
    def solve(b):
        print 'Starting iterative solve'
        x = ms.solve_mat_vec(A, b)
        cnt[0] += 1
        print 'Done with iterative solve %i' %cnt[0]
        return x
    return solve
#sigma_solve = sigma_solve_gen(S - sigma*M)
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2*100
print 'Sparse LU decomposition'
sigma_solve = scipy.linsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, nev, ncv=21, which='LR')


from AnalyticResults import accoustic_eigs, PEC_cavity
an_res = PEC_cavity['rect-white']

# norm_err = N.abs(N.sqrt(sorted(w)[0:10])-N.sqrt(an_res))/N.sqrt(an_res)

