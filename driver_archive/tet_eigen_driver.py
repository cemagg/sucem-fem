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
from NewCode.SystemMatrix import local_self_projection_matrix

mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh = mesh4

order = 2
mixed = True

print 'Mesh elements: ', len(mesh.elements)
#eigsys = NewCode.DiscretisedSystem.VectorWaveEigen(mesh, order, mixed=mixed)
eigsys = NewCode.DiscretisedSystem.VectorWaveEigen(mesh, order, mixed=mixed, BC=lambda x:True)

eigsys.disc.set_integrator(Integration.TetProdIntegrator)
eigsys.disc.setIntegrationRule(10)
eigsys.disc.D().set_integrator(Integration.TetProdIntegrator)
eigsys.disc.D().setIntegrationRule(10)

print 'Getting S'
S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()

from scipy.sparse.linalg.eigen.arpack import speigs

sigma = 0.01
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
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)
#w,v = scipy.linalg.eig(S.todense(), M.todense())
# w_m, v_m = scipy.linalg.eig(M.todense())
res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])
#res =  N.array(sorted(N.abs(w[w > 0.0000001])))

print res

from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']
# #ares = accoustic_eigs['rect1x0.75x0.5']
# #ares = PEC_cavity['rect1x0.25x5.1']

err = err_percentage(ares, res)
RMS_err = RMS(err)

print RMS_err
print err



