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
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Meshes import BrickMesh
from NewCode import DifferentialForm
from NewCode.DifferentialForm import BrickDiscretiser
import brick_cavity

mesh = BrickMesh.Mesh(
    brick_cavity.make_rect_cavity_brick_listmesh(1,0.75,0.5,[1/5., 0.75/5., 0.5/5.]))

print 'Mesh elements: ', len(mesh.elements)

class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser

class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
    DiscretiserModule = BrickDiscretiser

eigsys = BrickAudioEigen(mesh, 2, mixed=True)
#eigsys = BrickVectorWaveEigen(mesh, 1, mixed=True)
#eigsys = NewCode.DiscretisedSystem.VectorWaveEigen(mesh, 4, mixed=False, BC=lambda x:True)
print 'Getting A'
A = eigsys.stiffnessMatrix()
print 'Getting B'
B = eigsys.massMatrix()

from scipy.sandbox.arpack import speigs
sigma = 1.
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
#sigma_solve = sigma_solve_gen(A - sigma*B)
print "(nodofs, nnz, sparsity %)", B.shape[0], B.nnz, B.nnz/B.shape[0]**2*100
print 'Sparse LU decomposition'
sigma_solve = scipy.linsolve.factorized(A - sigma*B)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(B.matvec, sigma_solve, B.shape[0], sigma, 10, ncv=501)


#from AnalyticResults import accoustic_eigs


#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']
els = eigsys.disc.elements
p = eigsys.disc.permuter
disc = eigsys.disc
