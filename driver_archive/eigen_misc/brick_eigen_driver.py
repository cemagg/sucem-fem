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
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode import DifferentialForm, Utilities
from NewCode.DifferentialForm import BrickDiscretiser

h0 = 1/30.
a,b,c = 29,23,19

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

print 'Mesh elements: ', len(mesh.elements)

class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser

class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
    DiscretiserModule = BrickDiscretiser

#eigsys = BrickAudioEigen(mesh, 2, mixed=True)
eigsys = BrickVectorWaveEigen(mesh, 1, mixed=True)
from NewCode.Integration import BrickTrapzIntegrator

eigsys.disc.set_integrator(BrickTrapzIntegrator)
eigsys.disc.setIntegrationRule(0)
eigsys.disc.D().set_integrator(BrickTrapzIntegrator)
eigsys.disc.D().setIntegrationRule(0)

print 'Getting S'
S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()

from scipy.sandbox.arpack import speigs
sigma = 0.01
def sigma_solve_gen(S):
    from NewCode.MatrixUtils import MatrixSolver
    ms = MatrixSolver(5e-9)
    cnt = [0]
    def solve(b):
        print 'Starting iterative solve'
        x = ms.solve_mat_vec(S, b)
        cnt[0] += 1
        print 'Done with iterative solve %i' %cnt[0]
        return x
    return solve
#sigma_solve = sigma_solve_gen(S - sigma*M)
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2*100
print 'Sparse LU decomposition'
sigma_solve = scipy.linsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 21, ncv=201)

from AnalyticResults import PEC_cavity, err_percentage

ares = PEC_cavity['rect-white']
err = err_percentage(ares, sorted(w[w > 0.0000001])[0:10])
RMS_err = Utilities.RMS(err)
#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']
