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
from NewCode.tests.BrickMeshes import OneBrick as TestMesh
from NewCode import DifferentialForm, Utilities
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DifferentialForm import constrained_on_boundary, allfree
from NewCode.Integration import HexaGaussLobattoLobattoIntegrator, \
     HexaLobattoGaussGaussIntegrator, HexaLobattoIntegrator

#h0 = 1/8.
a,b,c = 29,23,19
#a,b,c = 1.,0.75,0.5
#a,b,c = 1, 0.25, 5.1
h = a/8.
# mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h]))

#mesh = BrickMesh.Mesh(TestMesh.listmesh)

order = 3

freeE = constrained_on_boundary

print 'Mesh elements: ', len(mesh.elements)

class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, 'cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, 'cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

#eigsys = BrickAudioEigen(mesh, order, mixed=True)
eigsys = BrickVectorWaveEigen(mesh, order, mixed=True, BC=freeE)

# eigsys.disc.set_integrator(HexaLobattoGaussGaussIntegrator if order > 1
#                            else HexaLobattoIntegrator)
# eigsys.disc.set_integrator(HexaGaussLobattoLobattoIntegrator if order > 1
#                            else HexaLobattoIntegrator)
# eigsys.disc.setIntegrationRule(order if order > 1 else 2)
eigsys.disc.diagonalise()

print 'Getting S'
S = eigsys.stiffnessMatrix()
print 'Getting M'
M = eigsys.massMatrix()

from scipy.sparse.linalg.eigen.arpack import speigs
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
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)
#w,v = scipy.linalg.eig(S.todense(), M.todense())
from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']
#ares = accoustic_eigs['rect1x0.75x0.5']
#ares = PEC_cavity['rect1x0.25x5.1']

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])
err = err_percentage(ares, sorted(w[w > 0.0000001])[0:10])
RMS_err = Utilities.RMS(err[0:4])
#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']

print RMS_err
print err
print res
