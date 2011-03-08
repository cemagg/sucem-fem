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
from NewCode.DifferentialForm import constrained_on_boundary
from NewCode.Integration import HexaGaussLobattoLobattoIntegrator, \
     HexaLobattoGaussGaussIntegrator, HexaLobattoIntegrator

h0 = 1/5.
a,b,c = 29,23,19
#a,b,c = 1.,0.75,0.5

mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

order = 3

print 'Mesh elements: ', len(mesh.elements)

class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, btype='cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
    DiscretiserModule = BrickDiscretiser
    def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
                 volBC=None):
        self.order = order
        self.mesh = mesh
        self.disc = self.DiscretiserModule.setup_PformDiscretiser(
            mesh, self.p, order, mixed, BC, volBC, btype='cohen98')
        self.massMatrix = self.disc.matrix.mass
        self.stiffnessMatrix = self.disc.matrix.stiffness

eigsys_B = BrickAudioEigen(mesh, order, mixed=True)
eigsys_E = BrickVectorWaveEigen(mesh, order, mixed=True)
E, B = eigsys_E.disc, eigsys_B.disc

B.set_integrator(HexaLobattoIntegrator)
B.setIntegrationRule(order+1)
E.set_integrator(HexaLobattoIntegrator)
E.setIntegrationRule(order+1)
E.D().set_integrator(HexaLobattoIntegrator)
E.D().setIntegrationRule(order+1)


# print 'Getting S'
# S = eigsys.stiffnessMatrix()
print "P_mu"
P_mu = B.matrix.projectionOnto(E.D())
M_e = E.matrix.mass()
# print "M_b"
# M_b = B.matrix.mass()
# M_b_inv = scipy.sparse.dia_matrix( ([1/M_b.diagonal()], [0]), shape=M_b.shape)
print "Constructing C"
#C = M_b_inv.matmat(P_mu.T)
C = E.matrix.partialExteriorDerivative(B)
print "Constructing S"
S = P_mu.matmat(C)


M = scipy.sparse.dia_matrix( ([M_e.diagonal()], [0]), shape=M_e.shape)




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
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=151)

from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']
#ares = accoustic_eigs['rect1x0.75x0.5']
res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])
err = err_percentage(ares, sorted(w[w > 0.0000001])[0:10])
RMS_err = Utilities.RMS(err)
#edges = eigsys.disc.geomEntities['edge']
#faces = eigsys.disc.geomEntities['face']
