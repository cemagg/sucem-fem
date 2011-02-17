"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc, MatrixUtils
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmark, CoupledFirstOrderSystem


mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

print 'Mesh elements: ', len(mesh.elements)

Discretiser.BasePformDiscretiser.defaultIntegrationOrder = 8

weights = 1.0
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
dt=0.01
# LT/QN E, H; Linear D, B
#orders = {'E':(2,True), 'H':(2,True), 'D':(1,False), 'B':(1,False)}
# All lowest order
orders = {'E':(4,True), 'H':(1,True), 'D':(1,True), 'B':(3,False)}
coupledSystem = CoupledFirstOrderSystem(
    mesh, dt=dt, BCs=Struct(E=cb, H=free, D=free, B=cb),
    disc_orders=orders)

totalDOFs_D = coupledSystem.discs.D.totalDOFs
totalDOFs_E = coupledSystem.discs.E.totalDOFs

E,D,H,B = [coupledSystem.discs[dn] for dn in ('E', 'D', 'H', 'B')]

M_e = E.matrix.mass()
#M_b = B.matrix.mass()
P_culreb = E.D().matrix.projectionOnto(B)

#import sys
#sys.path.append(ARPACK_modpath)
#import speigs
#sys.path.remove(ARPACK_modpath)
#ARPACK_modpath = '/home/nmarais/genugtig/akademie/scipy/ARPACK--simplify' 

from scipy.sandbox.arpack import speigs

#print "(nodofs, nnz, sparsity %)", B.shape[0], B.nnz, B.nnz/B.shape[0]**2*100

#M_e_LU = scipy.linsolve.factorized(M_e)
#M_b_LU = scipy.linsolve.factorized(M_b)

#matvec = lambda x: M_e_LU(P_culreb.T.matvec(M_b_LU(P_culreb.matvec(x))))

#matvec = lambda x: P_culreb.T.matvec(C.matvec(x))

C = E.matrix.exteriorDerivative(B)

S = P_culreb.T.matmat(C)
#S2 = E.matrix.stiffness()

sigma=3
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


#sigma_solve = sigma_solve_gen(Struct(matvec=lambda x: matvec(x) - sigma*M_e.matvec(x)))
print 'Constructing sigma_solve'
sigma_solve = scipy.linsolve.factorized(S - sigma*M_e)

print 'Solving Eigenproblem'

w,v = speigs.ARPACK_gen_eigs(M_e.matvec, sigma_solve, M_e.shape[0], sigma, 5, ncv=51)
#w, v = speigs.ARPACK_eigs(matvec, M_e.shape[0], 100)


from AnalyticResults import accoustic_eigs, PEC_cavity

