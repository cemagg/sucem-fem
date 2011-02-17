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
from scipy import linalg
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Utilities import Struct, partial
from NewCode import DifferentialForm, Waveforms, PostProc
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
dt = 0.01
cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree
# mixed 4th
orders = {'E':(4,True), 'H':(4,True), 'D':(3,False), 'B':(3,False)}
# mixed 3rd
#orders = {'E':(3,True), 'H':(3,True), 'D':(2,False), 'B':(2,False)}
# LT/QN E, H; Linear D, B
#orders = {'E':(2,True), 'H':(2,True), 'D':(1,False), 'B':(1,False)}
# All lowest order
#orders = {'E':(1,True), 'H':(1,True), 'D':(1,True), 'B':(1,True)}
coupledSystem = CoupledFirstOrderSystem(
    mesh, BCs=Struct(E=cb, H=free, D=free, B=cb),
    disc_orders=orders)

totalDOFs_D = coupledSystem.discs.D.totalDOFs
totalDOFs_E = coupledSystem.discs.E.totalDOFs

E,D,H,B = [coupledSystem.discs[dn] for dn in ('E', 'D', 'H', 'B')]

M_h = H.matrix.mass()
print "M_h dofs: ", M_h.shape[0]
M_e = E.matrix.mass()
print "M_e dofs: ", M_e.shape[0]
C_e = E.matrix.exteriorDerivative(B)
C_h = H.matrix.exteriorDerivative(D)
P_de = D.matrix.projectionOnto(E)
P_bh = B.matrix.projectionOnto(H)


#print "(nodofs, nnz, sparsity %)", B.shape[0], B.nnz, B.nnz/B.shape[0]**2*100

print "LU factorzing M_e"
M_e_LU = scipy.linsolve.factorized(M_e)
# print "LU factorzing M_h"
# M_h_LU = scipy.linsolve.factorized(M_h)
print "Inverting M_h"
M_inv = N.mat(linalg.inv(M_h.todense(), overwrite_a=True))
print "Constructing pseudo_S"
print 'P_bh*C_e'
pseudo_S = (P_bh.matmat(C_e)).todense() 
print 'M_inv*P_bh*C_e'
pseudo_S = M_inv*pseudo_S
del M_inv
print 'P_de*C_h*M_inv*P_bh*C_e'
pseudo_S = (P_de*C_h).todense()*pseudo_S
print "pseudo_S done"

def matvec(x):
    return M_e_LU(N.dot(pseudo_S.A, x))

from scipy.sandbox.arpack import speigs

# print 'Solving Eigenproblem'
# w, v = speigs.ARPACK_eigs(matvec, M_e.shape[0], 1, which='LM', ncv=31)

# max_dt = 2/N.sqrt(w)
# print orders
# print max_dt[0]
 
sigma = 0.0304 ; which = 'LM' ; nev = 5 # For finding modes close to TE101
#sigma = 0.0005 ; which = 'LR' ; nev = 11 # For finding modes > 0
print "LU decomposition for sigma_solve"
lu,piv = linalg.lu_factor(pseudo_S - sigma*M_e, overwrite_a=1)

sigma_solve = lambda b: linalg.lu_solve((lu, piv), b)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M_e.matvec, sigma_solve, M_e.shape[0],
                             sigma, nev, ncv=91, which=which)

from AnalyticResults import accoustic_eigs, PEC_cavity
an_res = PEC_cavity['rect-white']
 
#pickle.dump(dict((k,val) for k,val in zip(w,v.T)), file('+coupleda_4th_eigmodes.pickle', 'w'))
