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
import gc
#
# Local Imports
#
# import NewCode
# import NewCode.eMAGUSImport as eMAGUSImport
# from NewCode.Meshes import BrickMesh, BrickMeshGen
# from NewCode import DifferentialForm, Utilities
# from NewCode.DifferentialForm import BrickDiscretiser
# from NewCode.DifferentialForm import constrained_on_boundary


# h0 = 1/5.
# a,b,c = 29,23,19
# #a,b,c = 1.,0.75,0.5

# mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

# order = 3

# print 'Mesh elements: ', len(mesh.elements)

# class BrickVectorWaveEigen(NewCode.DiscretisedSystem.VectorWaveEigen):
#     DiscretiserModule = BrickDiscretiser
#     def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
#                  volBC=None):
#         self.order = order
#         self.mesh = mesh
#         self.disc = self.DiscretiserModule.setup_PformDiscretiser(
#             mesh, self.p, order, mixed, BC, volBC, btype='cohen98')
#         self.massMatrix = self.disc.matrix.mass
#         self.stiffnessMatrix = self.disc.matrix.stiffness

# class BrickAudioEigen(NewCode.DiscretisedSystem.AudioEigen):
#     DiscretiserModule = BrickDiscretiser
#     def __init__(self, mesh, order, mixed=True, BC=constrained_on_boundary,
#                  volBC=None):
#         self.order = order
#         self.mesh = mesh
#         self.disc = self.DiscretiserModule.setup_PformDiscretiser(
#             mesh, self.p, order, mixed, BC, volBC, btype='cohen98')
#         self.massMatrix = self.disc.matrix.mass
#         self.stiffnessMatrix = self.disc.matrix.stiffness

# #eigsys_B = BrickAudioEigen(mesh, order, mixed=True)
# eigsys_E = BrickVectorWaveEigen(mesh, order, mixed=True)
# E = eigsys_E.disc
# #E, B = eigsys_E.disc, eigsys_B.disc

# print 'Getting S'
# S = eigsys_E.stiffnessMatrix()
# print 'Getting M'
# M_e = E.matrix.mass()

# A = M_e + 10*S
from scipy.sparse.linalg import iterative

from scipy.io import mmread
A = mmread(file('+o3_h4_refined1m_A_imp_o3_dt32.mtx'))
from scipy.sparse.linalg.dsolve import factorized

A_csc = A.tocsc()
A_lu = factorized(A_csc)
umfcon = A_lu.func_closure[1].cell_contents

 # A = A.tocsr()
b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)



# def null_psolve(b):
#     #print "null psolve %d" % null_psolve.cnt
#     null_psolve.cnt += 1
#     return b
# null_psolve.cnt = 0

# A.psolve = null_psolve
# # print "Solving w/o preconditioner"
# # x_nopred, info = iterative.cg(A, b, tol=1e-10)
# # print "Done"


# A.ensure_sorted_indices(inplace=True)
# import petsc4py
# petsc4py.init(sys.argv)
# from petsc4py import PETSc

# A_petsc = PETSc.Mat().createAIJ(A.shape,
#                           csr=(A.indptr,
#                                A.indices,
#                                A.data))
# # obtain vectors for storing
# # the solution  and the rhs
# x_petsc, b_petsc = A_petsc.getVecs()
# # fill the rhs PETSc vector
# # from the rhs numpy array
# b_petsc[...] = b # note the syntax sugar
# PETSC_DECIDE      =   -1
# OptDB = PETSc.Options()
# for key in OptDB.getAll().keys(): OptDB.delValue(key)
# OptDB['ksp_type'] = 'cg'
# # OptDB['ksp_type']          = 'gmres'
# # OptDB['ksp_gmres_restart'] = 50
# # OptDB['pc_type'] = 'none'
# #OptDB['pc_type']  = 'icc'
# #OptDB['pc_type']  = 'ilu'
# #OptDB['pc_factor_levels'] = 5
# #OptDB['pc_factor_shift_positive_definite'] = 1
# #OptDB['pc_factor_shift_nonzero'] = PETSC_DECIDE
# OptDB['pc_type']        = 'hypre'
# OptDB['pc_hypre_type']  = 'boomeramg'
# OptDB['pc_hypre_boomeramg_strong_threshold'] = 0.5
# OptDB['pc_hypre_boomeramg_max_iter'] = 2
# ksp = PETSc.KSP().create()
# ksp.setFromOptions()  # configure from OptDB
# ksp.setOperators(A_petsc)   
# ksp.setTolerances(1e-10)
# ksp.solve(b_petsc, x_petsc) 

# x_ilu = x_petsc[...]
# print "   ilu:    ", N.max(N.abs(A.matvec(N.array(x_ilu)) - b))
# #print "nopred:    ", N.max(N.abs(A.matvec(N.array(x_nopred)) - b))
# print ksp.getPC().view()


# for i in range(10000):
#     #psolve.cnt = 0
#     #b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)
#     x_ilu = getfem.linsolve_cg(A_getfem, b_na, A_getfem_pred)
#     print "i: ", i, " ilu:    ", N.max(N.abs(A.matvec(N.array(x_ilu)) - b))
#     sys.stdout.flush()

