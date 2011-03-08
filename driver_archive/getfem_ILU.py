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
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode import DifferentialForm, Utilities
from NewCode.IOUtils import mm_write
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.DifferentialForm import constrained_on_boundary


# h0 = 1/10.
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


b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)
# def null_psolve(b):
#     #print "null psolve %d" % null_psolve.cnt
#     null_psolve.cnt += 1
#     return b
# null_psolve.cnt = 0

# A.psolve = null_psolve
# print "Solving w/o preconditioner"
# x_nopred, info = iterative.cg(A, b, tol=1e-10)
# print "Done"
import getfem, numarray
print "Filling getfem matrix"
fname = '/tmp/mat%d' % random.randint(0, sys.maxint)
#fname = '/
mm_write(file(fname, 'w'), A.tocoo())
#scipy.io.mmio.mmwrite(file(fname, 'w'), A, 'real', 16)
A_getfem = getfem.Spmat('load', 'mm', fname)

# print "A.to_coo()"
# coo_mat = A.tocoo()
# print "Filling getfem matrix"
# for  row, col, val in izip(coo_mat.row, coo_mat.col, coo_mat.data):
#     A_getfem[int(row), int(col)] = val
    
print "Calculating ILU"
A_getfem_pred = getfem.Precond('ilu', A_getfem)
#A_getfem_pred_fill = getfem.Precond('ilut', A_getfem)
A_getfem_pred_ildlt = getfem.Precond('ildlt', A_getfem)
#A_getfem_pred_ildlt_fill = getfem.Precond('ildltt', A_getfem)
#A_getfem_pred_spai = getfem.Precond('approx_inverse', A_getfem)
A_getfem_nullpred = getfem.Precond('identity')
def psolve(b):
    #print "getfem ILU psolve %d" % psolve.cnt
    psolve.cnt += 1
    b_na = numarray.array(b)
    del(b)
    M_inv_b_na = A_getfem_pred.mult(b_na)
    M_inv_b = N.array(M_inv_b_na).reshape(len(M_inv_b_na))
    #print "gc: ", gc.collect()
    del M_inv_b_na
    return M_inv_b

psolve.cnt = 0
A.psolve = psolve

# print "Solving with ILU"
# x_ilu, info = iterative.cg(A, b, tol=1e-10)

#print "nopred: ", N.max(N.abs(A.matvec(x_nopred) - b)), null_psolve.cnt
#print "ilu:    ", N.max(N.abs(A.matvec(x_ilu) - b))

b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)
b_na = numarray.array(b)



# for i in range(10000):
#     #psolve.cnt = 0
#     #b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)
#     x_ilu = getfem.linsolve_cg(A_getfem, b_na, A_getfem_pred, 'res', 1e-10)
#     print "i: ", i, " ilu:    ", N.max(N.abs(A.matvec(N.array(x_ilu)) - b))
#     sys.stdout.flush()

