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
from NewCode import DifferentialForm, Waveforms, PostProc
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh, TwoTets
from NewCode.DiscretisedSystem import CurlCurlNewmark


mesh1 = Mesh.Mesh(FlatTet.listmesh)
mesh2 = Mesh.Mesh(TwoTets.listmesh)
mesh3 = Mesh.Mesh(InscribedTetMesh.listmesh)
eMAGUSImport.init('workspace')
mesh4 = Mesh.Mesh(eMAGUSImport.get_listmesh())
mesh=mesh4

disc1f = Discretiser.setup_PformDiscretiser(mesh, 1, order=3, mixed=True)
#disc2f = Discretiser.setup_PformDiscretiser(mesh, 2, mixed=False)

print 'Mesh elements: ', len(mesh.elements)

import random
from scipy.linalg import iterative
from scipy.linsolve import splu

A = disc1f.matrix.mass()
b = N.array([random.random() for x in xrange(A.shape[0])], N.float64)
def null_psolve(b):
    #print "null psolve %d" % null_psolve.cnt
    null_psolve.cnt += 1
    return b
null_psolve.cnt = 0

A.psolve = null_psolve
print "Solving w/o preconditioner"
x_nopred, info = iterative.cg(A, b, tol=1e-10)
print "Done"
import getfem, numarray
A_getfem = getfem.Spmat('empty', A.shape[0])
print "Filling getfem matrix"
def mm_write(fd, coo_mat):
    fd.write('%%MatrixMarket matrix coordinate real general\n')
    fd.write('%d %d %d\n' % (coo_mat.shape[0], coo_mat.shape[1], coo_mat.nnz))
    for row, col, val in izip(coo_mat.row+1, coo_mat.col+1, coo_mat.data):
        fd.write('%d %d %f\n' % (row, col, val))
fname = '/tmp/mat%d' % random.randint(0, sys.maxint)
mm_write(file(fname, 'w'), A.tocoo())
A_getfem = getfem.Spmat('load', 'mm', fname)
print "Calculating ILU"
A_getfem_pred = getfem.Precond('ilu', A_getfem)
def psolve(b):
    #print "getfem ILU psolve %d" % psolve.cnt
    psolve.cnt += 1
    M_inv_b = A_getfem_pred.mult(numarray.array(b))
    return N.array(M_inv_b).reshape(len(M_inv_b))
psolve.cnt = 0
A.psolve = psolve

print "Solving with ILU"
x_ilu, info = iterative.cg(A, b, tol=1e-10)

print "nopred: ", N.max(N.abs(A.matvec(x_nopred) - b)), null_psolve.cnt
print "ilu:    ", N.max(N.abs(A.matvec(x_ilu) - b)), psolve.cnt
