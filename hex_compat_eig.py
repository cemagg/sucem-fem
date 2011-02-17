from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
from scipy import sparse, linalg
from numpy.testing import assert_almost_equal
#
# Local Imports
#
import NewCode
from NewCode import Utilities
from NewCode.Utilities import Struct
from NewCode.tests.PyramMeshes import SixPyram as TestMesh
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Meshes import PyramMesh, Conversions, BrickMeshGen
from NewCode.DifferentialForm import Discretiser, allfree

order = 1

eMAGUSImport.init('workspace')
mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())
h0 = 1.0001/20.
a,b,c = 29,23,19
pyr_mesh = PyramMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

#pyr_mesh = PyramMesh.Mesh(TestMesh.listmesh)

print 'Mesh elements: ', len(mesh.elements)

d_edgemap = Conversions.hexfaces2diagnodes(pyr_mesh.basefaces)
edgemap = dict((tuple(nds), i) for i, nds in enumerate(mesh.edges[:].nodes))
d_edges = set(edge for edge_nodes, edge in edgemap.iteritems()
                         if tuple(edge_nodes) in d_edgemap)

disc = Discretiser.setup_PformDiscretiser(
    mesh, form=1, order=order, mixed=False, freeFun=allfree)

edge_dofnos = disc.permuter.globalEntityPermutationTable('edge')

no_edges = len(mesh.edges)
no_d_edges = len(d_edgemap)

no_newdofs = no_edges - no_d_edges
no_olddofs = 2*no_edges

new_edge_dofnos = N.zeros((no_edges,1), dtype=N.int32)
dofno_n = 0
for i, edge in enumerate(mesh.edges):
    if i in d_edges: new_edge_dofnos[i,0] = -1
    else:
        new_edge_dofnos[i,0] = dofno_n
        dofno_n += 1



T = sparse.lil_matrix(shape=(no_newdofs, no_olddofs), dtype=N.float64)

for i, edge in enumerate(mesh.edges):
    if i in d_edges:
        (pen_1, pen_2), (nen_1, nen_2) = d_edgemap[tuple(edge.nodes)]
        penos = [edgemap[tuple(pen_1)], edgemap[tuple(pen_2)]]
        nenos = [edgemap[tuple(nen_1)], edgemap[tuple(nen_2)]]
        olddofno_d0, olddofno_d1 = edge_dofnos[i]
        pd1 = new_edge_dofnos[penos[0],0]
        pd2 = new_edge_dofnos[penos[1],0]
        nd1 = new_edge_dofnos[nenos[0],0]
        nd2 = new_edge_dofnos[nenos[1],0]
        T[[pd1, pd2, nd1, nd2], olddofno_d0] = 1/2
        T[[pd1, pd2], olddofno_d1] = 1/2
        T[[nd1, nd2], olddofno_d1] = -1/2
        continue
    T[new_edge_dofnos[i,0],edge_dofnos[i,0]] = 1
    
M_old = disc.matrix.mass()
S_old = disc.matrix.stiffness()

M = T.matmat(M_old.matmat(T.T))
S = T.matmat(S_old.matmat(T.T))

#w,v = scipy.linalg.eig(S.todense(), M.todense())
#w_old,v_old = scipy.linalg.eig(S_old.todense(), M_old.todense())

from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.01
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2*100
print 'Sparse LU decomposition'
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)

# print "(nodofs, nnz, sparsity %)", M.shape[0], M_old.nnz, M_old.nnz/M_old.shape[0]**2*100
# print 'Sparse LU decomposition'
# sigma_solve_old = scipy.sparse.linalg.dsolve.factorized(S_old - sigma*M_old)

# print 'Solving Eigenproblem'
# w_old,v_old = speigs.ARPACK_gen_eigs(M_old.matvec, sigma_solve_old, M_old.shape[0], sigma, 51, ncv=91)

# res_old =  N.array(sorted(N.abs(w_old[w_old > 0.0000001]))[0:10])
res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']

err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err)

print RMS_err
print err
