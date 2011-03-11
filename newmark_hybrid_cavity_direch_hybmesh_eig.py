"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from __future__ import division
import os
import pickle
import random
import numpy as N
import scipy
from scipy import sparse

#
# Local Imports
#
import NewCode
from NewCode import Utilities
import NewCode.eMAGUSImport as eMAGUSImport
import NewCode.Mesh as Mesh
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm, Waveforms, PostProc, Feeds, Runners
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser, \
     SubDimDiscretiserEntities, allconstrained, allfree
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet, BrickCurlCurlNewmarkDirechlet 

from NewCode.Feeds import WaveguideEigenMatcher

from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridSystem \
     import HybridMeshNewmarkHybridSystem
from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridMatrices \
     import HybridMeshHybridBlockMats

eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())

print 'Tet_mesh elements: ', len(tet_mesh.elements)

# h0 = 1.0001/1.
# ha = 1.0001/4.
a,b,c = 29,23,19
h = a/2, b, c
bfrac = 0.5
order = 3

# brick_mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(a*bfrac,b,c, [a*ha, b*h0, c*h0]))

brick_mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a*bfrac,b,c,h))

print 'Brick-Mesh elements: ', len(brick_mesh.elements)


g_eps = 1e-10                           # Geometrical tollerance
hybrid_boundary_p = close_to_point(a*bfrac, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(x) for (x,y,z) in ent.nodeCoords])
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
zero_p = close_to_point(0, g_eps)

def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))

#freeE = allfree

hyb_mesh = Struct(tet=tet_mesh, brick=brick_mesh, on_hbdry=on_hbdry)

system = HybridMeshNewmarkHybridSystem(tet_mesh, brick_mesh, implicit_beta=0.)
system.init_group_freefuns(on_hbdry, freeE)
system.init_discs(order)
system.init_block_matrices()
system.init_merged_mats()
system.init_dofs()
system.set_dt(1.)

M = system.merged_matrices.A()
S = 2*M - system.merged_matrices.B()

from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.01
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2.0*100
print 'Sparse LU decomposition'
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
#w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)

w,v = scipy.linalg.eig(S.todense(), M.todense())

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

print res


from postproc_code.AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']

err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err[0:4])

# err_mtet = err_percentage(ares, res_mtet)
# RMS_err_mtet = Utilities.RMS(err_mtet)

print RMS_err
print err
print 'Min real eigvalue: ', N.min(N.real(w))

d_e = system.discs.E.e
