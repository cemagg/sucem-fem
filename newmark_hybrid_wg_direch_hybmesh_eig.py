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
     SubDimDiscretiserEntities, allconstrained
from NewCode.DiscretisedSystem import CurlCurlNewmarkDirechlet, BrickCurlCurlNewmarkDirechlet 

from NewCode.Feeds import WaveguideEigenMatcher

from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridSystem \
     import HybridMeshNewmarkHybridSystem
from NewCode.ImplicitExplicit.HybridMeshNewmarkHybridMatrices \
     import HybridMeshHybridBlockMats

import wg_hybrid_disc_fake

eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())

print 'Tet_mesh elements: ', len(tet_mesh.elements)

bfrac = 0.45
h = 1/2.
a,b,c = 1, 0.25, 5.1
bfrac = N.ceil(bfrac*c/h)*h/c

brick_mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c*bfrac, [h, h, h]))

print 'Brick-Mesh elements: ', len(brick_mesh.elements)


order = 3
free = DifferentialForm.allfree

g_eps = 1e-10                           # Geometrical tollerance
hybrid_boundary_p = close_to_point(c*bfrac, g_eps)
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
zero_p = close_to_point(0, g_eps)

def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))

from analytical_WG_driver import WGT

drv_fun = WGT.discrete_drv_fn
z_measure = WGT.test_z
z_measure_p = close_to_point(z_measure, g_eps)
runtime = WGT.t_final
z_port = 0.
z_port_p = close_to_point(z_port, g_eps)
a=1. ; b=0.25                           # WG x/y dim
zero_p = close_to_point(0, g_eps)
on_port = lambda ent: N.all(z_port_p(ent.nodeCoords[:,2]))
def port_free(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)))
analytical_dt = WGT.dt

on_measurement_port = lambda ent: N.all(z_measure_p(ent.nodeCoords[:,2]))
measurement_port_free = port_free
direch_free = lambda ent: on_port(ent) and port_free(ent)  

hybrid_boundary_p = close_to_point(c*bfrac, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(z) for (x,y,z) in ent.nodeCoords])
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

M_cc_exp = system.block_matrices.M_cc_exp()
# M_dd = system.block_matrices.A_dd()
# M_ee = system.block_matrices.A_ee()

S_cc_exp = system.block_matrices.S_cc_exp()
# S_cd = - system.block_matrices.B_cd()
# S_dd = 2*M_dd - system.block_matrices.B_dd()
# S_de = - system.block_matrices.B_de()
# S_ee = 2*M_ee - system.block_matrices.B_ee()

# M = sparse.bmat([[M_cc, None, None],
#                  [None, M_dd, None],
#                  [None, None, M_ee]]).tocsc()
# S = sparse.bmat([[S_cc, S_cd, None],
#                  [S_cd.T, S_dd, S_de],
#                  [None, S_de.T, S_ee]]).tocsc()

M_aa = system.block_matrices.A_aa()
# M_ab = system.block_matrices.A_ab()
# M_bb = system.block_matrices.A_bb()
# M_bc = system.block_matrices.A_bc()
M_cc = system.block_matrices.M_cc_imp()

S_aa = 2*M_aa - system.block_matrices.B_aa()
# S_ab = 2*M_ab - system.block_matrices.B_ab()
# S_bb = 2*M_bb - system.block_matrices.B_bb()
# S_bc = 2*M_bc - system.block_matrices.B_bc()
S_cc = system.block_matrices.S_cc_imp()

# M = sparse.bmat([[M_aa,   M_ab,   None],
#                  [M_ab.T, M_bb,   M_bc],
#                  [None,   M_bc.T, M_cc]]).tocsc()

# S = sparse.bmat([[S_aa,   S_ab,   None],
#                  [S_ab.T, S_bb,   S_bc],
#                  [None,   S_bc.T, S_cc]]).tocsc()

# M = sparse.bmat([[M_aa,   M_ab],
#                  [M_ab.T, M_bb]]).tocsc()

# S = sparse.bmat([[S_aa,   S_ab],
#                  [S_ab.T, S_bb]]).tocsc()

# M = sparse.bmat([[M_bb,   M_bc],
#                  [M_bc.T, M_cc]]).tocsc()

# S = sparse.bmat([[S_bb,   S_bc],
#                  [S_bc.T, S_cc]]).tocsc()


# S = S_bb
# M = M_bb


from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.5
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2.0*100
print 'Sparse LU decomposition'
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)

#w,v = scipy.linalg.eig(S.todense(), M.todense())

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

# print res



from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect1x0.25x5.1']

err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err)

# err_mtet = err_percentage(ares, res_mtet)
# RMS_err_mtet = Utilities.RMS(err_mtet)

print RMS_err
print err
print 'Min real eigvalue: ', N.min(N.real(w))

