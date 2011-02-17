from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip, chain
import os, sys
import pickle
import random
import numpy as N
import scipy
from scipy import sparse, linalg
from scipy.sparse.linalg import iterative

from numpy.testing import assert_almost_equal
#
# Local Imports
#
import NewCode
import NewCode.eMAGUSImport as eMAGUSImport
from NewCode import Utilities, ProxyList
from NewCode.Utilities import Struct,  close_to_point
import NewCode.Mesh as Mesh
from NewCode.Meshes import BrickMesh, Conversions, BrickMeshGen
from NewCode import SubDimMesh, BrickSubDimMesh
from NewCode.DifferentialForm import Discretiser, BrickDiscretiser, allfree, HybridMesh
from NewCode.DifferentialForm import SubDimDiscretiserEntities, SubDimDiscretiser
from NewCode.DifferentialForm import BrickSubDimDiscretiser, BrickSubDimDiscretiserEntities
from NewCode.DifferentialForm import constrained_on_boundary
from NewCode.Integration import TriIntegrator
from NewCode.SystemMatrix import local_self_projection_matrix, insert_global, set_global, eps

h0 = 1.0001/2.
a,b,c = 29,23,19
bfrac = 0.25
brick_mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a*bfrac,b,c, [a*h0, b*h0, c*h0]))

mtet_order = 3
order = 2

# mtet_order = 1
# order = 1
eMAGUSImport.init('workspace')
tet_mesh = Mesh.Mesh(eMAGUSImport.get_listmesh())

print 'Tet_mesh elements: ', len(tet_mesh.elements)

print 'Hex-Mesh elements: ', len(brick_mesh.elements)

g_eps = 1e-10                           # Geometrical tollerance
hybrid_boundary_p = close_to_point(a*bfrac, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(x) for (x,y,z) in ent.nodeCoords])

zero_p = close_to_point(0, g_eps)
a_p, b_p, c_p = (close_to_point(xx, g_eps) for xx in (a,b,c))
def freeE(ent):
    x,y,z = ent.nodeCoords.T
    return not (N.all(a_p(x)) or N.all(zero_p(x)) or 
                N.all(b_p(y)) or N.all(zero_p(y)) or
                N.all(c_p(z)) or N.all(zero_p(z)))

#freeE = allfree
otet_freefun = lambda ent: not on_hbdry(ent) and freeE(ent)

mtet_disc = Discretiser.setup_PformDiscretiser( # mtet for Matching Tet
    tet_mesh, form=1, order=mtet_order, mixed=False, freeFun=freeE)
mtet_disc.setIntegrationRule(6)
mtet_disc.setFaceIntegrationRule(TriIntegrator(6))
otet_disc = Discretiser.setup_PformDiscretiser( # otet for Ordinary Tet
    tet_mesh, form=1, order=order, mixed=True, freeFun=otet_freefun)
brick_disc = BrickDiscretiser.setup_PformDiscretiser(
    brick_mesh, form=1, order=order, mixed=True, freeFun=freeE)
brick_disc.setIntegrationRule(3)

mtet_dofnos = Struct((k,mtet_disc.permuter.globalEntityPermutationTable(k))
                     for k in mtet_disc.geomEntities)
brick_dofnos = Struct((k,brick_disc.permuter.globalEntityPermutationTable(k))
                     for k in brick_disc.geomEntities)

mtet_submesh = SubDimMesh.SubSurface(tet_mesh, faceSelector=on_hbdry)
mtet_subGeomEntities = {
    'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
    SubDimDiscretiserEntities.Edge(mesh=mtet_submesh, freefun=allfree,
                                   attrs=mtet_submesh.edges.list_repr()))}
if 'face' in mtet_disc.basisSet.fns:
    mtet_subGeomEntities['face'] = SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
        SubDimDiscretiserEntities.Face(mesh=mtet_submesh, freefun=allfree,
                                       attrs=mtet_submesh.elements.list_repr()))

mtet_subdisc = SubDimDiscretiser.PformSubDimDiscretiser(
    1, mtet_submesh, mtet_subGeomEntities, Discretiser.Permuter, mtet_disc)

brick_submesh = BrickSubDimMesh.SubSurface(brick_mesh, faceSelector=on_hbdry)
brick_subGeomEntities = {
    'edge': BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
    BrickSubDimDiscretiserEntities.Edge(mesh=brick_submesh, freefun=allfree,
                                        attrs=brick_submesh.edges.list_repr()))}
if 'face' in brick_disc.basisSet.fns:
    brick_subGeomEntities['face'] = BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
        BrickSubDimDiscretiserEntities.Face(mesh=brick_submesh, freefun=allfree,
                                            attrs=brick_submesh.elements.list_repr()))
    
brick_subdisc = BrickSubDimDiscretiser.PformSubDimDiscretiser(
    1, brick_submesh, brick_subGeomEntities, Discretiser.Permuter, brick_disc)

mtet_nodofs = mtet_disc.totalDOFs
brick_nodofs = brick_disc.totalDOFs
hyb_nodofs = brick_subdisc.totalDOFs
otet_nodofs = otet_disc.totalDOFs 
tot_nodofs = otet_nodofs + brick_nodofs 

hc_T = sparse.lil_matrix(shape=(tot_nodofs, mtet_nodofs), dtype=N.float64)


def make_2face_permutation(tf_els):
    glob_perms = [el.permutation() for el in tf_els]
    unique_glob_dofs = N.unique(N.hstack([ep[1] for ep in glob_perms]))
    unique_glob_dofs = unique_glob_dofs[unique_glob_dofs >= 0]
    globto2face = dict((g,i) for i,g in enumerate(unique_glob_dofs))
    twotet_eldofs = [N.array([globto2face[g] for g in elperm[1]], N.int32)
                     for elperm in glob_perms]
    glob_perm = [N.arange(len(unique_glob_dofs)), unique_glob_dofs]
    return [(gp[0], tfp) for gp, tfp
            in zip(glob_perms, twotet_eldofs)], glob_perm
    

def get_d_nodes(hface_nodes):
    return hface_nodes[[0,2]]

def get_tface_nodes(hface_nodes):
    return (hface_nodes[[0,3,2]], hface_nodes[[0,1,2]])

subdim_ent_attr = dict(edge='edges', face='elements')
mtet_subdisc_glob_dofnos = Struct((k, mtet_dofnos[k][
    getattr(mtet_submesh, subdim_ent_attr[k])[:].superNo])
                                  for k in mtet_subdisc.geomEntities)
mtet_subdisc_dofnos = Struct((k,mtet_subdisc.permuter.globalEntityPermutationTable(k))
                     for k in mtet_subdisc.geomEntities)

mtet_subdim2glob_dofmap = dict((sd,gd) for sd,gd in izip(
    chain(*(nos.flat for nos in mtet_subdisc_dofnos.values())),
    chain(*(nos.flat for nos in mtet_subdisc_glob_dofnos.values()))))

brick_subdisc_glob_dofnos = Struct((k, brick_dofnos[k][
    getattr(brick_submesh, subdim_ent_attr[k])[:].superNo])
                                  for k in brick_subdisc.geomEntities)
brick_subdisc_dofnos = Struct((k,brick_subdisc.permuter.globalEntityPermutationTable(k))
                     for k in brick_subdisc.geomEntities)

brick_subdim2glob_dofmap = dict((sd,gd) for sd,gd in izip(
    chain(*(nos.flat for nos in brick_subdisc_dofnos.values())),
    chain(*(nos.flat for nos in brick_subdisc_glob_dofnos.values()))))

def rem_constr_perm(perm):
    pl, pg = perm
    pl_n = pl[pg >= 0]
    pg_n = pg[pg >= 0]
    return [pl_n, pg_n]

tsub_facemap = mtet_submesh.elements.nodemap
for hface in brick_subdisc.elements:
    tsub_facenodes = get_tface_nodes(hface.nodes)
    tsub_facenos = [tsub_facemap[tuple(nds)] for nds in tsub_facenodes]
    perm_2tf_els, twoface_to_tgsubdofs = make_2face_permutation(
        mtet_subdisc.elements[tsub_facenos])
    no_2tf_dofs = len(twoface_to_tgsubdofs[0])
    l2_M = N.zeros(shape=(no_2tf_dofs, no_2tf_dofs),
                   dtype=N.float64)
    l2_P = N.zeros(shape=(no_2tf_dofs, hface.noDOFs.element),
                   dtype=N.float64)
    hface_perm = hface.permutation()
    for el_i, el in enumerate(mtet_subdisc.elements[tsub_facenos]):
        tl_perm = el.permutation()[0]
        l_M = local_self_projection_matrix(el)
        insert_global(l2_M, l_M, perm_2tf_els[el_i])
        eval_r = el.physEvalPoints()
        eval_hface_l = [hface.global2local(r) for r in eval_r]
        hface_physvals = hface.physValsAtPoints(eval_hface_l)
        intg = el.rule.integrateFun
        l_P = N.array([[intg(N.sum(fn_i*fn_j, axis=1))
                        for fn_j in hface_physvals]
                       for fn_i in el.physVals()], N.float64)
        l_P *= el.size
        l_P[N.abs(l_P) < eps] = 0
        insert_global(l2_P, l_P, perm_2tf_els[el_i],
                      [N.arange(l_P.shape[1]),
                       N.arange(l_P.shape[1])])
    l2_M_inv = linalg.inv(l2_M)
    l2_hc_T = (N.dot(l2_M_inv, l2_P)).T
    l2_hc_T[N.abs(l2_hc_T) < eps] = 0
    l2_tet_perm = rem_constr_perm([twoface_to_tgsubdofs[0], N.array(
        [mtet_subdim2glob_dofmap[i] for i in twoface_to_tgsubdofs[1]], N.int32)])
    l2_brick_perm = rem_constr_perm([hface_perm[0], N.array(
        [brick_subdim2glob_dofmap[i] for i in hface_perm[1]], N.int32)])
    set_global(hc_T, l2_hc_T, l2_brick_perm, l2_tet_perm)


hct_dofnos = Struct((k,otet_disc.permuter.globalEntityPermutationTable(k))
                    for k in otet_disc.geomEntities)

for ent, hct_ent_dofnos in hct_dofnos.items():
    hct_ent_dofnos[hct_ent_dofnos >= 0] += brick_nodofs
    
for ent, entmap in hct_dofnos.items():
    for ent_i, ent_t_map in enumerate(entmap):
        hct_ent_dofno = hct_dofnos[ent][ent_i]
        for dof_j, hcf_dn in enumerate(hct_ent_dofno):
            if hcf_dn == -1: continue
            hc_T[hcf_dn, mtet_dofnos[ent][ent_i][dof_j]] = 1

hc_T = hc_T.tocsc()

brickglob_T = sparse.lil_matrix(shape=(tot_nodofs, brick_nodofs), dtype=N.int8)
for i in range(brick_nodofs):
    brickglob_T[i,i] = 1
brickglob_T = brickglob_T.tocsc()

M_mtet = mtet_disc.matrix.mass()
S_mtet = mtet_disc.matrix.stiffness()
M_brick_l = brick_disc.matrix.mass()
S_brick_l = brick_disc.matrix.stiffness()


M_brick = brickglob_T.matmat(M_brick_l.matmat(brickglob_T.T)).tocsc()
S_brick = brickglob_T.matmat(S_brick_l.matmat(brickglob_T.T)).tocsc()


M_hc = hc_T.matmat(M_mtet.matmat(hc_T.T)).tocsc()
S_hc = hc_T.matmat(S_mtet.matmat(hc_T.T)).tocsc()

M = M_hc + M_brick
S = S_hc + S_brick


from scipy.sparse.linalg.eigen.arpack import speigs
sigma = 0.01
print "(nodofs, nnz, sparsity %)", M.shape[0], M.nnz, M.nnz/M.shape[0]**2.0*100
print 'Sparse LU decomposition'
sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - sigma*M)

def sigma_solve_gen(A):
    from NewCode.MatrixUtils import MatrixSolver
    tol = 5e-9
    cnt = [0]
    itercnt = [0]
    import inspect
    def icb(x):
        f = inspect.currentframe()
        resid = f.f_back.f_locals['resid']
        print 'solve no: ', cnt[0], ' iteration no: ', itercnt[0],  ' resid: ', resid
        itercnt[0]+= 1
    def solve(b):
        print 'Starting iterative solve'
        x, info = iterative.cg(A, b, tol=tol, maxiter=A.shape[0]*10, callback=icb)
        if info > 0: raise Exception('Iterative solve did not converge')
        cnt[0] += 1
        print 'Done with iterative solve %i' %cnt[0]
        itercnt[0] = 0
        return x
    return solve


#sigma_solve = sigma_solve_gen(S - sigma*M)

print 'Solving Eigenproblem'
w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0], sigma, 51, ncv=91)

#w,v = scipy.linalg.eig(S.todense(), M.todense())

res =  N.array(sorted(N.abs(w[w > 0.0000001]))[0:10])

from AnalyticResults import PEC_cavity, accoustic_eigs, err_percentage

ares = PEC_cavity['rect-white']

err = err_percentage(ares, res)
RMS_err = Utilities.RMS(err)

# err_mtet = err_percentage(ares, res_mtet)
# RMS_err_mtet = Utilities.RMS(err_mtet)

print RMS_err
print err

