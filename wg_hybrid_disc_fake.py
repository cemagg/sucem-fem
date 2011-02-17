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
from NewCode.DifferentialForm import Discretiser, allfree, HybridMesh
from NewCode.DifferentialForm import constrained_on_boundary
from NewCode.SystemMatrix import set_global



def setup_PformDiscretiser(
    mesh, form, order=1, mixed=True, freeFun=allfree, vol_freeFun=None,
    dead_elements=None):
    assert(mixed==True)                 # Can't do full order hybrid discretisers
    hyb_disc = HybridMesh.HybridMeshOneformDiscretiser(mesh)
    hyb_disc.init_discs(order, freeFun, dead_elements=dead_elements,
                        constrained_mtet_vols=False)
    hyb_disc.init_subdiscs()
    hyb_disc.init_entity_dofnos()
    hyb_disc.init_nodofs()
    
    hc_T = sparse.lil_matrix(shape=(hyb_disc.nodofs.tot, hyb_disc.nodofs.mtet), dtype=N.float64)
    
    hyb_dofmaps = HybridMesh.HybridSubdim2GlobDofmaps()
    subdim_ent_attr = hyb_dofmaps.subdim_ent_attr
    
    hyb_dofmaps.init_dofmaps(hyb_disc.subdiscs, hyb_disc.entity_dofnos)
    
    
    mtet_subdim2glob_dofmap = hyb_dofmaps.dofmap.mtet
    brick_subdim2glob_dofmap = hyb_dofmaps.dofmap.brick
    hextrimatcher = HybridMesh.HexTriFaceMatcher(hyb_disc.submeshes)
    hface_tface_map = hextrimatcher.hex_tris

    print 'hexface supernos: ', hyb_disc.subdiscs.brick.elements[:].superNo
    #tsub_facemap = hyb_disc.subdiscs.mtet.mesh.elements.nodemap
    for hface in hyb_disc.subdiscs.brick.elements:
        #tsub_facenodes = HybridMesh.get_tface_nodes(hface.nodes)
        #tsub_facenos = [tsub_facemap[tuple(nds)] for nds in tsub_facenodes]
        tsub_facenos = hface_tface_map[hface.index]
        hface_perm = hface.permutation()
        l2_hc_T, twoface_to_tgsubdofs = HybridMesh.make_hface_locmat(
            hyb_disc.subdiscs.mtet.elements[tsub_facenos], hface)
        l2_tet_perm = HybridMesh.rem_constr_perm([twoface_to_tgsubdofs[0], N.array(
            [mtet_subdim2glob_dofmap[i] for i in twoface_to_tgsubdofs[1]], N.int32)])
        l2_brick_perm = HybridMesh.rem_constr_perm([hface_perm[0], N.array(
            [brick_subdim2glob_dofmap[i] for i in hface_perm[1]], N.int32)])
        set_global(hc_T, l2_hc_T, l2_brick_perm, l2_tet_perm)
    
    
    hct_dofnos = Struct((k,hyb_disc.discs.otet.permuter.globalEntityPermutationTable(k))
                        for k in hyb_disc.discs.otet.geomEntities)
    
    for ent, hct_ent_dofnos in hct_dofnos.items():
        hct_ent_dofnos[hct_ent_dofnos >= 0] += hyb_disc.nodofs.brick
        
    for ent, entmap in hct_dofnos.items():
        for ent_i, ent_t_map in enumerate(entmap):
            hct_ent_dofno = hct_dofnos[ent][ent_i]
            for dof_j, hcf_dn in enumerate(hct_ent_dofno):
                if hcf_dn == -1: continue
                hc_T[hcf_dn, hyb_disc.entity_dofnos.mtet[ent][ent_i][dof_j]] = 1
    
    hc_T = hc_T.tocsc()
    
    brickglob_T = sparse.lil_matrix(shape=(hyb_disc.nodofs.tot, hyb_disc.nodofs.brick), dtype=N.int8)
    for i in range(hyb_disc.nodofs.brick):
        brickglob_T[i,i] = 1
    brickglob_T = brickglob_T.tocsc()
    
    M_mtet = hyb_disc.discs.mtet.matrix.mass()
    S_mtet = hyb_disc.discs.mtet.matrix.stiffness()
    M_brick_l = hyb_disc.discs.brick.matrix.mass()
    S_brick_l = hyb_disc.discs.brick.matrix.stiffness()
    
    
    M_brick = brickglob_T.matmat(M_brick_l.matmat(brickglob_T.T)).tocsc()
    S_brick = brickglob_T.matmat(S_brick_l.matmat(brickglob_T.T)).tocsc()
    
    
    M_hc = hc_T.matmat(M_mtet.matmat(hc_T.T)).tocsc()
    S_hc = hc_T.matmat(S_mtet.matmat(hc_T.T)).tocsc()
    
    M = M_hc + M_brick
    S = S_hc + S_brick

    ###
    # Add various Discretiser class mimicking attributes
    ###
    def compoundProjection(target_disc):
        if target_disc.mesh is mesh.brick:
            P_o = hyb_disc.discs.brick.matrix.projectionOnto(target_disc)
            return P_o.matmat(brickglob_T.T)
        elif target_disc.mesh is mesh.tet:
            P_o = hyb_disc.discs.tet.matrix.projectionOnto(target_disc)
            return P_o.matmat(hc_T.T)

    def compoundProjection_D(target_disc):
        if target_disc.mesh is mesh.brick:
            P_o = hyb_disc.discs.brick.D().matrix.projectionOnto(target_disc)
            return P_o.matmat(brickglob_T.T)
        elif target_disc.mesh is mesh.tet:
            P_o = hyb_disc.discs.tet.D().matrix.projectionOnto(target_disc)
            return P_o.matmat(hc_T.T)
            

    hyb_disc.matrix = Struct(mass=lambda : M, stiffness=lambda : S,
                             compoundProjection=compoundProjection)
    hyb_disc.stuff = Struct(hc_T=hc_T, brickglob_T=brickglob_T,
                            M_brick_l=M_brick_l, S_brick_l=S_brick_l,
                            M_brick=M_brick, S_brick=S_brick,
                            M_mtet=M_mtet, S_mtet=S_mtet,
                            M_hc=M_hc, S_hc=S_hc
                            )
    hyb_disc.newDOFs = lambda : Struct(
        dofArray=N.zeros(hyb_disc.nodofs.tot, N.float64),
        solver=NewCode.MatrixUtils.MatrixSolver())

    def D():
        return Struct(matrix=Struct(compoundProjection=compoundProjection_D))
    hyb_disc.D = D
    hyb_disc.totalDOFs = hyb_disc.nodofs.tot
    return hyb_disc

