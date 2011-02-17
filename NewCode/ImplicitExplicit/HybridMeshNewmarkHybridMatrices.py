from __future__ import division

from itertools import izip
import numpy as N
from scipy import sparse
from NewCode.Utilities import CacheLast, Struct

from NewCode.SystemMatrix import set_global
from NewCode.DifferentialForm import HybridMesh
from NewCode.ImplicitExplicit.NewmarkHybridMatrices import HybridBlockMats

class HybridMeshHybridBlockMats(HybridBlockMats):
    no_dt_dep = HybridBlockMats.no_dt_dep | set(['hybrid_transform_mat'])
    def __init__(self, discs, implicit_beta, on_hbdry):
        self.discs = discs
        self.on_hbdry = on_hbdry
        self.implicit_beta = implicit_beta

    @CacheLast.CachedMethod
    def hybrid_transform_mat(self):
        hyb_nodofs = self.discs.E.c_exp.totalDOFs
        hybtet_nodofs = self.discs.E.c_imp.totalDOFs
        hc_T = sparse.lil_matrix(shape=(hyb_nodofs, hybtet_nodofs), dtype=N.float64)
        hyb_discs = Struct(mtet=self.discs.E.c_imp, brick=self.discs.E.c_exp)
        mt_perm = hyb_discs.mtet.permuter ; b_perm = hyb_discs.brick.permuter
        entity_dofnos = Struct()
        entity_dofnos.mtet = Struct((k,mt_perm.globalEntityPermutationTable(k))
                                    for k in hyb_discs.mtet.geomEntities)
        entity_dofnos.brick = Struct((k,b_perm.globalEntityPermutationTable(k))
                                     for k in hyb_discs.brick.geomEntities)
        submeshes = HybridMesh.make_hyb_submeshes(Struct(tet=hyb_discs.mtet.mesh,
                                                         brick=hyb_discs.brick.mesh,
                                                         on_hbdry=self.on_hbdry))
        subdiscs = HybridMesh.make_hybsubmesh_discs(submeshes, hyb_discs)
        hyb_dofmaps = HybridMesh.HybridSubdim2GlobDofmaps()
        hyb_dofmaps.init_dofmaps(subdiscs, entity_dofnos)

        mtet_subdim2glob_dofmap = hyb_dofmaps.dofmap.mtet
        brick_subdim2glob_dofmap = hyb_dofmaps.dofmap.brick

        hextrimatcher = HybridMesh.HexTriFaceMatcher(submeshes)
        hface_tface_map = hextrimatcher.hex_tris
        
        for hface in subdiscs.brick.elements:
            tsub_facenos = hface_tface_map[hface.index]
            hface_perm = hface.permutation()
            lh_hc_T, twoface_to_tgsubdofs = HybridMesh.make_hface_locmat(
                subdiscs.mtet.elements[tsub_facenos], hface)
            lh_tet_perm = HybridMesh.rem_constr_perm([twoface_to_tgsubdofs[0], N.array(
                [mtet_subdim2glob_dofmap[i] for i in twoface_to_tgsubdofs[1]], N.int32)])
            lh_brick_perm = HybridMesh.rem_constr_perm([hface_perm[0], N.array(
                [brick_subdim2glob_dofmap[i] for i in hface_perm[1]], N.int32)])
            set_global(hc_T, lh_hc_T, lh_brick_perm, lh_tet_perm)
        self.hybrid_transform_mat_debuginfo = Struct(
            entity_dofnos=entity_dofnos, submeshes=submeshes, subdiscs=subdiscs,
            hyb_dofmaps=hyb_dofmaps, hextrimatcher=hextrimatcher)
        return self._finalizeMatrix(hc_T)
        
    @CacheLast.CachedMethod
    def A_cc(self, dtype=N.float64):
        disc_imp = self.discs.E.c_imp ; disc_exp = self.discs.E.c_exp
        hc_T = self.hybrid_transform_mat()
        M_exp = disc_exp.matrix.mass()
        M_imp = disc_imp.matrix.mass()
        S_imp = disc_imp.matrix.stiffness()
        dtdt = self.dt**2
        return hc_T.matmat(
            (M_imp/dtdt + self.implicit_beta*S_imp).matmat(hc_T.T)
            ) + M_exp/dtdt
               
    def M_cc_imp(self):
        disc_imp = self.discs.E.c_imp ;
        hc_T = self.hybrid_transform_mat()
        M_imp = disc_imp.matrix.mass()
        return hc_T.matmat(M_imp.matmat(hc_T.T))

    def S_cc_imp(self):
        disc_imp = self.discs.E.c_imp ;
        hc_T = self.hybrid_transform_mat()
        S_imp = disc_imp.matrix.stiffness()
        return hc_T.matmat(S_imp.matmat(hc_T.T))

    def M_cc_exp(self):
        disc_exp = self.discs.E.c_exp
        M_exp = disc_exp.matrix.mass()
        return M_exp
        
    def S_cc_exp(self):
        disc_exp = self.discs.E.c_exp
        S_exp = disc_exp.matrix.stiffness()
        return S_exp

    @CacheLast.CachedMethod
    def A_bc(self, dtype=N.float64):
        disc_b = self.discs.E.b ; disc_c = self.discs.E.c_imp
        hc_T = self.hybrid_transform_mat()
        M_o = disc_c.matrix.projectionOnto(disc_b)
        S_o = disc_c.D().matrix.projectionOnto(disc_b.D())
        return (M_o/self.dt**2 + self.implicit_beta*S_o).matmat(hc_T.T)

    @CacheLast.CachedMethod
    def B_cc(self, dtype=N.float64):
        disc_imp = self.discs.E.c_imp ; disc_exp = self.discs.E.c_exp
        hc_T = self.hybrid_transform_mat()
        M_exp = disc_exp.matrix.mass()
        S_exp = disc_exp.matrix.stiffness()
        M_imp = disc_imp.matrix.mass()
        S_imp = disc_imp.matrix.stiffness()
        dtdt = self.dt**2
        return hc_T.matmat(
            (2*M_imp/dtdt - (1 - 2*self.implicit_beta)*S_imp).matmat(hc_T.T)
            ) + 2*M_exp/dtdt - S_exp
    
    @CacheLast.CachedMethod
    def B_cd(self):
        S = self.discs.E.d.D().matrix.projectionOnto(self.discs.E.c_exp.D())
        return self._finalizeMatrix(-S)

    @CacheLast.CachedMethod
    def B_bc(self, dtype=N.float64):
        disc_b = self.discs.E.b ; disc_c = self.discs.E.c_imp
        hc_T = self.hybrid_transform_mat()
        M_o = disc_c.matrix.projectionOnto(disc_b)
        S_o = disc_c.D().matrix.projectionOnto(disc_b.D())
        return (2*M_o/self.dt**2 - (1 - 2*self.implicit_beta)*S_o).matmat(hc_T.T)
    
    @CacheLast.CachedMethod
    def C_dc(self):
        mat = self.discs.E.c_exp.matrix.partialExteriorDerivative(
            self.discs.B.d)
        self.discs.E.c_exp.matrix.partialExteriorDerivative.clearCache()
        return mat
