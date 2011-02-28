from __future__ import division

from itertools import izip
import numpy as N
from scipy import sparse
from scipy.sparse.linalg import dsolve

from NewCode.Exceptions import AllZero
from NewCode.Utilities import CacheLast, Struct
from NewCode.MatrixUtils import merge_sparse_blocks, extract_diag_mat
from NewCode.SystemMatrix import local_self_projection_matrix, insert_global
from NewCode.DifferentialForm import BrickDiscretiser, DiscretiserMatrices


class HybridMergedMats(DiscretiserMatrices.Matrices):
    def __init__(self, hybrid_block_mats):
        self.hybrid_block_mats = hybrid_block_mats
        try: self.dt = hybrid_block_mats.dt
        except AttributeError: pass

    def set_dt(self, dt):
        self.clearAll()
        self.dt = dt

    @CacheLast.CachedMethod
    def A_imp(self):
        bmats = Struct((k, getattr(self.hybrid_block_mats, 'A_'+k)())
                     for k in ('aa', 'ab', 'bb', 'bc', 'cc'))
        ofs_b = bmats.aa.shape[0]
        ofs_c = ofs_b + bmats.bb.shape[0]

        return self._finalizeMatrix(merge_sparse_blocks((
            ((0,0), bmats.aa), ((0, ofs_b), bmats.ab), 
            ((ofs_b, 0), bmats.ab.T), ((ofs_b, ofs_b), bmats.bb),
            ((ofs_b, ofs_c), bmats.bc),
            ((ofs_c, ofs_b), bmats.bc.T), ((ofs_c, ofs_c), bmats.cc))))
                                   
    @CacheLast.CachedMethod
    def A_imp_LU(self):
        print "Sparse LU Decomposition of A_imp going DOOOOOOOWWWWWWNN"
        return Struct(solve=dsolve.factorized(self.A_imp()))
    
    @CacheLast.CachedMethod
    def A_exp(self):
        bmats = Struct((k, getattr(self.hybrid_block_mats, 'A_'+k)())
                     for k in ('dd', 'ee'))
        ofs_e = bmats.dd.shape[0]

        return merge_sparse_blocks((
            ((0,0), bmats.dd), ((ofs_e, ofs_e), bmats.ee)))

    @CacheLast.CachedMethod
    def A(self):
        bmats = Struct((k, getattr(self.hybrid_block_mats, 'A_'+k)())
                     for k in ('aa', 'ab', 'bb', 'bc', 'cc', 'dd'))
        ofs_b = bmats.aa.shape[0]
        ofs_c = ofs_b + bmats.bb.shape[0]
        ofs_d = ofs_c + bmats.cc.shape[0]
        ofs_e = ofs_d + bmats.dd.shape[0]

        if self.hybrid_block_mats.discs.E.e.totalDOFs > 0:
            bmats.ee = self.hybrid_block_mats.A_ee
            return self._finalizeMatrix(merge_sparse_blocks((
                ((0,0), bmats.aa), ((0, ofs_b), bmats.ab), 
                ((ofs_b, 0), bmats.ab.T), ((ofs_b, ofs_b), bmats.bb),
                ((ofs_b, ofs_c), bmats.bc),
                ((ofs_c, ofs_b), bmats.bc.T), ((ofs_c, ofs_c), bmats.cc),
                ((ofs_d, ofs_d), bmats.dd), ((ofs_e, ofs_e), bmats.ee))))
        else:
            return self._finalizeMatrix(merge_sparse_blocks((
                ((0,0), bmats.aa), ((0, ofs_b), bmats.ab), 
                ((ofs_b, 0), bmats.ab.T), ((ofs_b, ofs_b), bmats.bb),
                ((ofs_b, ofs_c), bmats.bc),
                ((ofs_c, ofs_b), bmats.bc.T), ((ofs_c, ofs_c), bmats.cc),
                ((ofs_d, ofs_d), bmats.dd))))


    @CacheLast.CachedMethod
    def B(self):
        bmats = Struct((k, getattr(self.hybrid_block_mats, 'B_'+k)())
                     for k in ('aa', 'ab', 'bb', 'bc', 'cc',
                               'cd', 'dd'))
        ofs_b = bmats.aa.shape[0]
        ofs_c = ofs_b + bmats.bb.shape[0]
        ofs_d = ofs_c + bmats.cc.shape[0]
        ofs_e = ofs_d + bmats.dd.shape[0]

        if self.hybrid_block_mats.discs.E.e.totalDOFs > 0:
            bmats.de = self.hybrid_block_mats.B_de
            bmats.ee = self.hybrid_block_mats.B_ee
            return self._finalizeMatrix(merge_sparse_blocks((
                ((0,0), bmats.aa), ((0, ofs_b), bmats.ab), 
                ((ofs_b, 0), bmats.ab.T), ((ofs_b, ofs_b), bmats.bb),
                ((ofs_b, ofs_c), bmats.bc),
                ((ofs_c, ofs_b), bmats.bc.T), ((ofs_c, ofs_c), bmats.cc),
                ((ofs_c, ofs_d), bmats.cd),
                ((ofs_d, ofs_c), bmats.cd.T), ((ofs_d, ofs_d), bmats.dd),
                ((ofs_d, ofs_e), bmats.de),
                ((ofs_e, ofs_d), bmats.de.T), ((ofs_e, ofs_e), bmats.ee)))
                                        )
        else:
            return self._finalizeMatrix(merge_sparse_blocks((
                ((0,0), bmats.aa), ((0, ofs_b), bmats.ab), 
                ((ofs_b, 0), bmats.ab.T), ((ofs_b, ofs_b), bmats.bb),
                ((ofs_b, ofs_c), bmats.bc),
                ((ofs_c, ofs_b), bmats.bc.T), ((ofs_c, ofs_c), bmats.cc),
                ((ofs_c, ofs_d), bmats.cd),
                ((ofs_d, ofs_c), bmats.cd.T), ((ofs_d, ofs_d), bmats.dd)))
                                        )
            
        
        
class HybridBlockMats(DiscretiserMatrices.Matrices):
    no_dt_dep = set(['C_dc', 'C_dd', 'C_ed', 'C_ee', 'P_de', 'P_dd', 'P_ee'])
    def __init__(self, discs, implicit_beta, imp_elset):
        self.discs = discs
        self.implicit_beta = implicit_beta
        self.imp_elset = imp_elset
        
    def set_dt(self, dt):
        for name, ca in self.cachedAttrs.iteritems():
            if name not in self.no_dt_dep: ca.clearCache()
        self.dt = dt

    @CacheLast.CachedMethod
    def A_aa(self):
        M = self.discs.E.a.matrix.mass() 
        S = self.discs.E.a.matrix.stiffness() 
        return M/self.dt**2 + self.implicit_beta*S

    @CacheLast.CachedMethod
    def A_ab(self):
        M = self.discs.E.b.matrix.projectionOnto(self.discs.E.a)
        S = self.discs.E.b.D().matrix.projectionOnto(self.discs.E.a.D())
        return M/self.dt**2 + self.implicit_beta*S

    @CacheLast.CachedMethod
    def A_bb(self):
        M = self.discs.E.b.matrix.mass()
        S = self.discs.E.b.matrix.stiffness() 
        return M/self.dt**2 + self.implicit_beta*S

    @CacheLast.CachedMethod
    def A_bc(self):
        M = self.discs.E.c.matrix.projectionOnto(self.discs.E.b)
        S = self.discs.E.c.D().matrix.projectionOnto(self.discs.E.b.D())
        return M/self.dt**2 + self.implicit_beta*S

    @CacheLast.CachedMethod
    def A_cc(self, dtype=N.float64):
        disc = self.discs.E.c
        no_dofs = disc.totalDOFs
        try: X_cc = sparse.lil_matrix((no_dofs, no_dofs), dtype=dtype)
        except TypeError:
            X_cc = sparse.lil_matrix((1,1), dtype=dtype)

        dtdt = self.dt**2
        for el, el_D in izip(disc.elements, disc.D().elements):
            try: perm = el.permutation()
            except AllZero: continue
            local_mat = local_self_projection_matrix(el)/dtdt
            if el.index in self.imp_elset:
                local_mat += self.implicit_beta*local_self_projection_matrix(el_D)
            insert_global(X_cc, local_mat, el.permutation())
        return X_cc

    @CacheLast.CachedMethod
    def A_dd(self):
        M = extract_diag_mat(self.discs.E.d.matrix.mass())
        return M/self.dt**2

    @CacheLast.CachedMethod
    def A_ee(self):
        M = extract_diag_mat(self.discs.E.e.matrix.mass())
        return M/self.dt**2

    @CacheLast.CachedMethod
    def B_aa(self):
        M = self.discs.E.a.matrix.mass()
        S = self.discs.E.a.matrix.stiffness()
        return self._finalizeMatrix(2*M/self.dt**2 - (1 - 2*self.implicit_beta)*S)

    @CacheLast.CachedMethod
    def B_bb(self):
        M = self.discs.E.b.matrix.mass()
        S = self.discs.E.b.matrix.stiffness() 
        return self._finalizeMatrix(2*M/self.dt**2 - (1 - 2*self.implicit_beta)*S)

    @CacheLast.CachedMethod
    def B_ab(self):
        M = self.discs.E.b.matrix.projectionOnto(self.discs.E.a)
        S = self.discs.E.b.D().matrix.projectionOnto(self.discs.E.a.D())
        return self._finalizeMatrix(2*M/self.dt**2 - (1 - 2*self.implicit_beta)*S)

    @CacheLast.CachedMethod
    def B_bc(self):
        M = self.discs.E.c.matrix.projectionOnto(self.discs.E.b)
        S = self.discs.E.c.D().matrix.projectionOnto(self.discs.E.b.D())
        return self._finalizeMatrix(2*M/self.dt**2 - (1 - 2*self.implicit_beta)*S)

    @CacheLast.CachedMethod
    def B_cc(self, dtype=N.float64):
        disc = self.discs.E.c
        no_dofs = disc.totalDOFs
        try: X_cc = sparse.lil_matrix((no_dofs, no_dofs), dtype=dtype)
        except TypeError:
            X_cc = sparse.lil_matrix((1,1), dtype=dtype)
        dtdt = self.dt**2
        for el, el_D in izip(disc.elements, disc.D().elements):
            try: perm = el.permutation()
            except AllZero: continue
            local_mat = 2*local_self_projection_matrix(el)/dtdt
            local_mat_D = local_self_projection_matrix(el_D)
            if el.index in self.imp_elset:
                local_mat -= (1 - 2*self.implicit_beta)*local_mat_D
            else: local_mat -= local_mat_D
            insert_global(X_cc, local_mat, el.permutation())
        return self._finalizeMatrix(X_cc)

    @CacheLast.CachedMethod
    def B_cd(self):
#        M = extract_diag_mat(self.discs.E.d.matrix.projectionOnto(self.discs.E.c))
        S = self.discs.E.d.D().matrix.projectionOnto(self.discs.E.c.D())
        return self._finalizeMatrix(-S)
            #2*M/self.dt**2 - S)

    @CacheLast.CachedMethod
    def B_dd(self):
        M = extract_diag_mat(self.discs.E.d.matrix.mass() )
        S = self.discs.E.d.matrix.stiffness()
        return self._finalizeMatrix(2*M/self.dt**2 - S)

    @CacheLast.CachedMethod
    def B_de(self):
#        M = self.discs.E.e.matrix.projectionOnto(self.discs.E.d)
        S = self.discs.E.e.D().matrix.projectionOnto(self.discs.E.d.D())
        return self._finalizeMatrix(-S)
            #2*M/self.dt**2 - S)

    @CacheLast.CachedMethod
    def B_ee(self):
        M = extract_diag_mat(self.discs.E.e.matrix.mass() )
        S = self.discs.E.e.matrix.stiffness()
        return self._finalizeMatrix(2*M/self.dt**2 - S)

    @CacheLast.CachedMethod
    def C_dc(self):
        mat = self.discs.E.c.matrix.partialExteriorDerivative(
            self.discs.B.d)
        self.discs.E.c.matrix.partialExteriorDerivative.clearCache()
        return mat

    @CacheLast.CachedMethod
    def C_dd(self):
        mat = self.discs.E.d.matrix.partialExteriorDerivative(
            self.discs.B.d)
        self.discs.E.d.matrix.partialExteriorDerivative.clearCache()
        return mat

    @CacheLast.CachedMethod
    def C_ed(self):
        mat = self.discs.E.d.matrix.partialExteriorDerivative(
            self.discs.B.e)
        self.discs.E.d.matrix.partialExteriorDerivative.clearCache()
        return mat

    @CacheLast.CachedMethod
    def C_ee(self):
        mat = self.discs.E.e.matrix.exteriorDerivative(
            self.discs.B.e)
        self.discs.E.e.matrix.exteriorDerivative.clearCache()
        return mat

    @CacheLast.CachedMethod
    def P_dd(self):
        mat = self.discs.B.d.matrix.projectionOnto(
            self.discs.E.d.D())
        self.discs.B.d.matrix.projectionOnto.clearCache()
        return mat

    @CacheLast.CachedMethod
    def P_de(self):
        mat = self.discs.B.e.matrix.projectionOnto(
            self.discs.E.d.D())
        self.discs.B.e.matrix.projectionOnto.clearCache()
        return mat

    @CacheLast.CachedMethod
    def P_ee(self):
        mat = self.discs.B.e.matrix.projectionOnto(
            self.discs.E.e.D())
        self.discs.B.e.matrix.projectionOnto.clearCache()
        return mat

    def B_allblocks(self):
        diag_blocks = ('aa', 'bb', 'cc', 'dd', 'ee')
        upper_blocks = ('ab', 'bc', 'cd', 'de')
        lower_blocks = tuple(x[::-1] for x in upper_blocks)
        allblocks = Struct((b, getattr(self, 'B_'+b)())
                           for b in diag_blocks + upper_blocks)
        allblocks.update(dict((b, allblocks[b[::-1]].T)
                              for b in lower_blocks))
        return allblocks

    def B_impblocks(self):
        diag_blocks = ('aa', 'bb', 'cc')
        upper_blocks = ('ab', 'bc', 'cd')
        lower_blocks = ('ba', 'cb')
        impblocks = Struct((b, getattr(self, 'B_'+b)())
                           for b in diag_blocks + upper_blocks)
        impblocks.update(dict((b, impblocks[b[::-1]].T)
                              for b in lower_blocks))
        return impblocks
