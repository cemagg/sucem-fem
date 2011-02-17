from __future__ import division

import numpy as N

from NewCode.Utilities import partial, Struct
from BrickCommon import extended_chebyshev_interp_pts, cov_coord_vecs, con_coord_vecs
import BrickCommon

def basis_set(order, mixed=True, btype='rieben04'):
    if btype == 'rieben04': bset = OneformBasisSet_rieben04(order, mixed)
    elif btype == 'cohen98': bset = OneformBasisSet_cohen98(order, mixed)
    else: raise ValueError('Unknown Brick basis type: '+ str(btype))
    return bset

class InterpBasisFuncs(BrickCommon.InterpBasisFuncs):
    @staticmethod
    def edgeFun(fn, tfn, local_edgeno, facedual_edge_numbering, name='edgefun_'):
        # The edge tangent direction is normal to the face normals of the two faces
        # that define the edge. Taking mod 3 of the face numbers gives only the
        # independent coords (i.e. lamda1 - lamda3).
        edge_faces = facedual_edge_numbering[local_edgeno]
        edge_direction_ind = (set([0,1,2]) - set(tuple(N.mod(edge_faces, 3)))).pop()
        edge_direction = cov_coord_vecs[edge_direction_ind]
        # If the edge is defined by e.g. faces 2 and 3, the function is zero'd
        # along the other edges by multiplying by the opposing faces'
        # zero-coordinate surface. This is using Graglia97's face/coordinate
        # numbering scheme.
        local_coords = edge_faces.tolist() + [edge_direction_ind]
        def edgefun_(coord):
            try:
                coords = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn(coords[0])*fn(coords[1])*tfn(coords[2])*edge_direction
        edgefun_.__name__ = name
        edgefun_.direction = edge_direction
        edgefun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fn, fn, tfn])
        return edgefun_


    @staticmethod
    def edgeFun_D(fn, tfn, local_edgeno, facedual_edge_numbering, name='edgefun_D_'):
        edge_faces = facedual_edge_numbering[local_edgeno]
        edge_direction_ind = (set([0,1,2]) - set(tuple(N.mod(edge_faces, 3)))).pop()
        con_coords = [con_coord_vecs[edge_face, edge_direction_ind]
                      for edge_face in edge_faces]
        local_coords = edge_faces.tolist() + [edge_direction_ind]
        fn_D = N.polyder(fn)
        def edgefun_D_(coord):
            try:
                coords = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn(coords[0])*fn_D(coords[1])*tfn(coords[2])*con_coords[1] +\
                   fn_D(coords[0])*fn(coords[1])*tfn(coords[2])*con_coords[0]    
        edgefun_D_.__name__ = name+'_D'
        return edgefun_D_

    @staticmethod
    def faceFuna(fn0, fn, tfn, local_faceno, dummy, name='facefuna_'):
        local_coords = [local_faceno, (local_faceno+1) % 3, (local_faceno+2) % 3]
        cov_a = cov_coord_vecs[local_coords[1]]
        def facefuna_(coord):
            try:
                coords = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn0(coords[0])*tfn(coords[1])*fn(coords[2])*cov_a
        
        facefuna_.__name__ = name+'a'
        facefuna_.direction = cov_a
        facefuna_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fn0, tfn, fn])
        return facefuna_

    @staticmethod
    def faceFunb(fn0, fn, tfn, local_faceno, dummy, name='facefunb_'):
        local_coords = [local_faceno, (local_faceno+1) % 3, (local_faceno+2) % 3]
        cov_b = cov_coord_vecs[local_coords[2]]
        def facefunb_(coord):
            try:
                coords = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn0(coords[0])*fn(coords[1])*tfn(coords[2])*cov_b

        facefunb_.__name__ = name+'b'
        facefunb_.direction = cov_b
        facefunb_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fn0, fn, tfn])
        return facefunb_

    @staticmethod
    def faceFuna_D(fn0, fn, tfn, local_faceno, dummy, name='facefuna_D_'):
        local_coords = [local_faceno, (local_faceno+1) % 3, (local_faceno+2) % 3]
        lc0, lc1, lc2 = local_coords
        con0_1 = con_coord_vecs[lc0, lc1]
        con2_1 = con_coord_vecs[lc2, lc1]
        fn0_D = N.polyder(fn0)
        fn_D = N.polyder(fn)
        def facefuna_D_(coord):
            try:
                c0,c1,c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn0_D(c0)*fn(c2)*tfn(c1)*con0_1 + \
                   fn0(c0)*fn_D(c2)*tfn(c1)*con2_1 
        
        facefuna_D_.__name__ = name+'a_D'
        return facefuna_D_

    @staticmethod
    def faceFunb_D(fn0, fn, tfn, local_faceno, dummy, name='facefunb_D_'):
        local_coords = [local_faceno, (local_faceno+1) % 3, (local_faceno+2) % 3]
        lc0, lc1, lc2 = local_coords
        con0_2 = con_coord_vecs[lc0, lc2]
        con1_2 = con_coord_vecs[lc1, lc2]
        fn0_D = N.polyder(fn0)
        fn_D = N.polyder(fn)
        def facefunb_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn0_D(c0)*fn(c1)*tfn(c2)*con0_2 + \
                   fn0(c0)*fn_D(c1)*tfn(c2)*con1_2 

        facefunb_D_.__name__ = name+'b_D'
        return facefunb_D_

    @staticmethod
    def volFun_a(fni, fnj, tfn, name='volfun'):
        lc0, lc1, lc2 = local_coords = [1, 2, 0]
        cov = cov_coord_vecs[lc2]
        def volfun_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni(c0)*fnj(c1)*tfn(c2)*cov
        volfun_.__name__ = name+'_a'
        volfun_.direction = cov
        volfun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fni, fnj, tfn])
        return volfun_

    @staticmethod
    def volFun_a_D(fni, fnj, tfn, name='volfun_D'):
        lc0, lc1, lc2 = local_coords = [1, 2, 0]
        con0_2 = con_coord_vecs[lc0, lc2]
        con1_2 = con_coord_vecs[lc1, lc2]
        fni_D = N.polyder(fni)
        fnj_D = N.polyder(fnj)
        def volfun_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni_D(c0)*fnj(c1)*tfn(c2)*con0_2 + \
                   fni(c0)*fnj_D(c1)*tfn(c2)*con1_2 
        volfun_D_.__name__ = name + '_a_D'
        return volfun_D_

    @staticmethod
    def volFun_b(fni, fnj, tfn, name='volfun'):
        lc0, lc1, lc2 = local_coords = [0, 2, 1]
        cov = cov_coord_vecs[lc2]
        def volfun_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni(c0)*fnj(c1)*tfn(c2)*cov
        volfun_.__name__ = name + '_b'
        volfun_.direction = cov
        volfun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fni, fnj, tfn])
        return volfun_

    @staticmethod
    def volFun_b_D(fni, fnj, tfn, name='volfun_D_'):
        lc0, lc1, lc2 = local_coords = [0, 2, 1]
        con0_2 = con_coord_vecs[lc0, lc2]
        con1_2 = con_coord_vecs[lc1, lc2]
        fni_D = N.polyder(fni)
        fnj_D = N.polyder(fnj)
        def volfun_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni_D(c0)*fnj(c1)*tfn(c2)*con0_2 + \
                   fni(c0)*fnj_D(c1)*tfn(c2)*con1_2 
        volfun_D_.__name__ = name + '_b_D'
        return volfun_D_

    @staticmethod
    def volFun_c(fni, fnj, tfn, name='volfun_c_'):
        lc0, lc1, lc2 = local_coords = [0, 1, 2]
        cov = cov_coord_vecs[lc2]
        def volfun_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni(c0)*fnj(c1)*tfn(c2)*cov
        volfun_.__name__ = name + '_c'
        volfun_.direction = cov
        volfun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fni, fnj, tfn])
        return volfun_

    @staticmethod
    def volFun_c_D(fni, fnj, tfn, name='volfun_D_'):
        lc0, lc1, lc2 = local_coords = [0, 1, 2]
        con0_2 = con_coord_vecs[lc0, lc2]
        con1_2 = con_coord_vecs[lc1, lc2]
        fni_D = N.polyder(fni)
        fnj_D = N.polyder(fnj)
        def volfun_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni_D(c0)*fnj(c1)*tfn(c2)*con0_2 + \
                   fni(c0)*fnj_D(c1)*tfn(c2)*con1_2 
        volfun_D_.__name__ = name + '_c_D'
        return volfun_D_

    
    def _get_edgeFuncs(self, efgen):
        fn = self.fInterpFns[0]; tfns = self.rInterpFns
        return [partial(efgen, fn, tfn, name='edgefun_%f_%d' %
                        (self.order, i))
                for i, tfn in enumerate(tfns)]

    def get_edgeFuncs(self):
        return self._get_edgeFuncs(self.edgeFun)

    def get_edgeFuncs_D(self):
        return self._get_edgeFuncs(self.edgeFun_D)

    def _get_faceFuncs(self, ffgen):
        fn0 = self.fInterpFns[0] ; fns = self.fInterpFns[1:-1]
        tfns = self.rInterpFns
        ffs = []
        for j, fn in enumerate(fns):
            for k, tfn in enumerate(tfns):
                ffs.append(partial(ffgen, fn0, fn, tfn,
                                   name='facefun_%f_%d_%d_' % (self.order, j+1, k)))
        return ffs

    def _get_volFuncs(self, vfgen):
        fns = self.fInterpFns[1:-1]
        tfns = self.rInterpFns
        vfs = []
        for i, fni in enumerate(fns):
            for j, fnj in enumerate(fns):
                for k, tfn in enumerate(tfns):
                    vfs.append(vfgen(fni, fnj, tfn, name='volfun_%f_%d_%d_%d_' %
                                       (self.order, i+1,j+1, k)))
        return vfs
        
    def get_faceFuncs(self):
        return self._get_faceFuncs(self.faceFuna) + self._get_faceFuncs(self.faceFunb)

    def get_faceFuncs_D(self):
        return self._get_faceFuncs(self.faceFuna_D) + self._get_faceFuncs(self.faceFunb_D)

    def get_volFuncs(self):
        return self._get_volFuncs(self.volFun_a) + \
               self._get_volFuncs(self.volFun_b) + \
               self._get_volFuncs(self.volFun_c)

    def get_volFuncs_D(self):
        return self._get_volFuncs(self.volFun_a_D) + \
               self._get_volFuncs(self.volFun_b_D) + \
               self._get_volFuncs(self.volFun_c_D)
            
class OneformBrickBasisSet(object):
    form = 1
    form_dof_order = Struct({'node': None, 'edge': 1, 'face': 2, 'vol': 2})
    InterpBasisFuncs = InterpBasisFuncs
    
class OneformBasisSet_rieben04(OneformBrickBasisSet,
                               BrickCommon.BasisSet_rieben04): pass
    
class OneformBasisSet_cohen98(OneformBrickBasisSet,
                              BrickCommon.BasisSet_cohen98): pass



ib1m = InterpBasisFuncs(extended_chebyshev_interp_pts, 1, blah=True)
edgefun0 = ib1m.get_edgeFuncs()[0]
edgefun0_D = ib1m.get_edgeFuncs_D()[0]
edgefuns0 = tuple(partial(edgefun0, i) for i in range(12))
edgefuns0_D = tuple(partial(edgefun0_D, i) for i in range(12))
