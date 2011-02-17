from __future__ import division
import numpy as N

from NewCode.Utilities import partial, Struct, gen_lagrange_polys

from BrickCommon import face_con_coord_vecs, vol_sign, con_coord_vecs
import BrickCommon

def facefun0(local_faceno, dummy):
    local_coord = (local_faceno + 3) % 6
    face_vec = face_con_coord_vecs[local_faceno % 3]
    def facefun0_(coord):
        try:
            coord = coord[local_coord]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return coord*face_vec
    return facefun0_

def facefun0_D(local_faceno, dummy):
    def facefun0_D_(coord):
        return N.array([-1. if local_faceno < 3 else 1.], N.float64)
    return facefun0_D_

facefuns0 = tuple(partial(facefun0,i) for i in range(6))
facefuns0_D = tuple(partial(facefun0_D,i) for i in range(6))

def basis_set(order, mixed=True, btype='rieben04'):
    if btype == 'rieben04': bset = TwoformBasisSet_rieben04(order, mixed)
    elif btype == 'cohen98': bset = TwoformBasisSet_cohen98(order, mixed)
    else: raise ValueError('Unknown Brick basis type: '+ btype)
    return bset

class InterpBasisFuncs(BrickCommon.InterpBasisFuncs):

    @staticmethod
    def faceFun(fn0, rfnj, rfnk, local_faceno, dummy, name='facefun_'):
        local_coords = [local_faceno, (local_faceno-2) % 3, (local_faceno-1) % 3]
        lc0,lc1,lc2 = local_coords
        conv = con_coord_vecs[lc1, lc2]
        def facefun_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fn0(c0)*rfnj(c1)*rfnk(c2)*conv
        
        facefun_.__name__ = name
        facefun_.direction = conv
        facefun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fn0, rfnj, rfnk])
        return facefun_

    @staticmethod
    def faceFun_D(fn0, rfnj, rfnk, local_faceno, dummy, name='facefun_D_'):
        local_coords = [local_faceno, (local_faceno-2) % 3, (local_faceno-1) % 3]
        lc0,lc1,lc2 = local_coords
        sign = vol_sign[lc0, lc1, lc2]
        fn0_d = N.polyder(fn0)
        def facefun_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return N.array([fn0_d(c0)*rfnj(c1)*rfnk(c2)*sign], N.float64)
        
        facefun_D_.__name__ = name
        return facefun_D_

    @staticmethod
    def volFun(fni, rfnj, rfnk, direction, name='volfun_'):
        lc0,lc1,lc2 = local_coords = [direction, (direction-2)%3, (direction-1)%3]
        conv = con_coord_vecs[lc1, lc2]
        def volfun_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return fni(c0)*rfnj(c1)*rfnk(c2)*conv
        volfun_.__name__ = name
        volfun_.direction = conv
        volfun_.interp_pt = BrickCommon.calc_interp_pt(
            local_coords, [fni, rfnj, rfnk])
        return volfun_

    @staticmethod
    def volFun_D(fni, rfnj, rfnk, direction, name='volfun_D_'):
        lc0,lc1,lc2 = local_coords = [direction, (direction-2)%3, (direction-1)%3]
        sign = vol_sign[lc0, lc1, lc2]
        fni_d = N.polyder(fni)
        def volfun_D_(coord):
            try:
                c0, c1, c2 = coord[local_coords]
            except TypeError:
                raise TypeError('coord must be a numpy array')
            return N.array([fni_d(c0)*rfnj(c1)*rfnk(c2)*sign], N.float64)
        volfun_D_.__name__ = name
        return volfun_D_

    def _get_faceFuncs(self, ffgen):
        fn0 = self.fInterpFns[0] ; rfns = self.rInterpFns
        ffs = []
        for j, rfnj in enumerate(rfns):
            for k, rfnk in enumerate(rfns):
                ffs.append(partial(ffgen, fn0, rfnj, rfnk,
                                   name='facefun_%f_%d_%d_' % (self.order, j, k)))
        return ffs

    def _get_volFuncs(self, vfgen):
        fns = self.fInterpFns[1:-1]
        rfns = self.rInterpFns
        vfs = []
        spacedim=3
        for i, fni in enumerate(fns):
            for j, rfnj in enumerate(rfns):
                for k, rfnk in enumerate(rfns):
                    vfs.extend(vfgen(fni, rfnj, rfnk,dr, name='volfun_%f_%d_%d_%d_' %
                                     (self.order, i+1,j, k))
                               for dr in range(spacedim))
        return vfs

    def get_faceFuncs(self):
        return self._get_faceFuncs(self.faceFun)

    def get_faceFuncs_D(self):
        return self._get_faceFuncs(self.faceFun_D)

    def get_volFuncs(self):
        return self._get_volFuncs(self.volFun)

    def get_volFuncs_D(self):
        return self._get_volFuncs(self.volFun_D)

class TwoformBrickBasisSet(object):
    form = 2
    form_dof_order = Struct({'node': None, 'edge': None, 'face': 1, 'vol': 2})
    InterpBasisFuncs = InterpBasisFuncs

class TwoformBasisSet_rieben04(TwoformBrickBasisSet,
                               BrickCommon.BasisSet_rieben04): pass

class TwoformBasisSet_cohen98(TwoformBrickBasisSet,
                              BrickCommon.BasisSet_cohen98): pass
