from __future__ import division

import numpy as N

from NewCode.Utilities import partial, Struct
from itertools import chain
import PyramCommon
from PyramCommon import con_coord_vecs, cov_coord_vecs, shapefuns, D_shapefuns
from NewCode.Meshes.PyramMesh import Element

def basis_set(order, mixed=True, btype='graglia99'):
    if btype == 'graglia99': bset = OneformBasisSet_graglia99(order, mixed)
    elif btype == 'coulomb97': bset = OneformBasisSet_coulomb97(order, mixed)
    else: raise ValueError('Unknown Pyramid basis type: '+ str(btype))
    return bset



class OneformBasisSet_Pyr(Struct):
    no_edges = 8
    no_base_edges = 4
    no_apex_edges = 4
    no_edges = no_base_edges + no_apex_edges
    no_faces = 5
    no_apexfaces = 4
    no_base_faces = 1
    form_dof_order = Struct(node=None, edge=1, apexface=2,
                            baseface=2, basevol=2, apexvol=3)
    form = 1
    def __init__(self, order, mixed=True):
        self.fns = Struct()
        self.fns_D = Struct()
        fdo = self.form_dof_order
        if order >= fdo.edge:
            self.fns.edge, self.fns_D.edge = self._get_edgefuns(order,mixed)
        if order >= fdo.apexface:
            self.fns.apexface, self.fns_D.apexface = self._get_apexfacefuns(order,mixed)
        if order >= fdo.baseface:
            self.fns.baseface, self.fns_D.baseface = self._get_basefacefuns(order,mixed)
        if order >= fdo.basevol:
            self.fns.vol, self.fns_D.vol = self._get_basevolfuns(order,mixed)
        if order >= fdo.apexvol:
            apexvol, apexvol_D = self._get_apexvolfuns(order,mixed)
#            self.fns.vol, self.fns_D.vol = apexvol, apexvol_D
#            self.fns.vol.append(apexvol), 
        if mixed: order -= 0.5
        self.info = Struct(form=self.form, type=self.btype, order=order )           

class OneformBasisSet_coulomb97(OneformBasisSet_Pyr):
    btype = 'coulomb97'
    def __init__(self, order, mixed=True):
        assert(order < 2 or (order <=2 and mixed))
        super(OneformBasisSet_coulomb97, self).__init__(order, mixed)
    
    def _get_edgefuns(self, order, mixed):
        # We assume the first four local edges are base edges and the last four
        # are apex-edges
        no_base, no_apex = self.no_base_edges, self.no_apex_edges
        edgefuns, edgefuns_D = [], []
        for i in range(no_base):
            edgefuns.append(partial(self.edgefun0_base, i))
            edgefuns_D.append(partial(self.edgefun0_base_D, i))
        for i in range(no_base, no_base+no_apex):
            edgefuns.append(partial(self.edgefun0_apex, i))
            edgefuns_D.append(partial(self.edgefun0_apex_D, i))

        if not mixed or order > 1:
            for i in range(no_base):
                edgefuns.append(partial(self.edgefun1_base, i))
                edgefuns_D.append(partial(self.edgefun1_base_D, i))
            for i in range(no_base, no_base+no_apex):
                edgefuns.append(partial(self.edgefun1_apex, i))
                edgefuns_D.append(partial(self.edgefun1_apex_D, i))

        return edgefuns, edgefuns_D

    def _get_apexfacefuns(self, order, mixed):
        facefuns, facefuns_D = [], []
        for i in range(self.no_apexfaces):
            facefuns.append(partial(self.facefun2m_apex1, i))
            facefuns_D.append(partial(self.facefun2m_apex1_D, i))
        for i in range(self.no_apexfaces):
            facefuns.append(partial(self.facefun2m_apex2, i))
            facefuns_D.append(partial(self.facefun2m_apex2_D, i))
        return facefuns, facefuns_D

    def _get_basefacefuns(self, order, mixed):
        facefuns, facefuns_D = [], []
        for i in range(self.no_base_edges):
            facefuns.append(partial(self.facefun2m_base, i))
            facefuns_D.append(partial(self.facefun2m_base_D, i))
        return facefuns, facefuns_D

    def _get_basevolfuns(self, order, mixed):
        return [self.volfun2], [self.volfun2_D]
        
            
    @staticmethod
    def edgefun0_base(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(fi != 4)                      # fi must be a triangle face
        assert(fj == 4)                      # fj must be the base face
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        sfi, sfj = shapefuns[[i,j]]
        D_sfi, D_sfj, D_sfi_do, D_sfj_do = D_shapefuns[[i,j,i_do,j_do]]
        def _edgefun0_base(coord):
            # coulomb97 table 1, base edges i-j
            return sfi(coord)*(D_sfj(coord)+D_sfi_do(coord)) \
                   - sfj(coord)*(D_sfi(coord)+D_sfj_do(coord))
        return _edgefun0_base

    @staticmethod
    def edgefun0_base_D(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(fi != 4)                      # fi must be a triangle face
        assert(fj == 4)                      # fj must be the base face
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        D_sfi, D_sfj, D_sfi_do, D_sfj_do = D_shapefuns[[i,j,i_do,j_do]]
        def _edgefun0_base_D(coord):
            # coulomb97 table 1, base edges i-j
            return N.cross(D_sfi(coord),(D_sfj(coord)+D_sfi_do(coord))) \
                   - N.cross(D_sfj(coord),(D_sfi(coord)+D_sfj_do(coord)))
        return _edgefun0_base_D

    @staticmethod
    def edgefun1_base(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(fi != 4)                      # fi must be a triangle face
        assert(fj == 4)                      # fj must be the base face
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        sfi, sfj = shapefuns[[i,j]]
        D_sfi, D_sfj, D_sfi_do, D_sfj_do = D_shapefuns[[i,j,i_do,j_do]]
        def _edgefun1_base(coord):
            # coulomb97 table 1, base edges i-j
            return sfi(coord)*(D_sfj(coord)+D_sfi_do(coord)) \
                   + sfj(coord)*(D_sfi(coord)+D_sfj_do(coord))
        return _edgefun1_base

    @staticmethod
    def edgefun1_base_D(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(fi != 4)                      # fi must be a triangle face
        assert(fj == 4)                      # fj must be the base face
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        D_sfi, D_sfj, D_sfi_do, D_sfj_do = D_shapefuns[[i,j,i_do,j_do]]
        def _edgefun1_base_D(coord):
            # coulomb97 table 1, base edges i-j
            return N.cross(D_sfi(coord),(D_sfj(coord)+D_sfi_do(coord))) \
                   + N.cross(D_sfj(coord),(D_sfi(coord)+D_sfj_do(coord)))
        return _edgefun1_base_D

    
    @staticmethod
    def edgefun0_apex(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(fj == (fi+1) % 4)
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        assert(j == 4)                  # Apex node
        sfi, sfj = shapefuns[[i,j]]
        D_sfi, D_sfj = D_shapefuns[[i,j]]
        def _edgefun0_apex(coord):
            # coulomb97 table 1, apex edges i-5
            return sfi(coord)*D_sfj(coord) - sfj(coord)*D_sfi(coord)
        return _edgefun0_apex

    @staticmethod
    def edgefun0_apex_D(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(fj == (fi+1) % 4)
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        assert(j == 4)                  # Apex node
        D_sfi, D_sfj = D_shapefuns[[i,j]]
        def _edgefun0_apex_D(coord):
            # coulomb97 table 1, apex edges i-5
            return N.cross(D_sfi(coord),D_sfj(coord)) \
                   - N.cross(D_sfj(coord),D_sfi(coord))
        return _edgefun0_apex_D

    @staticmethod
    def edgefun1_apex(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(fj == (fi+1) % 4)
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        assert(j == 4)                  # Apex node
        sfi, sfj = shapefuns[[i,j]]
        D_sfi, D_sfj = D_shapefuns[[i,j]]
        def _edgefun1_apex(coord):
            # coulomb97 table 1, apex edges i-5
            return sfi(coord)*D_sfj(coord) + sfj(coord)*D_sfi(coord)
        return _edgefun1_apex

    @staticmethod
    def edgefun1_apex_D(local_edgeno, edgefaces, edge_sense):
        fi,fj = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(fj == (fi+1) % 4)
        i,j = Element.LOCAL_EDGENODES[local_edgeno]
        assert(j == 4)                  # Apex node
        D_sfi, D_sfj = D_shapefuns[[i,j]]
        def _edgefun1_apex_D(coord):
            # coulomb97 table 1, apex edges i-5
            return N.cross(D_sfi(coord),D_sfj(coord)) \
                   + N.cross(D_sfj(coord),D_sfi(coord))
        return _edgefun1_apex_D

    @staticmethod
    def facefun2m_base(baseface_edgeno, *names):
        i,j = Element.LOCAL_EDGENODES[
            Element.LOCAL_BASEFACEEDGES[baseface_edgeno]]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        sfi, sfj = shapefuns[[i,j]]
        D_sfi_do, D_sfj_do = D_shapefuns[[i_do,j_do]]
        def _facefun2m_base(coord):
            # coulomb97 table 1, base face fns
            return sfi(coord)*sfj(coord)*(D_sfi_do(coord)+D_sfj_do(coord))
        return _facefun2m_base
        
    @staticmethod
    def facefun2m_base_D(baseface_edgeno, *names):
        i,j = Element.LOCAL_EDGENODES[
            Element.LOCAL_BASEFACEEDGES[baseface_edgeno]]
        i_do, j_do = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[[i,j]]
        sfi, sfj = shapefuns[[i,j]]
        D_sfi, D_sfj, D_sfi_do, D_sfj_do = D_shapefuns[[i,j,i_do,j_do]]
        def _facefun2m_base_D(c):
            # coulomb97 table 1, base face fns
            return N.cross(sfi(c)*D_sfj(c) + sfj(c)*D_sfi(c),
                           D_sfi_do(c) + D_sfj_do(c))
        return _facefun2m_base_D

    @staticmethod
    def facefun2m_apex1(local_face, *names):
        i,j = Element.LOCAL_APEXFACENODES[local_face][0:2]
        k = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[i]
        sfi, sfb = shapefuns[[i,4]]
        D_sfj, D_sfk = D_shapefuns[[j,k]]
        def _facefun2m_apex1(c):
            return sfi(c)*sfb(c)*(D_sfj(c) + D_sfk(c))
        return _facefun2m_apex1

    @staticmethod
    def facefun2m_apex1_D(local_face, *names):
        i,j = Element.LOCAL_APEXFACENODES[local_face][0:2]
        k = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[i]
        sfi, sfb = shapefuns[[i,4]]
        D_sfi, D_sfb, D_sfj, D_sfk = D_shapefuns[[i,4,j,k]]
        def _facefun2m_apex1_D(c):
            return N.cross(sfi(c)*D_sfb(c)+sfb(c)*D_sfi(c),
                           D_sfj(c) + D_sfk(c))
        return _facefun2m_apex1_D

    @staticmethod
    def facefun2m_apex2(local_face, *names):
        i,j = Element.LOCAL_APEXFACENODES[local_face][0:2]
        l = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[j]
        sfj, sfb = shapefuns[[j,4]]
        D_sfi, D_sfl = D_shapefuns[[i,l]]
        def _facefun2m_apex2(c):
            return sfj(c)*sfb(c)*(D_sfi(c) + D_sfl(c))
        return _facefun2m_apex2

    @staticmethod
    def facefun2m_apex2_D(local_face, *names):
        i,j = Element.LOCAL_APEXFACENODES[local_face][0:2]
        l = Element.LOCAL_BASEFACENODES_DIAG_OPPOSED[j]
        sfj, sfb = shapefuns[[j,4]]
        D_sfi, D_sfb, D_sfj, D_sfl = D_shapefuns[[i,4,j,l]]
        def _facefun2m_apex2_D(c):
            return N.cross(sfj(c)*D_sfb(c)+sfb(c)*D_sfj(c),
                           D_sfi(c) + D_sfl(c))
        return _facefun2m_apex2_D

    @staticmethod
    def volfun2(*names):
        i,k,b = 0,3,4
        sfi, sfk = shapefuns[[i,k]]
        D_sfb = D_shapefuns[b]
        def _volfun2(c):
            return sfi(c)*sfk(c)*D_sfb(c)
        return _volfun2
    
    @staticmethod
    def volfun2_D(*names):
        i,k,b = 0,3,4
        sfi, sfk = shapefuns[[i,k]]
        D_sfi, D_sfk, D_sfb = D_shapefuns[[i,k,b]]
        def _volfun2_D(c):
            return N.cross(sfi(c)*D_sfk(c) + sfk(c)*D_sfi(c),
                           D_sfb(c))
        return _volfun2_D
    
class OneformBasisSet_graglia99(OneformBasisSet_Pyr):
    btype = 'graglia99'
    # Base edges that have reversed interpolation-point ordering
    reversed_base_edges = (0,2)              
    def __init__(self, order, mixed=True):
        Struct.__init__(self)
        assert(order > 0)
        if order > 2 or not mixed: raise NotImplementedError 
        self.interpFuns = InterpFuncs(order-1)
        super(OneformBasisSet_graglia99, self).__init__(order, mixed=mixed)
        
    @staticmethod
    def combine_edgefun(edgefn, interpfn, local_edgeno, edgefaces, edge_sense):
        _ef = edgefn(local_edgeno, edgefaces, edge_sense)
        _if = interpfn(local_edgeno, edgefaces)
        def _interp_edgefun(coord):
            return _ef(coord)*_if(coord)
        _interp_edgefun.__name__ = _ef.__name__ + _if.__name__
        return _interp_edgefun

    @staticmethod
    def combine_edgefun_D(
        edgefn, D_edgefn, interpfn, local_edgeno, edgefaces, edge_sense):
        _ef = edgefn(local_edgeno, edgefaces, edge_sense)
        _D_ef = D_edgefn(local_edgeno, edgefaces, edge_sense)
        _if = interpfn(local_edgeno, edgefaces)
        _D_if = _if.D
        def _interp_edgefun_D(coord):
            return _if(coord)*_D_ef(coord) + N.cross(_D_if(coord), _ef(coord))
        _interp_edgefun_D.__name__ = _D_ef.__name__ + _D_if.__name__
        return _interp_edgefun_D
        
    @staticmethod
    def combine_apexfacefun(edgefn, interpfn, local_face,
                             faceedges, edgefaces, edge_sense):
        _if = interpfn(local_face)
        local_edgeno = faceedges[local_face, _if.local_edge]
        _ef = edgefn(local_edgeno, edgefaces, edge_sense)
        def _interp_apexfacefun(coord):
            return _ef(coord)*_if(coord)
        _interp_apexfacefun.__name__ = _ef.__name__ + _if.__name__
        return _interp_apexfacefun

    @staticmethod
    def combine_apexfacefun_D(edgefn, D_edgefn, interpfn, local_face,
                             faceedges, edgefaces, edge_sense):
        _if = interpfn(local_face)
        _D_if = _if.D
        local_edgeno = faceedges[local_face, _if.local_edge]
        _ef = edgefn(local_edgeno, edgefaces, edge_sense)
        _D_ef = D_edgefn(local_edgeno, edgefaces, edge_sense)
        def _interp_apexfacefun_D(coord):
            return _if(coord)*_D_ef(coord) + N.cross(_D_if(coord), _ef(coord))
        _interp_apexfacefun_D.__name__ = _D_ef.__name__ + _D_if.__name__
        return _interp_apexfacefun_D

    @staticmethod
    def combine_basefacefun(edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn(edgefaces[basedir][0])
        _ef = edgefn(basedir, edgefaces, edge_sense)
        def _interp_basefacefun(coord):
            return _ef(coord)*_if(coord)
        _interp_basefacefun.__name__ = _ef.__name__ + _if.__name__
        return _interp_basefacefun

    @staticmethod
    def combine_basefacefun_D(edgefn, D_edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn(edgefaces[basedir][0])
        _D_if = _if.D
        _ef = edgefn(basedir, edgefaces, edge_sense)
        _D_ef = D_edgefn(basedir, edgefaces, edge_sense)
        def _interp_basefacefun_D(coord):
            return _if(coord)*_D_ef(coord) + N.cross(_D_if(coord), _ef(coord))
        _interp_basefacefun_D.__name__ = _D_ef.__name__ + _D_if.__name__
        return _interp_basefacefun_D

    @staticmethod
    def combine_basevolfun(edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn(edgefaces[basedir][0])
        _ef = edgefn(basedir, edgefaces, edge_sense)
        def _interp_basevolfun(coord):
            return _ef(coord)*_if(coord)
        _interp_basevolfun.__name__ = _ef.__name__ + _if.__name__
        return _interp_basevolfun

    @staticmethod
    def combine_basevolfun_D(edgefn, D_edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn(edgefaces[basedir][0])
        _D_if = _if.D
        _ef = edgefn(basedir, edgefaces, edge_sense)
        _D_ef = D_edgefn(basedir, edgefaces, edge_sense)
        def _interp_basevolfun_D(coord):
            return _if(coord)*_D_ef(coord) + N.cross(_D_if(coord), _ef(coord))
        _interp_basevolfun_D.__name__ = _D_ef.__name__ + _D_if.__name__
        return _interp_basevolfun_D

    @staticmethod
    def combine_apexvolfun(edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn()
        _ef = edgefn(basedir, edgefaces, edge_sense)
        def _interp_apexvolfun(coord):
            return _ef(coord)*_if(coord)
        _interp_apexvolfun.__name__ = _ef.__name__ + _if.__name__
        return _interp_apexvolfun

    @staticmethod
    def combine_apexvolfun_D(edgefn, D_edgefn, interpfn, basedir,
                            edgefaces, edge_sense):
        _if = interpfn()
        _D_if = _if.D
        _ef = edgefn(basedir, edgefaces, edge_sense)
        _D_ef = D_edgefn(basedir, edgefaces, edge_sense)
        def _interp_apexvolfun_D(coord):
            return _if(coord)*_D_ef(coord) + N.cross(_D_if(coord), _ef(coord))
        _interp_apexvolfun_D.__name__ = _D_ef.__name__ + _D_if.__name__
        return _interp_apexvolfun_D

    def _get_edgefuns(self, *names):
        # We assume the first four local edges are base edges and the last four
        # are apex-edges
        no_base, no_apex = self.no_base_edges, self.no_apex_edges

        base_ifns = self.interpFuns.get_baseedgeFuncs()
        apex_ifns = self.interpFuns.get_apexedgeFuncs() 
        # We want the interpolated edge functions arranged as
        #
        # b_edge1*b_intfn1 ... b_edge4*b_intfn1,
        # a_edge1*a_intfn1 ... a_edge4*a_intfn1 ....
        # b_edge1*b_intfnn ... b_edge4*b_intfnn, 
        # a_edge1*a_intfnn ... a_edge4*a_intfnn
        #
        # where b_ == base, a_ == apex
        #
        # and there are 1 ... n apex and base interpolation functions
        
        def gen_gen(edgefun0, offs=0):
            return lambda ifn: [partial(
                self.combine_edgefun, edgefun0, ifn, i+offs)
                for i in range(no_base)]
        def gen_gen_D(edgefun0, D_edgefun0, offs=0):
            return lambda ifn: [partial(
                self.combine_edgefun_D, edgefun0, D_edgefun0, ifn, i+offs)
                for i in range(no_apex)]

        gen_base_edges = gen_gen(edgefun0_base)
        gen_D_base_edges = gen_gen_D(edgefun0_base, edgefun0_base_D)
        gen_apex_edges = gen_gen(edgefun0_apex, no_base)
        gen_D_apex_edges = gen_gen_D(edgefun0_apex, edgefun0_apex_D, no_base)

        base_ifns_r = base_ifns[::-1]
        no_ifns = max(len(base_ifns), len(apex_ifns))
        edgefuns_f = list(chain(*[
            gen_base_edges(base_ifns[i])+gen_apex_edges(apex_ifns[i])
            for i in range(no_ifns)]))
        edgefuns_r = list(chain(*[
            gen_base_edges(base_ifns_r[i])+gen_apex_edges(apex_ifns[i])
            for i in range(no_ifns)]))
        edgefuns = [eff if not (i%8) in self.reversed_base_edges else efr
                    for i, (eff, efr) in enumerate(zip(edgefuns_f, edgefuns_r))]
        
        edgefuns_D_f = list(chain(*[
            gen_D_base_edges(base_ifns[i])+gen_D_apex_edges(apex_ifns[i])
            for i in range(no_ifns)]))
        edgefuns_D_r = list(chain(*[
            gen_D_base_edges(base_ifns_r[i])+gen_D_apex_edges(apex_ifns[i])
            for i in range(no_ifns)]))
        edgefuns_D = [eff if not (i%8) in self.reversed_base_edges else efr
                    for i,(eff, efr) in enumerate(zip(edgefuns_D_f, edgefuns_D_r))]

        return edgefuns, edgefuns_D

    def _get_apexfacefuns(self, *names):
        ifns = self.interpFuns.get_apexfaceFuncs()
        no_f = self.no_apexfaces
        facefuns = [partial(self.combine_apexfacefun, edgefun0_apex, ifn, i)
                    for ifn in ifns for i in range(no_f)]
        facefuns_D = [partial(self.combine_apexfacefun_D, edgefun0_apex,
                              edgefun0_apex_D, ifn, i)
                    for ifn in ifns for i in range(no_f)]
        return facefuns, facefuns_D

    def _get_basefacefuns(self, *names):
        ifns = self.interpFuns.get_basefaceFuncs()
        basedirs = (1,3)                # Element local edgenos
        facefuns = [partial(self.combine_basefacefun, edgefun0_base, ifn, i)
                    for i in basedirs for ifn in ifns]
        facefuns_D = [partial(self.combine_basefacefun_D, edgefun0_base,
                              edgefun0_base_D, ifn, i)
                    for i in basedirs for ifn in ifns]
        assert (len(facefuns) == len(facefuns_D))
        return facefuns, facefuns_D

    def _get_basevolfuns(self, *names):
        ifns = self.interpFuns.get_basevolFuncs()
        basedirs = (1,3)                # Element local edgenos
        volfuns = [partial(self.combine_basevolfun, edgefun0_base, ifn, i)
                    for i in basedirs for ifn in ifns]
        volfuns_D = [partial(self.combine_basevolfun_D, edgefun0_base,
                             edgefun0_base_D, ifn, i)
                    for i in basedirs for ifn in ifns]
        return volfuns, volfuns_D

    def _get_apexvolfuns(self, *names):
        ifns = self.interpFuns.get_apexvolFuncs()
        volfuns = [partial(self.combine_apexvolfun, edgefun0_apex, ifn, 4)
                   for ifn in ifns]
        #for i in basedirs for ifn in ifns]
        volfuns_D = [partial(self.combine_apexvolfun_D, edgefun0_apex,
                             edgefun0_apex_D, ifn, 4)
                     for ifn in ifns]
                     #for i in basedirs for ifn in ifns]
        return volfuns, volfuns_D

def edgefun0_base(local_edgeno, edgefaces, edge_sense):
    i,j = efaces = edgefaces[local_edgeno]
    assert(i != 4)                      # i must be a triangle face
    assert(j == 4)                      # j must be the base face
    sense = edge_sense[local_edgeno]
    i1 = (i+1) % 4 ; i2 = (i+2) % 4; i3 = (i+3) % 4
    cvi1 = cov_coord_vecs[i1] ; cvi3 = cov_coord_vecs[i3]
    def _edgefun0_base(coord):
        ci1, ci2, ci3, cb = coord[[i1, i2, i3, 4]]
        # graglia99 eq (4), the last 4 functions
        return sense*((ci1*ci2*cvi3 - ci2*ci3*cvi1)/(1-cb))
    return _edgefun0_base

def edgefun0_base_D(local_edgeno, edgefaces, edge_sense):
    i,j = efaces = edgefaces[local_edgeno]
    assert(i != 4)                      # i must be a triangle face
    assert(j == 4)                      # j must be the base face
    sense = edge_sense[local_edgeno]
    i2 = (i+2) % 4 ; im1 = (i-1) % 4
    cvim1 = con_coord_vecs[im1,4] ; cvi2 = con_coord_vecs[i2,4]
    cvi2im1 = con_coord_vecs[i2, im1]
    def _edgefun0_base_D(coord):
        ci2, cim1, cb = coord[[i2, im1, 4]]
        # graglia99 eq (3), the second function
        return sense*(1/(1-cb)*(ci2*cvim1 + cim1*cvi2) + cvi2im1)
    return _edgefun0_base_D
        
def edgefun0_apex(local_edgeno, edgefaces, edge_sense):
    i,j = efaces = edgefaces[local_edgeno]
    assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
    assert(j == (i+1) % 4)
    sense = edge_sense[local_edgeno]
    ii = (i+2) % 4 ; jj = (j+2) % 4
    cvii = cov_coord_vecs[ii] ; cvjj = cov_coord_vecs[jj] ; cvb = cov_coord_vecs[4]
    def _edgefun0_apex(coord):
        cii, cjj, cb = coord[[ii, jj, 4]]
        # graglia99 eq (4), the first 4 functions
        return sense*((cii*cjj*cvb - cjj*cb*cvii - cii*cb*cvjj)/(1-cb) 
                      - cii*cjj*cb*cvb/(1-cb)**2)
    return _edgefun0_apex

def edgefun0_apex_D(local_edgeno, edgefaces, edge_sense):
    i,j = efaces = edgefaces[local_edgeno]
    assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
    assert(j == (i+1) % 4)
    sense = edge_sense[local_edgeno]
    i2 = (i+2) % 4 ; im1 = (i-1) % 4
    cvim1 = con_coord_vecs[im1,4] ; cvi2 = con_coord_vecs[i2,4]
    def _edgefun0_apex_D(coord):
        ci2, cim1, cb = coord[[i2, im1, 4]]
        # graglia99 eq (3), the first function
        return sense*2/(1-cb)*(ci2*cvim1 + cim1*cvi2)
    return _edgefun0_apex_D

class InterpFuncs(object):
    def __init__(self, order=1):
        self.order = order
        def one(*names):
            _one = lambda x: 1.
            _one.D = lambda x: N.array([0., 0., 0.], N.float64)
            return _one
        
        self.apex_efs = {0:[one],
                         1:[self.apex_p1_e1, self.apex_p1_e2]}
        self.base_efs = {0:[one],
                         1:[self.base_p1_e1, self.base_p1_e2]}
        self.apex_ffs = {1:[self.apex_p1_f1, self.apex_p1_f2]}
        self.base_ffs = {1:[self.base_p1_f1, self.base_p1_f2]}
        self.base_vfs = {1:[self.base_p1_v1]}
        self.apex_vfs = {1:[self.apex_p1_v1]}
        
    def get_apexedgeFuncs(self):
        return self.apex_efs[self.order]

    def get_baseedgeFuncs(self):
        return self.base_efs[self.order]

    def get_apexfaceFuncs(self):
        return self.apex_ffs[self.order]
    
    def get_basefaceFuncs(self):
        return self.base_ffs[self.order]
    
    def get_basevolFuncs(self):
        return self.base_vfs[self.order]
    
    def get_apexvolFuncs(self):
        return self.apex_vfs[self.order]

    @staticmethod
    def apex_p1_e1(local_edgeno, edgefaces):
        i,j = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(j == (i+1) % 4)
        def _apex_p1_e1(coord):
            cb = coord[4]
            return 3*cb-1
        cvb = cov_coord_vecs[4]
        def _apex_p1_e1_D(coord):
            return 3*cvb
        _apex_p1_e1.D = _apex_p1_e1_D

        return _apex_p1_e1
    
    @staticmethod
    def apex_p1_e2(local_edgeno, edgefaces):
        i,j = efaces = edgefaces[local_edgeno]
        assert(N.all(efaces != 4))          # Apex edges not on face 4 (pyramid base)
        assert(j == (i+1) % 4)
        i1 = (i+1) % 4; i2 = (i+2) % 4
        def _apex_p1_e2(coord):
            ci1,ci2,cb = coord[[i1, i2, 4]]
            return 3*ci2 - 3*ci1 - 1
        cvi1, cvi2 = cov_coord_vecs[[i1, i2]]
        def _apex_p1_e2_D(coord):
            return 3*cvi2 - 3*cvi1
        _apex_p1_e2.D = _apex_p1_e2_D
        return _apex_p1_e2
    
    
    @staticmethod
    def base_p1_e1(local_edgeno, edgefaces):
        i,j = efaces = edgefaces[local_edgeno]
        assert(i != 4)                      # i must be a triangle face
        assert(j == 4)                      # j must be the base face
        i2 = (i+2) % 4 ; i3 = (i+3) % 4
        def _base_p1_e1(coord):
            ci2, ci3, cb = coord[[i2, i3, 4]]
            return (3*ci3 - 1)*(2*ci2 + 2*cb - 1)
        cvi2, cvi3, cvb = cov_coord_vecs[[i2,i3,4]]
        def _base_p1_e1_D(coord):
            ci2, ci3, cb = coord[[i2, i3, 4]]
            return 2*(3*ci3-1)*cvb + 3*(2*cb+2*ci2-1)*cvi3 + 2*(3*ci3-1)*cvi2
        _base_p1_e1.D = _base_p1_e1_D
        return _base_p1_e1

    @staticmethod
    def base_p1_e2(local_edgeno, edgefaces):
        i,j = efaces = edgefaces[local_edgeno]
        assert(i != 4)                      # i must be a triangle face
        assert(j == 4)                      # j must be the base face
        i1 = (i+1) % 4; i2 = (i+2) % 4 
        def _base_p1_e2(coord):
            ci1, ci2, cb = coord[[i1, i2, 4]]
            return (3*ci1 - 1)*(2*ci2 + 2*cb - 1)
        cvi1, cvi2, cvb = cov_coord_vecs[[i1,i2,4]]
        def _base_p1_e2_D(coord):
            ci1, ci2, cb = coord[[i1, i2, 4]]
            return 2*(3*ci1-1)*cvb + 3*(2*cb+2*ci2-1)*cvi1 + 2*(3*ci1-1)*cvi2
        _base_p1_e2.D = _base_p1_e2_D
        return _base_p1_e2

    @staticmethod
    def apex_p1_f1(local_faceno):
        i1 = (local_faceno + 1) % 4 if local_faceno in (0,3) \
             else (local_faceno - 1) % 4
        def _apex_p1_f1(coord):
            return 3*coord[i1]
        cvi1 = cov_coord_vecs[i1]
        def _apex_p1_f1_D(coord):
            return 3*cvi1
        _apex_p1_f1.D = _apex_p1_f1_D
        _apex_p1_f1.local_edge = 1
        return _apex_p1_f1

    @staticmethod
    def apex_p1_f2(local_faceno):
        im1 = (local_faceno - 1) % 4 if local_faceno in (0,3) \
              else (local_faceno + 1) % 4
        def _apex_p1_f2(coord):
            return 3*coord[im1]
        cvim1 = cov_coord_vecs[im1]
        def _apex_p1_f2_D(coord):
            return 3*cvim1
        _apex_p1_f2.local_edge = 2
        _apex_p1_f2.D = _apex_p1_f2_D
        return _apex_p1_f2
        
    @staticmethod
    def apex_p1_v1():
        def _apex_p1_v1(coord):
            c0,c1,c2,c3,cb = coord
            return c0*c1*c2*c3/(1-cb)**2

        cv0, cv1, cvb = cov_coord_vecs[[0,1,4]]
        def _apex_p1_v1_D(coord):
            lam1, lam2, lam5 = coord[[0,1,4]]
            return (((lam2*(lam5+2*lam1-1)*(lam5+lam2-1))*cv0
                     + (lam1*(lam5+lam1-1)*(lam5+2*lam2-1))*cv1)/(lam5-1)**2
                    -(lam1*lam2*(lam2*lam5+lam1*lam5+2*lam1*lam2-lam2-lam1)
                      )*cvb/(lam5-1)**3)
        _apex_p1_v1.D = _apex_p1_v1_D
        return _apex_p1_v1
        
    @staticmethod
    def base_p1_f1(basedir):
        #basedir is the number of the other face connected to the base face,
        #thereby defining the edge function that needs to be multiplied by this
        #interpolatation func
        i = basedir
        im1 = (i-1) % 4
        def _base_p1_f1(coord):
            ci, cim1 = coord[[i, im1]]
            return 2*ci*(3*cim1-1)
        cvim1, cvi = cov_coord_vecs[[im1, i]]
        def _base_p1_f1_D(coord):
            ci, cim1 = coord[[i, im1]]
            return 2*(3*cim1-1)*cvi+6*ci*cvim1
        _base_p1_f1.D = _base_p1_f1_D
        return _base_p1_f1

    @staticmethod
    def base_p1_f2(basedir):
        i = basedir
        i1 = (i+1) % 4
        def _base_p1_f2(coord):
            ci, ci1 = coord[[i, i1]]
            return 2*ci*(3*ci1-1)
        cvi, cvi1 = cov_coord_vecs[[i, i1]]
        def _base_p1_f2_D(coord):
            ci, ci1 = coord[[i, i1]]
            return 2*(3*ci1-1)*cvi+6*ci*cvi1
        _base_p1_f2.D = _base_p1_f2_D
        return _base_p1_f2

    @staticmethod
    def base_p1_v1(basedir):
        i = basedir
        def _base_p1_v1(coord):
            ci, cb = coord[[i, 4]]
            return 9*ci*cb
        cvi, cvb = cov_coord_vecs[[i, 4]]
        def _base_p1_v1_D(coord):
            ci, cb = coord[[i, 4]]
            return 9*(cb*cvi + ci*cvb)
        _base_p1_v1.D = _base_p1_v1_D
        return _base_p1_v1
    
