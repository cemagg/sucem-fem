from __future__ import division

import numpy as N
import scipy

from NewCode.Utilities import gen_lagrange_polys, partial, Struct
from NewCode.Integration import gauss_lobatto_coeffs

cov_coord_vecs = N.array([[1,0,0], [0,1,0], [0,0,1]], N.float64)
v = [cv for cv in cov_coord_vecs.copy()]
con_coord_vecs = {(0,1):  v[2], (1,2):  v[0], (2,0):  v[1],
                  (1,0): -v[2], (2,1): -v[0], (0,2): -v[1],
                  (0,4): -v[2], (1,5): -v[0], (2,3): -v[1],
                  (4,0):  v[2], (5,1):  v[0], (3,2):  v[1],
                  (1,3):  v[2], (2,4):  v[0], (0,5):  v[1],
                  (3,1): -v[2], (4,2): -v[0], (5,0): -v[1]}

face_con_coord_vecs = v

vol_sign = {(0,1,2):  1., (1,2,0):  1., (2,0,1):  1.,
            (3,1,2): -1., (4,2,0): -1., (5,0,1): -1.}

def extended_chebyshev_interp_pts_fn(p): 
    return N.array([-N.cos((2*i+1)*N.pi/(2*p+2))/(2*N.cos(N.pi/(2*p+2))) + 1/2
                    for i in range(p+1)], N.float64)

def gaussian_interp_pts_fn(p):
    return N.real(scipy.special.orthogonal.ps_roots(p+1)[0])

def lobatto_interp_pts_fn(p):
    if p > 0: return gauss_lobatto_coeffs(p)[0]/2 + 0.5
    else: return N.array([0.5])

extended_chebyshev_interp_pts = Struct(
    full=extended_chebyshev_interp_pts_fn,
    reduced=extended_chebyshev_interp_pts_fn)

def calc_interp_pt(local_coordnos, fns):
    ipt = N.zeros(6, N.float64)
    for lc, fn in zip(local_coordnos, fns):
        ipt[lc] = fn.pt 
        ipt[(lc+3)%6] = 1-fn.pt
    return ipt
    

class InterpBasisFuncs(object):
    def __init__(self, int_pts_fns, order, mixed=True, blah=None):
        # Full interpolation order (e.g. normal components for one-form)
        self.fOrder = forder = order    
        # Reduced interpolation order (e.g. Tangential components for one-form)
        self.rOrder = rorder = forder-1 if mixed else forder 
        self.order = (rorder+forder)/2.
        self.fInterpPts = int_pts_fns.full(forder)
        self.rInterpPts = int_pts_fns.reduced(rorder)
        self.fInterpFns = gen_lagrange_polys(self.fInterpPts)
        self.rInterpFns = gen_lagrange_polys(self.rInterpPts)

class BasisSet(Struct):
    btype = 'btype'
    no_edges = 12
    no_faces = 6
    form_dof_order = Struct({'node': None, 'edge': None, 'face': None, 'vol': None})
    form = None
    InterpBasisFuncs = None
    
    @staticmethod
    def _get_interpBasisFuncs(order, mixed):
        raise NotImplementedError

    def __init__(self, order, mixed=True):
        Struct.__init__(self)
        assert(order > 0)
        # Arb order should work, but an error should be raised if an unreasonably
        # high order is requested.
        if order >= 10: raise NotImplementedError 
        ib = self._get_interpBasisFuncs(order, mixed)
        self.fns = Struct()
        self.fns_D = Struct()
        fdo = self.form_dof_order
        if fdo.edge is not None and order >= fdo.edge:
            self.fns.edge, self.fns_D.edge = self._get_edgefuns(ib)
        if fdo.face is not None and order >= fdo.face:
            self.fns.face, self.fns_D.face = self._get_facefuns(ib)
        if fdo.vol is not None and order >= fdo.vol:
            self.fns['vol'] = ib.get_volFuncs()
            self.fns_D['vol'] = ib.get_volFuncs_D()

        if mixed: order -= 0.5
        self.info = Struct(form=self.form, type=self.btype, order=order)

    def _get_edgefuns(self, ib):
        edgefuns = [] ; edgefuns_D = [] ; no_edges = self.no_edges
        for ef, ef_D in zip(ib.get_edgeFuncs(), ib.get_edgeFuncs_D()):
            edgefuns.extend(partial(ef, i) for i in range(no_edges))
            edgefuns_D.extend(partial(ef_D, i) for i in range(no_edges))
        return edgefuns, edgefuns_D

    def _get_facefuns(self, ib):
        facefuns = [] ; facefuns_D = [] ; no_faces = self.no_faces
        for ff, ff_D in zip(ib.get_faceFuncs(), ib.get_faceFuncs_D()):
            for i in range(no_faces):
                facefuns.append(partial(ff, i)) ; facefuns_D.append(partial(ff_D, i))        
        return facefuns, facefuns_D

class BasisSet_rieben04(BasisSet):
    btype = 'rieben04'

    def _get_interpBasisFuncs(self, order, mixed):
        return self.InterpBasisFuncs(extended_chebyshev_interp_pts, order, mixed)


class BasisSet_cohen98(BasisSet):
    btype = 'cohen98'
    
    def _get_interpBasisFuncs(self, order, mixed):
        if mixed: int_pts_fns = Struct(
            full=lobatto_interp_pts_fn, reduced=gaussian_interp_pts_fn)
        else: int_pts_fns = Struct(
            full=lobatto_interp_pts_fn, reduced=lobatto_interp_pts_fn)
        return self.InterpBasisFuncs(int_pts_fns, order, mixed)
