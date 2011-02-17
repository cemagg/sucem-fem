import numpy as N

from NewCode.Utilities import Struct
from NewCode import Mesh
from NewCode.DifferentialForm.DiscretiserEntities import geom_entity_names, \
     DiscretiserEntity
from NewCode.DifferentialForm.BoundaryConditions import allfree

class NdimPformElement(DiscretiserEntity):
    dead_elements = set()
    
    def __init__(self, *names, **kwargs):
        """
        Init as for Mesh.Element, with the addition of the freefun kwarg.

        When the DiscretiserElement is to be used as a volume DiscretiserEntity
        the freefun kwarg is used. When there aren't any volume basis
        functions, specify freefun=None
        """
        self.mesh = kwargs['mesh']
        if kwargs['freefun'] is None: kwargs['freefun'] = allfree
        super(NdimPformElement, self).__init__(*names, **kwargs)

    def entityNo(self):
        if isinstance(self.index, slice):
            raise NotImplementedError('Slice handling not supported for volume entities')
        return self.index

    def setIntegrationRule(self, intg_rule):
        self.rule = intg_rule
        try: del self._refVals               # Kill any cached values
        except AttributeError: pass

    def refValsAtPoints(self, points):
        return N.array(
            [[bf(pt) for pt in points]
             for bf in self.basisSet],
            N.float64)

    def refVals(self):
        """
        Evaluate all basis functions on the reference element at all integration points

        Return Value
        ============

        Assuming n basisfunctions, bf1, ... bfn, a list is returned with one
        element per bf, containing the given bf evaluated at each integration
        point pt1, .., ptm eg:

        [[bf1(pt1), ... bf1(ptm)], ... [bfn(pt1), ... bfn(ptm)]]
        
        """
        
        try: return self._refVals
        except AttributeError: pass

        intg_points = self.rule.evalPoints()
        self._refVals = self.refValsAtPoints(intg_points)

        return self._refVals

    def physValsFromRefVals(self, refVals):
        coordVecs = self.coordVecs()
        return N.array([[N.dot(coordVecs, val_at_pt)
                       for val_at_pt in bf_vals]
                      for bf_vals in refVals],
                     N.float64)

    def physVals(self):
        """
        Evaluate all basis functions at all integration points in xyz coordinates

        Assumes the coordinate transform is constant throughout the element
        """
        refVals = self.refVals()
        return self.physValsFromRefVals(refVals)

    def physValsAtPoints(self, points):
        refVals = self.refValsAtPoints(points)
        return self.physValsFromRefVals(refVals)

    def bfPhysValAtPoint(self, bf, point):
        coordVecs = self.coordVecs()
        return N.dot(coordVecs, bf(point))
    
    def physEvalPoints(self):
        return N.array(map(self.local2global, self.rule.evalPoints()))
    
class BasePformElement(NdimPformElement):
    entityNames=geom_entity_names
    volnos = property(NdimPformElement.entityNo) 
    
    def setFaceIntegrationRule(self, face_rule):
        self.faceRule = face_rule
        
    def _get_local_numbering(self):
        return dict(edge=(self.LOCAL_EDGENODES,), face=(self.LOCAL_FACENODES,))

    def setBasisFunctions(self, basisSet):
        """
        Set the basis functions for the reference element. 

        Function Definitions
        ===================

        With m edge functions, n face functions, o volume functions:

        exfy
          local edge no x, edge function no y
        fxfy
          local face no x, face function no y
        vfy
          volume function no y

        Basis Function Signature
        ------------------------

        A given basisfunction generator functio, bf_g,
        should have the signature::

          bf_g(local_edgenodes) for edge functions,
          bf_g(local_facenodes) for face functions,
          bf(coord) for volume functions

          bf_g() should return bf(coord), the actual basis function. Volume
          functions don't need this treatment.

        local_edgenodes/facenodes describe the edge/face numbering convention
        used.  coord is an array of the four simplex coordinates [l1, l2, l3,
        l4].

        For 1-forms bf(coord) should return the value of the basis function
        evaluated at coord i.t.o. covariant simplex basis vectors coefiecients,
        i.e. a return value of [a,b,c,d] implies::

          bf = a*grad(l1) + b*grad(l2) + ... + d*grad(l4)

        For 2-forms, WRITEME!

        Arguments
        =========
        
        basisSet Nodal, Edge, face and volume basis generation functions are
        expected in a dict with keys 'node', 'edge', 'face' and 'vol'

        basisSet['edge']
          [e1f1 ... e6f1 e1f2 ... e6f2 ... e6fm]

        basisSet['face']
          [f1f1 ... f4f1 f1f2 ... f4fn]

        basisSet['vol']
          [vf1 ... vfo]

          
        """
        self.noDOFs = Struct()
        self.entityBasisFuns = Struct()
        self.basisSet = []
        local_numbering = self._get_local_numbering()
        for ent_name in (g for g in self.entityNames if g in basisSet):
            # 6 edges per tet,so the number of edge functions should be a
            # multiple of 6 per tet, multiple of facefuns should be 4, etc.
            basis_fns = basisSet[ent_name]
            assert(len(basis_fns)%self.COUNT[ent_name] == 0)    
            self.noDOFs[ent_name] = Struct(
            perElement=len(basis_fns), perEntity=len(basis_fns)/self.COUNT[ent_name])
            if ent_name in local_numbering:
                ent_basis_funs = [bf_g(*local_numbering[ent_name])
                                  for bf_g in basis_fns]
            else:
                ent_basis_funs = basis_fns
            self.entityBasisFuns[ent_name] = ent_basis_funs
            self.basisSet += ent_basis_funs

        self.noDOFs.element = len(self.basisSet)

class PformElement(BasePformElement, Mesh.Element):
    proxy_attrs = BasePformElement.proxy_attrs + Mesh.Element.proxy_attrs

class OneformElement(PformElement):
    entityNames = PformElement.entityNames[1:]
    coordVecs = PformElement.covBaseVecs

class OneformElement_D(OneformElement):
    coordVecs = PformElement.conBaseVecs

class TwoformElement(PformElement):
    entityNames = PformElement.entityNames[2:]
    coordVecs = PformElement.conBaseVecs

class TwoformElement_D(TwoformElement):
    def physValsFromRefVals (self, refVals):
        return N.array(refVals, N.float64)/N.linalg.det(self.J())

    def bfPhysValAtPoint(self, bf, point):
        return bf(point)/N.linalg.det(self.J())
