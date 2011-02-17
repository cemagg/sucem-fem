import numpy as N
from numpy.testing import assert_almost_equal

from NewCode.Utilities import Struct
from NewCode import SubDimMesh
from NewCode.DifferentialForm.DiscretiserEntities import geom_entity_names
from NewCode.DifferentialForm import DiscretiserElement

class BasePformSubDimElement(object):
    entityNames = geom_entity_names[:-1] # Remove 3D enitity name
    def __init__(self, attrs, superMesh, mesh, *names, **kwargs):
        self.superMesh = superMesh
        nodeCoords = attrs['nodeCoords']
        super(BasePformSubDimElement, self).__init__(attrs=attrs, mesh=mesh, *names, **kwargs)
        # Because of some hackery, 'nodeCoords' is removed from the attrs hash
        # while initting PformSubDimElement. This should prolly be fixed
        # properly sometime!  FIXME
        attrs['nodeCoords'] = nodeCoords
        self.SubDimMeshModule.FaceElement.__init__(self, attrs)
        # A hack, assming that subdimg faces are always free FIXME
        self._freeNo = N.arange(self.numberof, dtype=N.int32)
        
    def setBasisFunctions(self, superElement):
        self.superParentElement = superElement
        if hasattr(superElement, 'faceRule'):
            self.rule = superElement.faceRule
        bf_ent_names = (g for g in self.entityNames
                        if g in superElement.entityBasisFuns)
        def face_ent_bfs_from_super_edges(super_bfs):
            bfs = []
            fpe = self.getFuncsPerEntity(super_bfs, superElement.COUNT.edge)
            bf_per_edge = len(fpe[0])   # Assuming each edge is identical
            for face_edges in superElement.LOCAL_FACEEDGES:
                bfs.append([])
                bfs[-1].extend(fpe[e][i] for i in range(bf_per_edge)
                               for e in face_edges)
            return bfs

        def face_ent_bfs_from_super_faces(super_bfs):
            return 
        face_ent_bfs_from_super = {
            'edge': face_ent_bfs_from_super_edges,
            'face': lambda sbf: self.getFuncsPerEntity(sbf, superElement.COUNT.face)}
        self.noDOFs = Struct()
        self.noDOFs.element = 0
        self.superEntityBasisFuns = Struct()
        self._basisSet=[list() for i in range(superElement.COUNT.face)]
        for ent_name in (en for en in self.entityNames if en in bf_ent_names):
            sEBF = face_ent_bfs_from_super[ent_name](
                superElement.entityBasisFuns[ent_name])
            for i, bfs in enumerate(sEBF):
                self._basisSet[i].extend(bfs)
            self.superEntityBasisFuns[ent_name] = sEBF
            # We're assuming that every super face is identical
            self.noDOFs[ent_name] = Struct(
                perElement=len(sEBF[0]),
                perEntity=superElement.noDOFs[ent_name].perEntity)
            self.noDOFs.element += len(sEBF[0])

    @staticmethod
    def getFuncsPerEntity(superEntityFuncs, EntitiesPerSuperEl):
        return [superEntityFuncs[i::EntitiesPerSuperEl]
                for i in range(EntitiesPerSuperEl)]
            
    @property
    def basisSet(self):
        return self._basisSet[self.superLocalFaceno]
    
    @property
    def superElement(self):
        # For basis function evaluation any one of the elements connected to
        # the face will do due to tangential continuity, so use the first
        return self.connect2SuperElems.T[0]

    @property
    def superLocalFaceno(self):
        # Element-local faceno to go with the element selected in superElement
        return self.superLocalFacenos.T[0]

class PformSubDimElement(BasePformSubDimElement,
                         DiscretiserElement.NdimPformElement,
                         SubDimMesh.FaceElement):
    SubDimMeshModule = SubDimMesh
    facenos = property(DiscretiserElement.NdimPformElement.entityNo)
    

class BaseOneformSubDimElement(object):
    """
    Oneforms have tangential continuity between elements, and the tangential
    component of the 3D basis function on a given face gives the definition of
    the 2D Sub-basis. The tangential component of a function F on a given face
    is -n x (n x F) where n is a unit normal vector. It does not matter wether n
    is inward or outward normal (this follows from the properties of the cross product)
    """
    entityNames = PformSubDimElement.entityNames[1:] # No nodal 1-form elements

    def refValsAtPoints(self, points, superLocalFaceno=None):
        """
        The actual values depends on which local face of the super element this
        face belongs to. This can be determined implicitly for the current
        face, or it can be specified by superLocalFaceno.
        """
        if superLocalFaceno is None: superLocalFaceno = self.superLocalFaceno
        spel = self.superParentElement
        superPoints = N.array([spel.face_coords2vol_coords(superLocalFaceno, p)
                               for p in points], N.float64)
        return N.array(
            [[bf(pt) for pt in superPoints]
             for bf in self._basisSet[int(superLocalFaceno)]],
            N.float64)

    def bfPhysValAtPoint(self, bf, point):
        coordVecs = self.coordVecs()
        spel = self.superParentElement
        super_point = spel.face_coords2vol_coords(self.superLocalFaceno, point)
        return N.dot(coordVecs, bf(super_point))

    def coordVecs(self):
        se = self.superMesh.elements[self.superElement]
        normal = se.normals()[self.superLocalFaceno]
        # Tangential component is calculated as -n x (n x vec))
        return N.array([-N.cross(normal, N.cross(normal, cv))
                        for cv in se.covBaseVecs().T], N.float64).T

    def normal(self):
        se = self.superMesh.elements[self.superElement]
        return se.normals()[self.superLocalFaceno]
    
    def refVals(self):
        try: return self._refVals[int(self.superLocalFaceno)]
        except AttributeError: pass

        intg_points = self.rule.evalPoints()
        self._refVals = [self.refValsAtPoints(intg_points, superLocalFaceno=i)
                         for i in range(self.superParentElement.COUNT.face)]
        return self._refVals[int(self.superLocalFaceno)]

    def local2global(self, lam):
        #TESTME !!
        spel = self.superParentElement
        spel.index = self.superElement
        return spel.local2global(spel.face_coords2vol_coords(self.superLocalFaceno, lam))

    def global2local(self,r):
        #TESTME !!
        spel = self.superParentElement
        spel.index = self.superElement
        slfno = self.superLocalFaceno
        local = spel.vol_coords2face_coords(slfno, spel.global2local(r))
        assert_almost_equal(
            spel.local2global(spel.face_coords2vol_coords(slfno, local)),
            r, decimal=10)
        return local

class BaseOneformOutwardSubDimElement(BaseOneformSubDimElement):
    @property
    def superElement(self):
        return self._superElement[self.index]

    def outward_normal(self):
        se = self.superMesh.elements[self.superElement]
        return se.outwardNormals()[self.superLocalFaceno]

    @property
    def superLocalFaceno(self):
        # Element-local faceno to go with the element selected in superElement
        return self._superLocalFaceno[self.index]

    def coordVecs(self):
        se = self.superMesh.elements[self.superElement]
        normal = se.outwardNormals()[self.superLocalFaceno]
        # Outward normal cross value:  n x vec)
        return N.array([N.cross(normal, cv)
                        for cv in se.covBaseVecs().T], N.float64).T


class OneformSubDimElement(BaseOneformSubDimElement, PformSubDimElement):
    pass

class OutwardOneformSubDimElement(BaseOneformOutwardSubDimElement,
                                  PformSubDimElement):
    pass

class OneformSubDimElement_DPhysvals(object):
    def coordVecs(self):
        se = self.superMesh.elements[self.superElement]
        normal = se.normals()[self.superLocalFaceno]
        # Normal component is calculated as n.vec*n
        return N.array([N.dot(normal, cv)*normal
                        for cv in se.conBaseVecs().T], N.float64).T

class OneformSubDimElement_D(OneformSubDimElement_DPhysvals,
                             OneformSubDimElement): pass
