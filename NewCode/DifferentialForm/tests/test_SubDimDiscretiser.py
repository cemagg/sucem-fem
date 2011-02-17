from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.TestMeshes import InscribedTetMesh, TwoTets
from NewCode.tests import xfail
from NewCode import Mesh, SubDimMesh
from NewCode.DifferentialForm import SubDimDiscretiserEntities, Discretiser
from NewCode.DifferentialForm import SubDimDiscretiser, allfree

class _basetest_PformSubDimDiscretiser(NumpyTestCase):
    testMesh = None                     
    p = None                            # p-form
    fs = None                           # Face selector function for submesh
    subFreefun = staticmethod(allfree)
    superOrder = 1
    def setUp(self):
        self.superMesh = Mesh.Mesh(self.testMesh.listmesh)
        self.superDisc = Discretiser.setup_PformDiscretiser(
            self.superMesh, self.p, self.superOrder)
        
        self.mesh = SubDimMesh.SubSurface(self.superMesh, faceSelector=self.fs)
        freefun = self.subFreefun 
        self.geomEntities = {'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            SubDimDiscretiserEntities.Edge(mesh=self.mesh, freefun=freefun,
                                           attrs=self.mesh.edges.list_repr())),}
        if 'face' in self.superDisc.basisSet.fns:
            self.geomEntities['face'] = SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
                SubDimDiscretiserEntities.Face(mesh=self.mesh, freefun=freefun,
                                               attrs=self.mesh.elements.list_repr()))

    def _initPformSubDimDiscretiser(self):
        return SubDimDiscretiser.PformSubDimDiscretiser(
           self.p, self.mesh, self.geomEntities,
           Discretiser.Permuter, self.superDisc)


class test_PformSubDimDiscretiserInscribedInit(_basetest_PformSubDimDiscretiser):
    testMesh = InscribedTetMesh
    fs = staticmethod(lambda face: 10 in face.connect2elem) # Faces of the center-tet
    p = 1
    def test_init(self):
       inst = self._initPformSubDimDiscretiser()
       assert_equal([el.permutation() for el in inst.elements],
                    zip([range(3)]*4, [[0, 2, 5,], [3, 0, 1], [1, 2, 4], [3, 5, 4]]))

class _basetest_PformSubDimDiscretiserCoords(_basetest_PformSubDimDiscretiser):
    fs = staticmethod(allfree)
    p = 1

    def test_local2global(self):
        inst = self._initPformSubDimDiscretiser()
        desired = N.array([N.average(f.nodeCoords, axis=0)
                           for f in self.superMesh.faces], N.float64)
        assert_almost_equal(
            [el.local2global(N.array([1/3.]*3, N.float64))
             for el in inst.elements],
            desired)

class test_PformSubDimDiscretiserCoordsInscribed(_basetest_PformSubDimDiscretiserCoords):
    testMesh = InscribedTetMesh

class test_PformSubDimDiscretiserPermutation(_basetest_PformSubDimDiscretiserCoords):
    testMesh = TwoTets
    superOrder = 2
    fs = staticmethod(allfree)
    
    def test_Permutation(self):
        inst = self._initPformSubDimDiscretiser()
        # This result is from an assumed-good run, was not hand verified
        desired = [(N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 0,  6,  2,  1,  7,  3, 18, 19], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 0,  8,  4,  1,  9,  5, 20, 21], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 2, 10,  4,  3, 11,  5, 22, 23], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 6, 10,  8,  7, 11,  9, 24, 25], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 2, 14, 12,  3, 15, 13, 26, 27], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([ 4, 16, 12,  5, 17, 13, 28, 29], N.int32)),
                   (N.array([0, 1, 2, 3, 4, 5, 6, 7], N.int32),
                    N.array([10, 16, 14, 11, 17, 15, 30, 31], N.int32))]
        assert_equal([el.permutation() for el in inst.elements], desired)
