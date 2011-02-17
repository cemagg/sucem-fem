from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.TestMeshes import InscribedTetMesh
from NewCode.tests import xfail
from NewCode import Mesh, SubDimMesh
from NewCode.DifferentialForm import SubDimDiscretiserEntities, DiscretiserEntities

class test_SubDimDiscretiserEntities(NumpyTestCase):
    def setUp(self):
        self.superMesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        fs = lambda face: 10 in face.connect2elem # Faces of the center-tet
        self.mesh = SubDimMesh.SubSurface(self.superMesh, faceSelector=fs)

    def test_Edge(self):
        edges = SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            SubDimDiscretiserEntities.Edge(mesh=self.mesh, freefun=lambda x: True,
                                           attrs=self.mesh.edges.list_repr()))
        self.assert_(len(edges) == 6)
        assert_equal(edges[:].freeNo, range(6))
        assert_equal(edges[:].superNo, [5, 10, 14, 18, 21, 23])
