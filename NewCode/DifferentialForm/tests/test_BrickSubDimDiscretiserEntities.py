from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.BrickMeshes import FourBricks
from NewCode.tests import xfail
from NewCode import BrickSubDimMesh
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import BrickSubDimDiscretiserEntities, BrickDiscretiserEntities

class test_SubDimDiscretiserEntities(TestCase):
    def setUp(self):
        self.superMesh = BrickMesh.Mesh(FourBricks.listmesh)
        fs = lambda face: face.index == 15
        self.mesh = BrickSubDimMesh.SubSurface(self.superMesh, faceSelector=fs)

    def test_Edge(self):
        edges = BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            BrickSubDimDiscretiserEntities.Edge(mesh=self.mesh, freefun=lambda x: True,
                                           attrs=self.mesh.edges.list_repr()))
        self.assert_(len(edges) == 4)
        assert_equal(edges[:].freeNo, range(4))
        assert_equal(edges[:].superNo, [1,4,10,16])
