from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.tests.PyramMeshes import SixPyram
from NewCode.Meshes import PyramMesh
from NewCode.DifferentialForm import PyramDiscretiserEntities

class test_DiscretiserEntities(TestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        self.mesh = PyramMesh.Mesh(self.listmesh)

    def test_Edge(self):
        edges = PyramDiscretiserEntities.DiscretiserEntityList(
            PyramDiscretiserEntities.Edge(mesh=self.mesh, freefun=lambda x: True,
                                             attrs=self.mesh.edges.list_repr(partial=True)
                                             ))
        noEdges = self.mesh.noEdges
        list_repr = edges.list_repr(partial=True)
        desired_list_repr = {'nodes': self.testMesh.edgeNodes,
                             'nodeCoords': self.testMesh.nodes,
                             'freeNo' : N.arange(noEdges, dtype=N.int32),
                             }
        assert_equal(desired_list_repr, list_repr)

