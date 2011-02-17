from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.BrickMeshes import FourBricks
from NewCode.tests import xfail
from NewCode import BrickSubDimMesh
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import BrickSubDimDiscretiserEntities, BrickDiscretiser
from NewCode.DifferentialForm import BrickSubDimDiscretiser

class test_PformSubDimDiscretiser(NumpyTestCase):
    def setUp(self):
        self.superMesh = BrickMesh.Mesh(FourBricks.listmesh)
        self.superDisc = BrickDiscretiser.setup_PformDiscretiser(
            self.superMesh, 1)
        fs = lambda face: face.index == 15
        self.mesh = BrickSubDimMesh.SubSurface(self.superMesh, faceSelector=fs)
        freefun = lambda x: True
        self.geomEntities = {
            'edge': BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            BrickSubDimDiscretiserEntities.Edge(mesh=self.mesh, freefun=freefun,
                                                attrs=self.mesh.edges.list_repr())),
                             }
    def test_init(self):
       inst = BrickSubDimDiscretiser.PformSubDimDiscretiser(
           1, self.mesh, self.geomEntities, BrickDiscretiser.Permuter, self.superDisc)
       assert_equal([el.permutation() for el in inst.elements],
                    zip([range(4)], [[2,1,3,0]]))

