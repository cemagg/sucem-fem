from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.tests.PyramMeshes import SixPyram
from NewCode.tests.TetMeshes import SixPyramConverted

from NewCode.Meshes import Conversions

class test_pyramid2tet(TestCase):
    pyram_testmesh = SixPyram
    tet_testmesh = SixPyramConverted

    def test_pyramid_els2tet_els(self):
        assert_equal(Conversions.pyramid_els2tet_els(self.pyram_testmesh.elementNodes),
                     self.tet_testmesh.listmesh['ElementNodes'])
