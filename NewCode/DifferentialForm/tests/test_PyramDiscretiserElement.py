from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode import ProxyList
from NewCode.Utilities import Struct
from NewCode.tests import xfail
from NewCode.tests.PyramMeshes import SixPyram
from NewCode.DifferentialForm.BasisFunction import PyramOneform
from NewCode.Meshes import PyramMesh
from NewCode.DifferentialForm import PyramDiscretiserElement, PyramDiscretiserEntities
import pyram_physvals

class _basetest_PformElement(TestCase):
    PformElementClass = PyramDiscretiserElement.PformElement
    TestMesh = SixPyram
    class Intg(object):
        def evalPoints(self):
            return self.lams

    def setUp(self):
        self.Intg.lams = self.TestMesh.test_local_coords
        self.testMesh = self.TestMesh()
        self.mesh = PyramMesh.Mesh(self.testMesh.listmesh)
        self.inst = self.PformElementClass(self.mesh.elements.list_repr(partial=True),
                                           mesh=self.mesh, freefun=None)
        self.inst.setIntegrationRule(self.Intg())
        self.inst_list = ProxyList.ProxyList(self.inst)
        
class test_PformElement(_basetest_PformElement):
    def test_init(self):
        assert self.inst.mesh is self.mesh

class test_OneformElement(test_PformElement):
    PformElementClass = PyramDiscretiserElement.OneformElement    

    def test_physVals_efs0(self):
        self.inst.setBasisFunctions(PyramOneform.basis_set(1, mixed=True).fns)
        assert_almost_equal([el.physVals() for el in self.inst_list],
                            pyram_physvals.oneform_edge0, decimal=15)

class test_OneformElement_D(test_PformElement):
    PformElementClass = PyramDiscretiserElement.OneformElement_D

    def test_physVals_D_efs0(self):
        self.inst.setBasisFunctions(PyramOneform.basis_set(1, mixed=True).fns_D)
        assert_almost_equal([el.physVals() for el in self.inst_list],
                            pyram_physvals.D_oneform_edge0, decimal=15)
