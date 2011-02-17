from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests import xfail
from NewCode.tests.BrickMeshes import FourBricks, OneBrick
from NewCode.DifferentialForm import BasisFunction, BrickDiscretiserElement
from NewCode import Integration, BrickSubDimMesh, ProxyList
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import BrickSubDimDiscretiserElement

class _basetest_PformSubDimElement(NumpyTestCase):
    testMesh = None
    SubDimElementClass = BrickSubDimDiscretiserElement.PformSubDimElement
    desiredEntityNames = ('node', 'edge', 'face')
    faceSelector = staticmethod(lambda face: face.index == 15)

    class Intg(object):
        no_steps = 3
        steps = N.linspace(0,1,no_steps)
        lams = N.array([[i,j] for i in steps for j in steps], N.float64) 
        def evalPoints(self):
            return  self.lams

    def setUp(self):
        self.superMesh = BrickMesh.Mesh(self.testMesh.listmesh)
        # Select only faces made up of nodes 1,2,3,5, all on one side of larger tet
        fs = self.faceSelector
        self.mesh = BrickSubDimMesh.SubSurface(self.superMesh, faceSelector=fs)
        self.inst = self.SubDimElementClass(
            self.mesh.elements.list_repr(), self.superMesh, self.mesh, freefun=None)
        self.instList = ProxyList.ProxyList(self.inst)
    
    def test_rule(self):
        intg_rule = self.Intg()
        self.inst.setIntegrationRule(intg_rule)
        self.assert_(intg_rule is self.inst.rule)

class basetest_PformSubDimElementFourBrick(_basetest_PformSubDimElement):
    testMesh = FourBricks
    def test_init(self):
        inst = self.inst
        self.assert_(inst.mesh is self.mesh)
        self.assert_(inst.superMesh is self.superMesh)
        assert_equal(inst.entityNames, self.desiredEntityNames)
        assert_equal(self.instList[:].superFacenos, [15])
        assert_equal(self.instList[:].superLocalFaceno, [5])
        assert_equal(self.instList[:].superElement, [0])

class _basetest_OneformSubDimElement(_basetest_PformSubDimElement):
    SubDimElementClass = BrickSubDimDiscretiserElement.OneformSubDimElement
    desiredEntityNames = ('edge', 'face')

    def setUp(self):
        _basetest_PformSubDimElement.setUp(self)
        self.superElement = BrickDiscretiserElement.OneformElement(
            self.superMesh.elements.list_repr(), mesh=self.superMesh,freefun=None)
        self.superElement.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.BrickOneform.edgefuns0))

class test_OneformSubDimElementFourBrick(
    _basetest_OneformSubDimElement, basetest_PformSubDimElementFourBrick):
    pass

class test_OneformSubDimElementBasisFuns(_basetest_OneformSubDimElement):
    testMesh = OneBrick
    faceSelector = staticmethod(lambda face: True)
    def setUp(self):
        _basetest_OneformSubDimElement.setUp(self)
        self.inst.setBasisFunctions(self.superElement)


    def test_setBasisFunctions(self):
        bfs = N.array(self.superElement.basisSet)
        desired = [bfs[bf_nos].tolist() for bf_nos in
                      N.array([[9,6,10,5], 
                               [9,2,11,1], 
                               [5,3,7,1],  
                               [11,8,12,7],
                               [10,4,12,3],
                               [6,4,8,2]], N.int32)-1]
        assert_equal(self.inst.superEntityBasisFuns['edge'], desired)
        assert_equal(self.inst._basisSet, desired)
        assert_equal([el.basisSet for el in self.instList], [
            desired[el.superLocalFaceno] for el in self.instList] )

        assert_equal(self.inst.noDOFs.edge.perElement, 4)
        assert_equal(self.inst.noDOFs.edge.perEntity, 1)

#     def test_oneform_CTLN_refValsAtPoints(self):
#         assert_almost_equal([
#             self.inst.refValsAtPoints(self.Intg.lams, superLocalFaceno=i)
#             for i in range(4)],
#                             onef_efs_R1_face_refvals, decimal=16)
    
#     def test_refVals(self):
#         self.inst.setIntegrationRule(self.Intg())
#         ref_vals = [el.refVals() for el in self.instList[:]]
#         for i, ref_val in enumerate(ref_vals):
#             self.assert_(ref_val is self.instList[i].refVals())
#         assert_almost_equal(ref_vals, onef_efs_R1_face_refvals, decimal=16)

    def test_oneform_CTLN_physValsAtPoints(self):
        assert_almost_equal([el.physValsAtPoints(self.Intg.lams)
                      for el in self.instList[:]],
                     onef_efs_face_physvals,
                            decimal=15)

onef_efs_face_physvals = N.array([[[[0,0,1/3],[0,0,1/3],[0,0,1/3],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,1/4,0],[0,1/2,0],[0,0,0],[0,1/4,0],[0,1/2,0],[0,0,0],[0,1/4,0],[0,1/2,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,1/3],[0,0,1/3],[0,0,1/3]],[[0,1/2,0],[0,1/4,0],[0,0,0],[0,1/2,0],[0,1/4,0],[0,0,0],[0,1/2,0],[0,1/4,0],[0,0,0]]],[[[0,0,1/3],[0,0,1/3],[0,0,1/3],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,1/4,0],[0,1/2,0],[0,0,0],[0,1/4,0],[0,1/2,0],[0,0,0],[0,1/4,0],[0,1/2,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,1/3],[0,0,1/3],[0,0,1/3]],[[0,1/2,0],[0,1/4,0],[0,0,0],[0,1/2,0],[0,1/4,0],[0,0,0],[0,1/2,0],[0,1/4,0],[0,0,0]]],[[[0,0,1/3],[0,0,1/3],[0,0,1/3],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,1/3],[0,0,1/3],[0,0,1/3]],[[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0]]],[[[0,0,1/3],[0,0,1/3],[0,0,1/3],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,1/6],[0,0,1/6],[0,0,1/6],[0,0,1/3],[0,0,1/3],[0,0,1/3]],[[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0]]],[[[0,1/2,0],[0,1/2,0],[0,1/2,0],[0,1/4,0],[0,1/4,0],[0,1/4,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,1/4,0],[0,1/4,0],[0,1/4,0],[0,1/2,0],[0,1/2,0],[0,1/2,0]],[[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0]]],[[[0,1/2,0],[0,1/2,0],[0,1/2,0],[0,1/4,0],[0,1/4,0],[0,1/4,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0],[0,0,0],[1/2,0,0],[1,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,1/4,0],[0,1/4,0],[0,1/4,0],[0,1/2,0],[0,1/2,0],[0,1/2,0]],[[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0],[1,0,0],[1/2,0,0],[0,0,0]]]], N.float64)
