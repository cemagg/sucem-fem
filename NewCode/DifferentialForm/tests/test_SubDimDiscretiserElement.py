from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests import xfail
from NewCode.tests.TestMeshes import FlatTet, InscribedTetMesh
from NewCode.DifferentialForm import BasisFunction, DiscretiserElement
from NewCode import Integration, Mesh, SubDimMesh, ProxyList
from NewCode.DifferentialForm import SubDimDiscretiserElement

class _basetest_PformSubDimElement(TestCase):
    testMesh = None
    SubDimElementClass = SubDimDiscretiserElement.PformSubDimElement
    desiredEntityNames = ('node', 'edge', 'face')
    faceSelector = staticmethod(lambda face: face.onBoundary and \
                   len(set([1,2,3,5]) - set(face.nodes+1)) == 1)

    class Intg(object):
        lams = [(i,j,1-i-j)
                for i in N.arange(0,1.01, 1/3.)
                for j in N.arange(0,1.01-i, 1/3.)
                ]
        lams.reverse()          # To match the maxima order
        lams = N.array(lams, N.float64) # Basis fn requires ndarray
        def evalPoints(self):
            return  self.lams

    def setUp(self):
        self.superMesh = Mesh.Mesh(self.testMesh.listmesh)
        # Select only faces made up of nodes 1,2,3,5, all on one side of larger tet
        fs = self.faceSelector
        self.mesh = SubDimMesh.SubSurface(self.superMesh, faceSelector=fs)
        self.inst = self.SubDimElementClass(
            self.mesh.elements.list_repr(), self.superMesh, self.mesh, freefun=None)
        self.instList = ProxyList.ProxyList(self.inst)

    def test_rule(self):
        intg_rule = self.Intg()
        self.inst.setIntegrationRule(intg_rule)
        self.assert_(intg_rule is self.inst.rule)

class basetest_PformSubDimElementInscribed(_basetest_PformSubDimElement):
    testMesh = InscribedTetMesh
    def test_init(self):
        inst = self.inst
        self.assert_(inst.mesh is self.mesh)
        self.assert_(inst.superMesh is self.superMesh)
        assert_equal(inst.entityNames, self.desiredEntityNames)
        assert_equal(self.instList[:].superFacenos+1, [5,13,17])
        assert_equal(self.instList[:].superLocalFaceno, [0,0,0])
        assert_equal(self.instList[:].superElement+1, [2,4, 5])

class test_getFuncsPerEntity(TestCase):
    def test_getFuncsPerEntity(self):
        assert_equal(SubDimDiscretiserElement.OneformSubDimElement.getFuncsPerEntity(
            ['e1f1', 'e2f1', 'e1f2', 'e2f2', 'e1f3', 'e2f3'], 2),
                     [['e1f1', 'e1f2', 'e1f3'],
                      ['e2f1', 'e2f2', 'e2f3']])

class _basetest_OneformSubDimElement(_basetest_PformSubDimElement):
    SubDimElementClass = SubDimDiscretiserElement.OneformSubDimElement
    desiredEntityNames = ('edge', 'face')
    basisSet = BasisFunction.Oneform.basis_set(1, mixed=True)

    def setUp(self):
        _basetest_PformSubDimElement.setUp(self)
        self.superElement = DiscretiserElement.OneformElement(
            self.superMesh.elements.list_repr(), mesh=self.superMesh, freefun=None)
        self.superElement.setBasisFunctions(basisSet=self.basisSet.fns)

class test_OneformSubDimElementInscibed(
    _basetest_OneformSubDimElement, basetest_PformSubDimElementInscribed):
    pass

class test_OneformSubDimElementBasisFuns(_basetest_OneformSubDimElement):
    testMesh = FlatTet
    faceSelector = staticmethod(lambda face: True)
    def setUp(self):
        _basetest_OneformSubDimElement.setUp(self)
        self.inst.setBasisFunctions(self.superElement)


    def test_setBasisFunctions(self):
        self.inst.setBasisFunctions(self.superElement)
        bfs = N.array(self.superElement.basisSet)
        desired = map(list,
                      [bfs[[0,3,1]], bfs[[0,4,2]], bfs[[1,5,2]], bfs[[3,5,4]]])
        assert_equal(self.inst.superEntityBasisFuns['edge'], desired)
        assert_equal(self.inst._basisSet, desired)
        assert_equal([el.basisSet for el in self.instList], desired)
        assert_equal(self.inst.noDOFs.edge.perElement, 3)
        assert_equal(self.inst.noDOFs.edge.perEntity, 1)

    def test_oneform_CTLN_refValsAtPoints(self):
        assert_almost_equal([
            self.inst.refValsAtPoints(self.Intg.lams, superLocalFaceno=i)
            for i in range(4)],
                            onef_efs_R1_face_refvals, decimal=16)
    
    def test_refVals(self):
        self.inst.setIntegrationRule(self.Intg())
        ref_vals = [el.refVals() for el in self.instList[:]]
        for i, ref_val in enumerate(ref_vals):
            self.assert_(ref_val is self.instList[i].refVals())
        assert_almost_equal(ref_vals, onef_efs_R1_face_refvals, decimal=16)

    def test_oneform_CTLN_physValsAtPoints(self):
        assert_almost_equal([el.physValsAtPoints(self.Intg.lams)
                      for el in self.instList[:]],
                     onef_efs_R1_face_physvals,
                            decimal=15)

class test_OneformSubDimElementBasisFunsQTCuN(_basetest_OneformSubDimElement):
    testMesh = FlatTet
    basisSet = BasisFunction.Oneform.basis_set(3, mixed=True)
    faceSelector = staticmethod(lambda face: True)
    def test_setBasisFunctions(self):
        self.inst.setBasisFunctions(self.superElement)
        bfs = N.array(self.superElement.basisSet)
        edge_desired = [bfs[bf_nos].tolist() for bf_nos in
                        N.array([[0, 3, 1, 6, 9, 7, 12, 15, 13],
                                 [0, 4, 2, 6, 10, 8, 12, 16, 14],
                                 [1, 5, 2, 7, 11, 8, 13, 17, 14],
                                 [3, 5, 4, 9, 11, 10, 15, 17, 16]], N.int32)]
        assert_equal(self.inst.superEntityBasisFuns['edge'], edge_desired)
        no_edge_funs = 3*6
        face_desired = [bfs[no_edge_funs+bf_nos].tolist() for bf_nos in
                        N.array([[0, 4, 8, 12, 16, 20],
                                 [1, 5, 9, 13, 17, 21],
                                 [2, 6, 10, 14, 18, 22],
                                 [3, 7, 11, 15, 19, 23]], N.int32)]
        assert_equal(self.inst.superEntityBasisFuns['face'], face_desired)
        desired = [ed+fd for ed, fd in zip(edge_desired,face_desired)]
        assert_equal(self.inst._basisSet, desired)
        assert_equal([el.basisSet for el in self.instList], desired)
        per_edge, per_face = 3,6
        assert_equal(self.inst.noDOFs.edge.perElement, 3*per_edge)
        assert_equal(self.inst.noDOFs.edge.perEntity, per_edge)
    
    def test_bfPhysValAtPoint(self):
        pts = self.Intg().evalPoints()
        self.inst.setBasisFunctions(self.superElement)
        for el in self.instList:
            assert_equal([[el.bfPhysValAtPoint(bf, pt) for pt in pts]
                          for bf in el.basisSet],
                         el.physValsAtPoints(pts))

onef_efs_R1_face_refvals = N.array([[[[0,1,0,0],[-1/3,2/3,0,0],[0,2/3,0,0],[-2/3,1/3,0,0],[-1/3,1/3,0,0],[0,1/3,0,0],[-1,0,0,0],[-2/3,0,0,0],[-1/3,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,1/3,0],[0,-1/3,0,0],[0,0,2/3,0],[0,-1/3,1/3,0],[0,-2/3,0,0],[0,0,1,0],[0,-1/3,2/3,0],[0,-2/3,1/3,0],[0,-1,0,0]],[[0,0,1,0],[0,0,2/3,0],[-1/3,0,2/3,0],[0,0,1/3,0],[-1/3,0,1/3,0],[-2/3,0,1/3,0],[0,0,0,0],[-1/3,0,0,0],[-2/3,0,0,0],[-1,0,0,0]]],[[[0,1,0,0],[-1/3,2/3,0,0],[0,2/3,0,0],[-2/3,1/3,0,0],[-1/3,1/3,0,0],[0,1/3,0,0],[-1,0,0,0],[-2/3,0,0,0],[-1/3,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,1/3],[0,-1/3,0,0],[0,0,0,2/3],[0,-1/3,0,1/3],[0,-2/3,0,0],[0,0,0,1],[0,-1/3,0,2/3],[0,-2/3,0,1/3],[0,-1,0,0]],[[0,0,0,1],[0,0,0,2/3],[-1/3,0,0,2/3],[0,0,0,1/3],[-1/3,0,0,1/3],[-2/3,0,0,1/3],[0,0,0,0],[-1/3,0,0,0],[-2/3,0,0,0],[-1,0,0,0]]],[[[0,0,1,0],[-1/3,0,2/3,0],[0,0,2/3,0],[-2/3,0,1/3,0],[-1/3,0,1/3,0],[0,0,1/3,0],[-1,0,0,0],[-2/3,0,0,0],[-1/3,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,1/3],[0,0,-1/3,0],[0,0,0,2/3],[0,0,-1/3,1/3],[0,0,-2/3,0],[0,0,0,1],[0,0,-1/3,2/3],[0,0,-2/3,1/3],[0,0,-1,0]],[[0,0,0,1],[0,0,0,2/3],[-1/3,0,0,2/3],[0,0,0,1/3],[-1/3,0,0,1/3],[-2/3,0,0,1/3],[0,0,0,0],[-1/3,0,0,0],[-2/3,0,0,0],[-1,0,0,0]]],[[[0,0,1,0],[0,-1/3,2/3,0],[0,0,2/3,0],[0,-2/3,1/3,0],[0,-1/3,1/3,0],[0,0,1/3,0],[0,-1,0,0],[0,-2/3,0,0],[0,-1/3,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,1/3],[0,0,-1/3,0],[0,0,0,2/3],[0,0,-1/3,1/3],[0,0,-2/3,0],[0,0,0,1],[0,0,-1/3,2/3],[0,0,-2/3,1/3],[0,0,-1,0]],[[0,0,0,1],[0,0,0,2/3],[0,-1/3,0,2/3],[0,0,0,1/3],[0,-1/3,0,1/3],[0,-2/3,0,1/3],[0,0,0,0],[0,-1/3,0,0],[0,-2/3,0,0],[0,-1,0,0]]]], N.float64)

onef_efs_R1_face_physvals = N.array([[[[-1,-1,0],[-1/3,-2/3,1/3],[-2/3,-2/3,0],[1/3,-1/3,2/3],[0,-1/3,1/3],[-1/3,-1/3,0],[1,0,1],[2/3,0,2/3],[1/3,0,1/3],[0,0,0]],[[0,0,0],[2/3,1/3,1/3],[1/3,1/3,0],[4/3,2/3,2/3],[1,2/3,1/3],[2/3,2/3,0],[2,1,1],[5/3,1,2/3],[4/3,1,1/3],[1,1,0]],[[2,1,1],[4/3,2/3,2/3],[5/3,2/3,1],[2/3,1/3,1/3],[1,1/3,2/3],[4/3,1/3,1],[0,0,0],[1/3,0,1/3],[2/3,0,2/3],[1,0,1]]],[[[-1,0,1],[-1/3,-1/3,2/3],[-2/3,0,2/3],[1/3,-2/3,1/3],[0,-1/3,1/3],[-1/3,0,1/3],[1,-1,0],[2/3,-2/3,0],[1/3,-1/3,0],[0,0,0]],[[0,0,0],[2/3,-1/3,-1/3],[1/3,0,-1/3],[4/3,-2/3,-2/3],[1,-1/3,-2/3],[2/3,0,-2/3],[2,-1,-1],[5/3,-2/3,-1],[4/3,-1/3,-1],[1,0,-1]],[[2,-1,-1],[4/3,-2/3,-2/3],[5/3,-1,-2/3],[2/3,-1/3,-1/3],[1,-2/3,-1/3],[4/3,-1,-1/3],[0,0,0],[1/3,-1/3,0],[2/3,-2/3,0],[1,-1,0]]],[[[3/11,12/11,21/11],[4/11,5/11,17/11],[2/11,8/11,14/11],[5/11,-2/11,13/11],[3/11,1/11,10/11],[1/11,4/11,7/11],[6/11,-9/11,9/11],[4/11,-6/11,6/11],[2/11,-3/11,3/11],[0,0,0]],[[0,0,0],[1/11,-7/11,-4/11],[-1/11,-4/11,-7/11],[2/11,-14/11,-8/11],[0,-1,-1],[-2/11,-8/11,-14/11],[3/11,-21/11,-12/11],[1/11,-18/11,-15/11],[-1/11,-15/11,-18/11],[-3/11,-12/11,-21/11]],[[3/11,-21/11,-12/11],[2/11,-14/11,-8/11],[4/11,-17/11,-5/11],[1/11,-7/11,-4/11],[3/11,-10/11,-1/11],[5/11,-13/11,2/11],[0,0,0],[2/11,-3/11,3/11],[4/11,-6/11,6/11],[6/11,-9/11,9/11]]],[[[3/11,21/11,12/11],[4/11,17/11,5/11],[2/11,14/11,8/11],[5/11,13/11,-2/11],[3/11,10/11,1/11],[1/11,7/11,4/11],[6/11,9/11,-9/11],[4/11,6/11,-6/11],[2/11,3/11,-3/11],[0,0,0]],[[0,0,0],[1/11,-4/11,-7/11],[-1/11,-7/11,-4/11],[2/11,-8/11,-14/11],[0,-1,-1],[-2/11,-14/11,-8/11],[3/11,-12/11,-21/11],[1/11,-15/11,-18/11],[-1/11,-18/11,-15/11],[-3/11,-21/11,-12/11]],[[3/11,-12/11,-21/11],[2/11,-8/11,-14/11],[4/11,-5/11,-17/11],[1/11,-4/11,-7/11],[3/11,-1/11,-10/11],[5/11,2/11,-13/11],[0,0,0],[2/11,3/11,-3/11],[4/11,6/11,-6/11],[6/11,9/11,-9/11]]]], N.float64)
