from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests.BrickMeshes import OneBrick
from NewCode.DifferentialForm import BasisFunction
from NewCode import Integration
from NewCode.DifferentialForm import BrickDiscretiserElement

import BrickPhysvals

class basetest_PformElement(TestCase):
    TestMesh = OneBrick
    PformElementClass = BrickDiscretiserElement.PformElement
    class Intg(object):
        def evalPoints(self):
            return self.evalpts
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.Intg.evalpts = self.testMesh.test_local_coords
        self.mesh=Struct()
        self.inst = self.PformElementClass(
            {'facenos' : self.testMesh.elementFaces,
             'edgenos' : self.testMesh.elementEdges,
             'nodes' : self.testMesh.elementNodes,
             'nodeCoords' : self.testMesh.nodes,
             'connect2elem' : self.testMesh.elementConnect2Elem,
             'connect2face' : self.testMesh.elementConnect2Face,
             'gridStepSize' : self.testMesh.listmesh['GridStepSize'],
             }, mesh=self.mesh, freefun=None)

class test_PformElement(basetest_PformElement):
    def test_init(self):
        assert self.inst.mesh is self.mesh

    def test_rule(self):
        tet_rule = self.Intg
        self.inst.setIntegrationRule(tet_rule)
        self.assert_(tet_rule is self.inst.rule)

    def test_setBasisFunctions(self):
        edgeFuns = [lambda l: lambda x: [0,0,1],
                    lambda l: lambda x: [0,0,2],
                    lambda l: lambda x: [0,0,3],
                    lambda l: lambda x: [0,0,4],
                    lambda l: lambda x: [0,1,0],
                    lambda l: lambda x: [0,2,0],
                    lambda l: lambda x: [0,3,0],
                    lambda l: lambda x: [0,4,0],
                    lambda l: lambda x: [0,0,1],
                    lambda l: lambda x: [0,0,2],
                    lambda l: lambda x: [0,0,3],
                    lambda l: lambda x: [0,0,4]]
        faceFuns = [lambda l: lambda x: [10,0,0],
                    lambda l: lambda x: [0,10,0],
                    lambda l: lambda x: [0,0,10],
                    lambda l: lambda x: [20,0,0],
                    lambda l: lambda x: [0,20,0],
                    lambda l: lambda x: [0,0,20]]
        basisSet = dict(edge=edgeFuns, face=faceFuns)
        self.inst.setBasisFunctions(basisSet=basisSet)
        assert_equal(len(edgeFuns), self.inst.noDOFs.edge.perElement)
        assert_equal(len(faceFuns), self.inst.noDOFs.face.perElement)
        assert_equal(self.inst.noDOFs.edge.perEntity, 1)
        assert_equal(self.inst.noDOFs.face.perEntity, 1)
        assert_equal(self.inst.noDOFs.element, 18)
        assert_equal([f(None) for f in self.inst.basisSet],
                     [f(None)(None) for f in edgeFuns+faceFuns])
        self.assertRaises(AssertionError, self.inst.setBasisFunctions,
                          dict(edge=edgeFuns[0:4])) 
        self.assertRaises(AssertionError, self.inst.setBasisFunctions,
                          dict(face=faceFuns[0:2])) 
        
    def test_refVals(self):
        edgeFuns = [lambda e: lambda x: x,
                    lambda e: lambda x: x*10,
                    lambda e: lambda x: x*100,
                    lambda e: lambda x: x*1000,
                    lambda e: lambda x: x*10000,
                    lambda e: lambda x: x*100000,
                    lambda e: lambda x: x*2,
                    lambda e: lambda x: x*20,
                    lambda e: lambda x: x*200,
                    lambda e: lambda x: x*2000,
                    lambda e: lambda x: x*20000,
                    lambda e: lambda x: x*200000                    
                    ]
        self.inst.setBasisFunctions(basisSet=dict(edge=edgeFuns))
        tet_rule = self.Intg()
        coords = tet_rule.evalPoints()
        self.inst.setIntegrationRule(tet_rule)
        desired = [[fac*pt for pt in coords]
                   for fac in (1, 10, 100, 1000, 10000, 100000,
                               2, 20, 200, 2000, 20000, 200000,)]
        actual = self.inst.refVals()
        assert_equal(actual, desired)
        # again to check the caching
        self.assert_(self.inst.refVals() is actual)


    def test_physEvalPoints(self):
        tetrule = self.Intg()
        desired = self.testMesh.test_xyz_coords[0]
        self.inst.setIntegrationRule(tetrule)
        assert_almost_equal(self.inst.physEvalPoints(), desired,
                            decimal=15)          

class test_OneformElement(basetest_PformElement):
    PformElementClass = BrickDiscretiserElement.OneformElement

    def test_init(self):
        assert_equal(self.inst.entityNames, ('edge', 'face', 'vol'))

    def test_physVals_efs_R1(self):
        oneform_edge_vals0 = BrickPhysvals.oneform_edge_vals0
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.BrickOneform.edgefuns0))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            oneform_edge_vals0,
                            decimal=15)


class test_OneformElement_D(basetest_PformElement):
    PformElementClass = BrickDiscretiserElement.OneformElement_D
    def test_init(self):
        assert_equal(self.inst.entityNames, ('edge', 'face', 'vol'))

    def test_physVals_efs_R1(self):
        oneform_edge_vals0_D = BrickPhysvals.oneform_edge_vals0_D
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.BrickOneform.edgefuns0_D))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            oneform_edge_vals0_D,
                            decimal=16)

    def test_physVals_efs_genR1(self):
        oneform_edge_vals0_D = BrickPhysvals.oneform_edge_vals0_D
        BrickOneform = BasisFunction.BrickOneform
        from NewCode.Utilities import partial
        ib = BrickOneform.InterpBasisFuncs(BrickOneform.extended_chebyshev_interp_pts, 1)
        edgefuns0_D = tuple(partial(ib.get_edgeFuncs_D()[0], i) for i in range(12))
        self.inst.setBasisFunctions(basisSet=dict(edge=edgefuns0_D))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            oneform_edge_vals0_D,
                            decimal=16)

        
class test_TwoformElement(basetest_PformElement):
    PformElementClass = BrickDiscretiserElement.TwoformElement

    def test_init(self):
        assert_equal(self.inst.entityNames, ('face', 'vol'))

    def test_physVals(self):
        twoform_face_vals0 = BrickPhysvals.twoform_face_vals0
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.BrickTwoform.facefuns0))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            twoform_face_vals0,
                            decimal=16)


class test_TwoformElement(basetest_PformElement):
    PformElementClass = BrickDiscretiserElement.TwoformElement_D

    def test_init(self):
        assert_equal(self.inst.entityNames, ('face', 'vol'))

    def test_physVals(self):
        twoform_face_vals0_D = BrickPhysvals.twoform_face_vals0_D
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.BrickTwoform.facefuns0_D))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(
            self.inst.physVals(),
            twoform_face_vals0_D.reshape(twoform_face_vals0_D.shape+(1,)),
            decimal=16)

