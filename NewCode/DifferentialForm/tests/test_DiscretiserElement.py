from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests.TestMeshes import FlatTet
from NewCode.DifferentialForm import BasisFunction
from NewCode import Integration
from NewCode.DifferentialForm import DiscretiserElement, DiscretiserEntities
from NewCode.DifferentialForm.BoundaryConditions import allconstrained

class test_volFreenos(NumpyTestCase):
    PformElementClass = DiscretiserElement.PformElement
    testMesh = FlatTet

    def test_noFree(self):
        mesh=Struct()
        inst = DiscretiserEntities.DiscretiserEntityListNoBoundary(
            self.PformElementClass(
            {'facenos' : self.testMesh.listmesh['ElementFaces'],
             'edgenos' : self.testMesh.listmesh['ElementEdges'],
             'nodes' : self.testMesh.listmesh['ElementNodes'],
             'nodeCoords' : self.testMesh.listmesh['Nodes'],
             'connect2elem' : self.testMesh.listmesh['ElementConnect2Elem'],
             'connect2face' : self.testMesh.listmesh['ElementConnect2Face']
             }, mesh=mesh, freefun=allconstrained))
        assert_equal(inst.noFree, 0)

        inst = DiscretiserEntities.DiscretiserEntityListNoBoundary(
            self.PformElementClass(
            {'facenos' : self.testMesh.listmesh['ElementFaces'],
             'edgenos' : self.testMesh.listmesh['ElementEdges'],
             'nodes' : self.testMesh.listmesh['ElementNodes'],
             'nodeCoords' : self.testMesh.listmesh['Nodes'],
             'connect2elem' : self.testMesh.listmesh['ElementConnect2Elem'],
             'connect2face' : self.testMesh.listmesh['ElementConnect2Face']
             }, mesh=mesh, freefun=None))
        assert_equal(inst.noFree, 1)
        
class _basetest_PformElement(NumpyTestCase, FlatTet):
    PformElementClass = DiscretiserElement.PformElement
    class Intg(object):
        def evalPoints(self):
            lams = [(i,j,k,1-i-j-k)
                    for i in N.arange(0,1.01, 1/3.)
                    for j in N.arange(0,1.01-i, 1/3.)
                    for k in N.arange(0, 1.01-i-j, 1/3.)
                    ]
            lams.reverse()          # To match the maxima order
            return N.array(lams, N.float64) # Basis fn requires ndarray

    def setUp(self):
        self.mesh=Struct()
        self.inst = self.PformElementClass(
            {'facenos' : self.listmesh['ElementFaces'],
             'edgenos' : self.listmesh['ElementEdges'],
             'nodes' : self.listmesh['ElementNodes'],
             'nodeCoords' : self.listmesh['Nodes'],
             'connect2elem' : self.listmesh['ElementConnect2Elem'],
             'connect2face' : self.listmesh['ElementConnect2Face']
             }, mesh=self.mesh, freefun=None)

class test_PformElement(_basetest_PformElement):
      
    def test_init(self):
        assert self.inst.mesh is self.mesh

    def test_rule(self):
        tet_rule = Integration.TetIntegrator(2)
        self.inst.setIntegrationRule(tet_rule)
        self.assert_(tet_rule is self.inst.rule)

    def test_setBasisFunctions(self):
        edgeFuns = [lambda l: lambda x: [1,0,0,0],
                    lambda l: lambda x: [0,1,0,0],
                    lambda l: lambda x: [0,0,1,0],
                    lambda l: lambda x: [0,0,0,1],
                    lambda l: lambda x: [2,0,0,0],
                    lambda l: lambda x: [0,2,0,0]]
        faceFuns = [lambda l: lambda x: [10,0,0,0],
                    lambda l: lambda x: [0,10,0,0],
                    lambda l: lambda x: [0,0,10,0],
                    lambda l: lambda x: [0,0,0,10]]
        basisSet = dict(edge=edgeFuns, face=faceFuns)
        self.inst.setBasisFunctions(basisSet=basisSet)
        assert_equal(len(edgeFuns), self.inst.noDOFs.edge.perElement)
        assert_equal(len(faceFuns), self.inst.noDOFs.face.perElement)
        assert_equal(self.inst.noDOFs.edge.perEntity, 1)
        assert_equal(self.inst.noDOFs.face.perEntity, 1)
        assert_equal(self.inst.noDOFs.element, 10)
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
                    lambda e: lambda x: x*100000
                    ]
        self.inst.setBasisFunctions(basisSet=dict(edge=edgeFuns))
        tet_rule = Integration.TetIntegrator(2)
        coords = tet_rule.evalPoints()
        self.inst.setIntegrationRule(tet_rule)
        desired = [[fac*pt for pt in coords]
                   for fac in (1, 10, 100, 1000, 10000, 100000)
                   ]
        actual = self.inst.refVals()
        assert_equal(actual, desired)
        # again to check the caching
        self.assert_(self.inst.refVals() is actual)


    def test_physEvalPoints(self):
        tetrule = Struct(evalPoints=lambda :self.test_local_coords)
        desired = self.test_xyz_coords
        self.inst.setIntegrationRule(tetrule)
        assert_almost_equal(self.inst.physEvalPoints(), desired,
                            decimal=15)                     

class _basetest_PformElementPhysval(_basetest_PformElement):
    def test_bfPhysValAtPoint(self):
        pts = self.Intg().evalPoints()
        self.inst.setBasisFunctions(basisSet=self.testBFset)
        assert_equal([[self.inst.bfPhysValAtPoint(bf, pt) for pt in pts]
                      for bf in self.inst.basisSet],
                     self.inst.physValsAtPoints(pts))


class test_OneformElement(_basetest_PformElementPhysval):
    PformElementClass = DiscretiserElement.OneformElement
    testBFset = BasisFunction.Oneform.basis_set(1)['fns']

    def test_init(self):
        assert_equal(self.inst.entityNames, ('edge', 'face', 'vol'))

    def test_physVals_efs_R1(self):
        global oneform_edge_vals0
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.Oneform.edgefuns0))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            oneform_edge_vals0,
                            decimal=15)

    def test_physVals_efs_G1(self):
        global oneform_edge_vals1
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.Oneform.edgefuns1))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            oneform_edge_vals1,
                            decimal=15)        

    def test_physVals_ffs_R2(self):
        global onef_ffs_R2_physvals
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.Oneform.facefunsR2))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            onef_ffs_R2_physvals,
                            decimal=15)

class test_OneformElement_D(_basetest_PformElementPhysval):
    PformElementClass = DiscretiserElement.OneformElement_D
    testBFset = BasisFunction.Oneform.basis_set(1)['fns_D']
    def test_init(self):
        assert_equal(self.inst.entityNames, ('edge', 'face', 'vol'))

    class Intg1(object):
        def evalPoints(self):
            return N.array([[.25,.25,.25,.25]], N.float64)

    def test_physVals_efs_R1(self):
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.Oneform.edgefuns0_D))
        self.inst.setIntegrationRule(self.Intg1())
        assert_almost_equal(self.inst.physVals(),
                             [[[0,3,3]],[[3,3,-6]],[[-3,-6,3]],
                              [[-3,6,-3]],[[3,-3,6]],[[0,9,-9]]],
                            decimal=14)

    def test_physVals_efs_G1(self):
        self.inst.setBasisFunctions(basisSet=dict(
            edge=BasisFunction.Oneform.edgefuns1_D))
        self.inst.setIntegrationRule(self.Intg1())
        assert_equal(self.inst.physVals(),
                     N.zeros((6,1,3), N.float64))
        
    def test_physVals_efs_R2(self):
        global onef_ffs_R2_physvals_D
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.Oneform.facefunsR2_D))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            onef_ffs_R2_physvals_D,
                            decimal=14)

        
class test_TwoformElement(_basetest_PformElementPhysval):
    PformElementClass = DiscretiserElement.TwoformElement
    testBFset = BasisFunction.Twoform.basis_set(1)['fns']
    
    def test_init(self):
        assert_equal(self.inst.entityNames, ('face', 'vol'))

    def test_physVals(self):
        global twoform_face_vals0
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.Twoform.facefuns0))

        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            twoform_face_vals0,
                            decimal=14)

class test_TwoformElement_D(test_TwoformElement):
    PformElementClass = DiscretiserElement.TwoformElement_D
    testBFset = BasisFunction.Twoform.basis_set(1)['fns_D']
    class Intg(object):
        def evalPoints(self):
            return N.array([[.25,.25,.25,.25]], N.float64)
    def test_physVals(self):
        self.inst.setBasisFunctions(basisSet=dict(
            face=BasisFunction.Twoform.facefuns0_D))
        self.inst.setIntegrationRule(self.Intg())
        assert_almost_equal(self.inst.physVals(),
                            [[[27]], [[-27]], [[27]], [[-27]]],
                            decimal=15)

# Standard CT/LN basis fns evaluated on FlatTet element.


oneform_edge_vals0 = N.array([[[-3/2,-1/2,1/2],[-1/2,-1/2,1/2],[-1,-1/3,1/3],[-1,-1/3,1/3],[1/2,-1/2,1/2],[0,-1/3,1/3],[0,-1/3,1/3],[-1/2,-1/6,1/6],[-1/2,-1/6,1/6],[-1/2,-1/6,1/6],[3/2,-1/2,1/2],[1,-1/3,1/3],[1,-1/3,1/3],[1/2,-1/6,1/6],[1/2,-1/6,1/6],[1/2,-1/6,1/6],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[3/2,3/2,3/2],[1,1,1],[3/2,5/6,7/6],[1,1,1],[1/2,1/2,1/2],[1,1/3,2/3],[1/2,1/2,1/2],[3/2,1/6,5/6],[1,1/3,2/3],[1/2,1/2,1/2],[0,0,0],[1/2,-1/6,1/6],[0,0,0],[1,-1/3,1/3],[1/2,-1/6,1/6],[0,0,0],[3/2,-1/2,1/2],[1,-1/3,1/3],[1/2,-1/6,1/6],[0,0,0]],[[3/2,-3/2,-3/2],[1,-1,-1],[1,-1,-1],[3/2,-7/6,-5/6],[1/2,-1/2,-1/2],[1/2,-1/2,-1/2],[1,-2/3,-1/3],[1/2,-1/2,-1/2],[1,-2/3,-1/3],[3/2,-5/6,-1/6],[0,0,0],[0,0,0],[1/2,-1/6,1/6],[0,0,0],[1/2,-1/6,1/6],[1,-1/3,1/3],[0,0,0],[1/2,-1/6,1/6],[1,-1/3,1/3],[3/2,-1/2,1/2]],[[0,0,0],[1/2,1/2,1/2],[1/2,1/6,-1/6],[0,0,0],[1,1,1],[1,2/3,1/3],[1/2,1/2,1/2],[1,1/3,-1/3],[1/2,1/6,-1/6],[0,0,0],[3/2,3/2,3/2],[3/2,7/6,5/6],[1,1,1],[3/2,5/6,1/6],[1,2/3,1/3],[1/2,1/2,1/2],[3/2,1/2,-1/2],[1,1/3,-1/3],[1/2,1/6,-1/6],[0,0,0]],[[0,0,0],[1/2,-1/2,-1/2],[0,0,0],[1/2,1/6,-1/6],[1,-1,-1],[1/2,-1/2,-1/2],[1,-1/3,-2/3],[0,0,0],[1/2,1/6,-1/6],[1,1/3,-1/3],[3/2,-3/2,-3/2],[1,-1,-1],[3/2,-5/6,-7/6],[1/2,-1/2,-1/2],[1,-1/3,-2/3],[3/2,-1/6,-5/6],[0,0,0],[1/2,1/6,-1/6],[1,1/3,-1/3],[3/2,1/2,-1/2]],[[0,0,0],[0,0,0],[1/2,-1/2,-1/2],[-1/2,-1/2,-1/2],[0,0,0],[1/2,-1/2,-1/2],[-1/2,-1/2,-1/2],[1,-1,-1],[0,-1,-1],[-1,-1,-1],[0,0,0],[1/2,-1/2,-1/2],[-1/2,-1/2,-1/2],[1,-1,-1],[0,-1,-1],[-1,-1,-1],[3/2,-3/2,-3/2],[1/2,-3/2,-3/2],[-1/2,-3/2,-3/2],[-3/2,-3/2,-3/2]]]                           , N.float64)
oneform_edge_vals1 = N.array([[[-3/2,-1/2,1/2],[-3/2,-1/6,1/6],[-1,-1/3,1/3],[-1,-1/3,1/3],[-3/2,1/6,-1/6],[-1,0,0],[-1,0,0],[-1/2,-1/6,1/6],[-1/2,-1/6,1/6],[-1/2,-1/6,1/6],[-3/2,1/2,-1/2],[-1,1/3,-1/3],[-1,1/3,-1/3],[-1/2,1/6,-1/6],[-1/2,1/6,-1/6],[-1/2,1/6,-1/6],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[3/2,3/2,3/2],[1,1,1],[1/2,7/6,5/6],[1,1,1],[1/2,1/2,1/2],[0,2/3,1/3],[1/2,1/2,1/2],[-1/2,5/6,1/6],[0,2/3,1/3],[1/2,1/2,1/2],[0,0,0],[-1/2,1/6,-1/6],[0,0,0],[-1,1/3,-1/3],[-1/2,1/6,-1/6],[0,0,0],[-3/2,1/2,-1/2],[-1,1/3,-1/3],[-1/2,1/6,-1/6],[0,0,0]],[[3/2,-3/2,-3/2],[1,-1,-1],[1,-1,-1],[1/2,-5/6,-7/6],[1/2,-1/2,-1/2],[1/2,-1/2,-1/2],[0,-1/3,-2/3],[1/2,-1/2,-1/2],[0,-1/3,-2/3],[-1/2,-1/6,-5/6],[0,0,0],[0,0,0],[-1/2,1/6,-1/6],[0,0,0],[-1/2,1/6,-1/6],[-1,1/3,-1/3],[0,0,0],[-1/2,1/6,-1/6],[-1,1/3,-1/3],[-3/2,1/2,-1/2]],[[0,0,0],[1/2,1/2,1/2],[-1/2,-1/6,1/6],[0,0,0],[1,1,1],[0,1/3,2/3],[1/2,1/2,1/2],[-1,-1/3,1/3],[-1/2,-1/6,1/6],[0,0,0],[3/2,3/2,3/2],[1/2,5/6,7/6],[1,1,1],[-1/2,1/6,5/6],[0,1/3,2/3],[1/2,1/2,1/2],[-3/2,-1/2,1/2],[-1,-1/3,1/3],[-1/2,-1/6,1/6],[0,0,0]],[[0,0,0],[1/2,-1/2,-1/2],[0,0,0],[-1/2,-1/6,1/6],[1,-1,-1],[1/2,-1/2,-1/2],[0,-2/3,-1/3],[0,0,0],[-1/2,-1/6,1/6],[-1,-1/3,1/3],[3/2,-3/2,-3/2],[1,-1,-1],[1/2,-7/6,-5/6],[1/2,-1/2,-1/2],[0,-2/3,-1/3],[-1/2,-5/6,-1/6],[0,0,0],[-1/2,-1/6,1/6],[-1,-1/3,1/3],[-3/2,-1/2,1/2]],[[0,0,0],[0,0,0],[1/2,-1/2,-1/2],[1/2,1/2,1/2],[0,0,0],[1/2,-1/2,-1/2],[1/2,1/2,1/2],[1,-1,-1],[1,0,0],[1,1,1],[0,0,0],[1/2,-1/2,-1/2],[1/2,1/2,1/2],[1,-1,-1],[1,0,0],[1,1,1],[3/2,-3/2,-3/2],[3/2,-1/2,-1/2],[3/2,1/2,1/2],[3/2,3/2,3/2]]], N.float64)
onef_ffs_R2_physvals = N.array([[[0,0,0],[-2/3,-2/3,-2/3],[-1/3,-1/9,1/9],[0,0,0],[-2/3,-2/3,-2/3],[-2/3,-1/3,-1/3],[-1/3,-1/3,-1/3],[-1/3,-1/9,1/9],[-1/6,-1/18,1/18],[0,0,0],[0,0,0],[-1/3,1/9,-1/9],[0,0,0],[-1/3,1/9,-1/9],[-1/6,1/18,-1/18],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[-2/3,2/3,2/3],[0,0,0],[-1/3,-1/9,1/9],[-2/3,2/3,2/3],[-1/3,1/3,1/3],[-2/3,1/3,1/3],[0,0,0],[-1/6,-1/18,1/18],[-1/3,-1/9,1/9],[0,0,0],[0,0,0],[-1/3,1/9,-1/9],[0,0,0],[-1/6,1/18,-1/18],[-1/3,1/9,-1/9],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[-2/3,2/3,2/3],[1/3,1/3,1/3],[0,0,0],[-1/3,1/3,1/3],[1/6,1/6,1/6],[-2/3,2/3,2/3],[-1/3,5/9,4/9],[1/3,1/3,1/3],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[-1/6,1/18,-1/18],[0,0,0],[0,0,0],[-1/3,1/9,-1/9],[-1/3,1/9,-1/9],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[-1/3,1/3,1/3],[1/6,1/6,1/6],[0,0,0],[-1/6,-1/18,1/18],[0,0,0],[0,0,0],[-2/3,2/3,2/3],[1/3,1/3,1/3],[-2/3,2/3,2/3],[-1/3,4/9,5/9],[1/3,1/3,1/3],[0,0,0],[-1/3,-1/9,1/9],[-1/3,-1/9,1/9],[0,0,0]],[[0,0,0],[1/3,1/3,1/3],[-1/3,-1/9,1/9],[0,0,0],[1/3,1/3,1/3],[1/3,0,1/3],[1/6,1/6,1/6],[-1/3,-1/9,1/9],[-1/6,-1/18,1/18],[0,0,0],[0,0,0],[2/3,-2/9,2/9],[0,0,0],[2/3,-2/9,2/9],[1/3,-1/9,1/9],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[1/3,-1/3,-1/3],[0,0,0],[-1/3,-1/9,1/9],[1/3,-1/3,-1/3],[1/6,-1/6,-1/6],[1/3,-1/3,0],[0,0,0],[-1/6,-1/18,1/18],[-1/3,-1/9,1/9],[0,0,0],[0,0,0],[2/3,-2/9,2/9],[0,0,0],[1/3,-1/9,1/9],[2/3,-2/9,2/9],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,0,0],[1/3,-1/3,-1/3],[1/3,1/3,1/3],[0,0,0],[1/6,-1/6,-1/6],[1/6,1/6,1/6],[1/3,-1/3,-1/3],[2/3,-1/9,1/9],[1/3,1/3,1/3],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[1/3,-1/9,1/9],[0,0,0],[0,0,0],[2/3,-2/9,2/9],[2/3,-2/9,2/9],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[1/6,-1/6,-1/6],[1/6,1/6,1/6],[0,0,0],[1/3,1/9,-1/9],[0,0,0],[0,0,0],[1/3,-1/3,-1/3],[1/3,1/3,1/3],[1/3,-1/3,-1/3],[2/3,1/9,-1/9],[1/3,1/3,1/3],[0,0,0],[2/3,2/9,-2/9],[2/3,2/9,-2/9],[0,0,0]]], N.float64)

onef_ffs_R2_physvals_D = N.array([[[9/2,-9,9/2],[3/2,-15/2,6],[3,-6,3],[3,-6,3],[-3/2,-6,15/2],[0,-9/2,9/2],[0,-9/2,9/2],[3/2,-3,3/2],[3/2,-3,3/2],[3/2,-3,3/2],[-9/2,-9/2,9],[-3,-3,6],[-3,-3,6],[-3/2,-3/2,3],[-3/2,-3/2,3],[-3/2,-3/2,3],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[-9/2,9/2,-9],[-3/2,6,-15/2],[-3,3,-6],[-3,3,-6],[3/2,15/2,-6],[0,9/2,-9/2],[0,9/2,-9/2],[-3/2,3/2,-3],[-3/2,3/2,-3],[-3/2,3/2,-3],[9/2,9,-9/2],[3,6,-3],[3,6,-3],[3/2,3,-3/2],[3/2,3,-3/2],[3/2,3,-3/2],[0,0,0],[0,0,0],[0,0,0],[0,0,0]],[[0,-27/2,27/2],[0,-9,9],[3/2,-6,15/2],[0,-9,9],[0,-9/2,9/2],[3/2,-3/2,3],[0,-9/2,9/2],[3,3/2,3/2],[3/2,-3/2,3],[0,-9/2,9/2],[0,0,0],[3/2,3,-3/2],[0,0,0],[3,6,-3],[3/2,3,-3/2],[0,0,0],[9/2,9,-9/2],[3,6,-3],[3/2,3,-3/2],[0,0,0]],[[0,0,0],[0,-9/2,9/2],[-3/2,3/2,-3],[0,0,0],[0,-9,9],[-3/2,-3,3/2],[0,-9/2,9/2],[-3,3,-6],[-3/2,3/2,-3],[0,0,0],[0,-27/2,27/2],[-3/2,-15/2,6],[0,-9,9],[-3,-3/2,-3/2],[-3/2,-3,3/2],[0,-9/2,9/2],[-9/2,9/2,-9],[-3,3,-6],[-3/2,3/2,-3],[0,0,0]],[[0,0,0],[3/2,3/2,-3],[0,3/2,3/2],[0,0,0],[3,3,-6],[3/2,3,-3/2],[3/2,3/2,-3],[0,3,3],[0,3/2,3/2],[0,0,0],[9/2,9/2,-9],[3,9/2,-9/2],[3,3,-6],[3/2,9/2,0],[3/2,3,-3/2],[3/2,3/2,-3],[0,9/2,9/2],[0,3,3],[0,3/2,3/2],[0,0,0]],[[0,0,0],[-3/2,-3,3/2],[0,0,0],[0,3/2,3/2],[-3,-6,3],[-3/2,-3,3/2],[-3/2,-3/2,3],[0,0,0],[0,3/2,3/2],[0,3,3],[-9/2,-9,9/2],[-3,-6,3],[-3,-9/2,9/2],[-3/2,-3,3/2],[-3/2,-3/2,3],[-3/2,0,9/2],[0,0,0],[0,3/2,3/2],[0,3,3],[0,9/2,9/2]],[[0,0,0],[0,0,0],[-3/2,-3,3/2],[3/2,3/2,-3],[0,0,0],[-3/2,-3,3/2],[3/2,3/2,-3],[-3,-6,3],[0,-3/2,-3/2],[3,3,-6],[0,0,0],[-3/2,-3,3/2],[3/2,3/2,-3],[-3,-6,3],[0,-3/2,-3/2],[3,3,-6],[-9/2,-9,9/2],[-3/2,-9/2,0],[3/2,0,-9/2],[9/2,9/2,-9]],[[0,0,0],[0,0,0],[3/2,-3/2,3],[-3/2,3,-3/2],[0,0,0],[3/2,-3/2,3],[-3/2,3,-3/2],[3,-3,6],[0,3/2,3/2],[-3,6,-3],[0,0,0],[3/2,-3/2,3],[-3/2,3,-3/2],[3,-3,6],[0,3/2,3/2],[-3,6,-3],[9/2,-9/2,9],[3/2,0,9/2],[-3/2,9/2,0],[-9/2,9,-9/2]]], N.float64)

twoform_face_vals0 = N.array([[[-3,6,-3],[-3,3,0],[-2,5,-1],[-2,4,-2],[-3,0,3],[-2,2,2],[-2,1,1],[-1,4,1],[-1,3,0],[-1,2,-1],[-3,-3,6],[-2,-1,5],[-2,-2,4],[-1,1,4],[-1,0,3],[-1,-1,2],[0,3,3],[0,2,2],[0,1,1],[0,0,0]],[[3,-3,6],[3,0,3],[2,-2,4],[2,-1,5],[3,3,0],[2,1,1],[2,2,2],[1,-1,2],[1,0,3],[1,1,4],[3,6,-3],[2,4,-2],[2,5,-1],[1,2,-1],[1,3,0],[1,4,1],[0,0,0],[0,1,1],[0,2,2],[0,3,3]],[[0,9,-9],[0,6,-6],[1,8,-7],[1,7,-8],[0,3,-3],[1,5,-4],[1,4,-5],[2,7,-5],[2,6,-6],[2,5,-7],[0,0,0],[1,2,-1],[1,1,-2],[2,4,-2],[2,3,-3],[2,2,-4],[3,6,-3],[3,5,-4],[3,4,-5],[3,3,-6]],[[0,0,0],[0,3,-3],[-1,1,-2],[-1,2,-1],[0,6,-6],[-1,4,-5],[-1,5,-4],[-2,2,-4],[-2,3,-3],[-2,4,-2],[0,9,-9],[-1,7,-8],[-1,8,-7],[-2,5,-7],[-2,6,-6],[-2,7,-5],[-3,3,-6],[-3,4,-5],[-3,5,-4],[-3,6,-3]]], N.float64)
