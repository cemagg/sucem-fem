from __future__ import division

from numpy.testing import NumpyTestCase, assert_array_equal,\
     assert_array_almost_equal, assert_almost_equal, assert_equal
import numpy as N
#
# Local Imports
#
from NewCode import Coordinates, ProxyList
from NewCode.tests.TestMeshes import FlatTet
from NewCode.tests.BrickMeshes import OneBrick
from NewCode.tests.PyramMeshes import SixPyram

class test_SimplexCoord(NumpyTestCase, FlatTet):
    def setUp(self):
        nodeCoords = self.listmesh['Nodes'] 
        self.instance = Coordinates.SimplexCoord(nodeCoords=nodeCoords)
        self.instance.LOCAL_FACENODES = N.array([[1,2,3],
                                               [1,2,4],
                                               [1,3,4],
                                               [2,3,4]], N.int32) - 1
        self.instance.LOCAL_FACE_OPPOSING_NODE = N.array([4,3,2,1], N.int32) - 1
        
    def test_covBaseVecs(self):
        assert_array_almost_equal(self.instance.covBaseVecs().transpose(),
                                  [[-1.5, 0.5,-0.5], # Calculated with Maxima
                                   [-1.5,-0.5, 0.5],
                                   [ 1.5, 1.5, 1.5],
                                   [ 1.5,-1.5,-1.5]],
                                  decimal=15)
        
    def test_FaceBasisVecs(self):
        # Calculated with Maxima 
        facebases = [[[-1.5, 3.0,-1.5],[-1.5,-1.5, 3.0],[ 0.0, 1.5, 1.5]],
                     [[ 1.5,-1.5, 3.0],[ 1.5, 3.0,-1.5],[ 0.0, 1.5, 1.5]],
                     [[ 0.0, 4.5,-4.5],[ 1.5, 3.0,-1.5],[ 1.5, 1.5,-3.0]],
                     [[ 0.0, 4.5,-4.5],[-1.5, 1.5,-3.0],[-1.5, 3.0,-1.5]]]

        assert_array_almost_equal([self.instance.FaceBasisVecs(i)
                                   for i in range(4)],
                                  facebases, decimal=15)

    def test_conBaseVecs(self):
        desired = N.array([[0,3/2,3/2],[3/2,3/2,-3],[-3/2,-3,3/2],
                         [-3/2,3,-3/2],[3/2,-3/2,3],[0,9/2,-9/2]],
                        N.float64).transpose()
        assert_array_almost_equal(self.instance.conBaseVecs(), desired)
                   
    def test_local2global2local(self):
        conv_globals = [self.instance.local2global(local_coord)
                        for local_coord in self.test_local_coords]
        assert_almost_equal([self.instance.global2local(global_coord)
                             for global_coord in conv_globals
                             ],
                            self.test_local_coords, decimal=15)

    def test_local2global(self):
        el = self.instance
        assert_almost_equal([el.local2global(lam)
                             for lam in self.test_local_coords],
                            self.test_xyz_coords)


    def test_face_coords2vol_coords(self):
        el = self.instance
        desired = N.array([[ 1.,  2.,  3.,  0.],
                         [ 1.,  2.,  0.,  3.],
                         [ 1.,  0.,  2.,  3.],
                         [ 0.,  1.,  2.,  3.]])
        assert_array_equal([
            el.face_coords2vol_coords(i, N.array([1,2,3.])) for i in range(4)],
            desired)

    def test_vol_coords2face_coords(self):
        el = self.instance
        test_inputs = N.array([[ 1.,  2.,  3.,  0.],
                               [ 1.,  2.,  0.,  3.],
                               [ 1.,  0.,  2.,  3.],
                               [ 0.,  1.,  2.,  3.]])
        for i, ti in enumerate(test_inputs):
            assert_equal(el.vol_coords2face_coords(i, ti), N.array([1,2,3.]))
                
    def test_face_coords2vol_coords2face_coords(self):
        for i in range(4):
            el = self.instance
            assert_equal(el.vol_coords2face_coords(i, el.face_coords2vol_coords(
                i, N.array([1,2,3.]))),
                         N.array([1,2,3.]))

    def test_InElement(self):
        # This test is far from rigourous, but the results seemed OK to me :)
        # At least we'll know when the results change. Ideally a proper
        # analytical calculation for points inside/outside the triangle should
        # be done for this testcase.

        # While the coords no longer sum to 1, naively converting them to
        # global will still yield some points inside the element
        el = self.instance
        local_coords = self.test_local_coords + 0.001
        test_globals = N.array([el.local2global(l)
                      for l in local_coords], N.float64)
        inel = [i for i, xyz in enumerate(test_globals)
                if el.InElement(xyz)]
        outel = [i for i, xyz in enumerate(test_globals)
                 if not el.InElement(xyz)]
        assert_equal(inel, [8, 14, 17, 18])
        assert_equal(outel, [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 15,
                             16, 19])

class test_BrickCoord(NumpyTestCase):
    TestMesh = OneBrick
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        nodeCoords = self.testMesh.nodes
        self.instance = Coordinates.BrickCoord(nodeCoords=nodeCoords)
        self.instance.gridStepSize = self.testMesh.listmesh['GridStepSize']

    def test_covBaseVecs(self):
        self.instance.gridStepSize = N.array([1,2,3], N.float64)
        assert_array_almost_equal(self.instance.covBaseVecs(),
                                  [[1,  0,  0],
                                   [0, 1/2, 0],
                                   [0,  0, 1/3]], decimal=15)
        
    def test_conBaseVecs(self):
        self.instance.gridStepSize = N.array([1,2,3], N.float64)
        assert_array_almost_equal(self.instance.conBaseVecs(), [[1/2/3,0,0],
                                                                [0,1/3,0],
                                                                [0,0,1/2]],
                                  decimal=15)
                   
    def test_J(self):
        self.instance.gridStepSize = N.array([1,2,3], N.float64)
        assert_array_almost_equal(self.instance.J(),
                                  [[1,0,0],
                                   [0,2,0],
                                   [0,0,3]], decimal=15)

    def test_local2global2local(self):
        conv_globals = [self.instance.local2global(local_coord)
                        for local_coord in self.testMesh.test_local_coords]
        assert_almost_equal([self.instance.global2local(global_coord)
                             for global_coord in conv_globals
                             ],
                            self.testMesh.test_local_coords, decimal=15)

    def test_local2global(self):
        el = self.instance
        assert_almost_equal([el.local2global(lam)
                             for lam in self.testMesh.test_local_coords],
                            self.testMesh.test_xyz_coords[0], decimal=15)
    

    def test_face_coords2vol_coords(self):
        el = self.instance
        desired = N.array([[ 0, 1, 2, 1, 0,-1],
                           [ 1, 0, 2, 0, 1,-1],
                           [ 1, 2, 0, 0,-1, 1],

                           [ 1, 1, 2, 0, 0,-1],
                           [ 1, 1, 2, 0, 0,-1],
                           [ 1, 2, 1, 0,-1, 0]], N.float64)
        assert_array_equal([
            el.face_coords2vol_coords(i, N.array([1,2.])) for i in range(6)],
            desired)

    def test_vol_coords2face_coords(self):
        el = self.instance
        test_inputs = N.array([[ 0, 1, 2, 1, 0,-1],
                               [ 1, 0, 2, 0, 1,-1],
                               [ 1, 2, 0, 0,-1, 1],
                               [ 1, 1, 2, 0, 0,-1],
                               [ 1, 1, 2, 0, 0,-1],
                               [ 1, 2, 1, 0,-1, 0]], N.float64)
        assert_array_equal([
            el.vol_coords2face_coords(i, ti) for i,ti in enumerate(test_inputs)],
            N.array([[1,2.]]*len(test_inputs)))

class test_PyramCoord(NumpyTestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        nodeCoords = self.testMesh.elementNodeCoords
        class ElementItemClass(ProxyList.ItemClassFactory(
            ('nodeCoords',), 'ElementItemClass'),
                               Coordinates.PyramCoord): pass
        self.instance = ProxyList.ProxyList(
            ElementItemClass({'nodeCoords':nodeCoords}))

    def test_covBaseVecs(self):
        assert_array_almost_equal(
            [el.covBaseVecs().T for el in self.instance],
            [[[-1,1/2,0],[-1,0,1/3],[2,0,0]],[[1,1/2,0],[1,0,1/3],[-2,0,0]],[[1,-1/2,0],[0,-1/2,1/3],[0,1,0]],[[1,1/2,0],[0,1/2,1/3],[0,-1,0]],[[1,0,-1/3],[0,1/2,-1/3],[0,0,2/3]],[[1,0,1/3],[0,1/2,1/3],[0,0,-2/3]]],
            decimal=15)
        
    def test_conBaseVecs(self):
        assert_array_almost_equal(
            [el.conBaseVecs().T for el in self.instance],
            [[[0,2/3,0],[0,0,1],[1/6,1/3,1/2]],[[0,-2/3,0],[0,0,-1],[1/6,-1/3,-1/2]],[[-1/3,0,0],[0,0,-1],[-1/6,-1/3,-1/2]],[[1/3,0,0],[0,0,1],[1/6,-1/3,1/2]],[[1/3,0,0],[0,2/3,0],[1/6,1/3,1/2]],[[-1/3,0,0],[0,-2/3,0],[-1/6,-1/3,1/2]]], decimal=15)
                   
