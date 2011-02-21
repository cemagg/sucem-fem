from __future__ import division

from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from numpy import array, linspace, ones, float64, zeros, newaxis
from itertools import izip

from NewCode import PostProc, Mesh
from NewCode import DifferentialForm
from NewCode.DifferentialForm import Discretiser
from NewCode.tests.TestMeshes import InscribedTetMesh, TwoTets

class test_LocatePoints(TestCase, InscribedTetMesh):
    postpro_coords = zeros((21,3), float64)
    postpro_coords[:,2] = linspace(-1/2, 1/2, 21)

    def setUp(self):
        self.mesh = Mesh.Mesh(self.listmesh)

    def test_LocatePoints(self):
        (elnos, coords) = PostProc.LocatePoints(self.mesh, self.postpro_coords)
        elements = self.mesh.elements
        for elno, el_coord, global_coord in izip(elnos, coords, self.postpro_coords):
            assert_almost_equal(elements[elno].local2global(el_coord),
                                global_coord, decimal=16)
            
            
class test_ReconstructPoints(TestCase, TwoTets):
    def setUp(self):
        self.mesh = Mesh.Mesh(self.listmesh)
        self.disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        self.disc1_dofs = self.disc1.newDOFs()
        self.disc1_dofs.dofArray += 1
        self.disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2)
        self.disc2_dofs = self.disc2.newDOFs()
        self.disc2_dofs.dofArray += 1
        (nodeA, nodeB) = self.mesh.nodes[[4,1]] # A line from node5 to node1
        delta = nodeB - nodeA
        self.test_coords = nodeA + linspace(0,1,11)[:,newaxis]*delta 
        self.oneform_desired_vals = array(
            [[-1,-1,-1],[-1,-6/5,-4/5],[-1,-7/5,-3/5],[-1,-8/5,-2/5],
             [-1/10,-9/10,1/2],[0,-1,1/2],[1/10,-11/10,1/2],[1/5,-6/5,1/2],
             [3/10,-13/10,1/2],[2/5,-7/5,1/2],[1/2,-3/2,1/2]], float64)
        self.twoform_desired_vals = array(
            [[2,-2,-2],[2,-2,-2],[2,-2,-2],
             [2,-2,-2],[0,0,-2],[0,0,-2],
             [0, 0,-2],[0,0,-2],[0,0,-2],
             [0,0,-2],[0,0,-2]], float64)
            
    def test_ReconstructPoints(self):
        assert_almost_equal(PostProc.ReconstructPoints(self.disc1_dofs, self.test_coords),
                            self.oneform_desired_vals, decimal=15)
        assert_almost_equal(PostProc.ReconstructPoints(self.disc2_dofs, self.test_coords),
                            self.twoform_desired_vals, decimal=15)

class test_misc(TestCase):
    def test_MakeLine(self):
        startpt = array([-0.5, -0.5, -0.75])
        endpt = array([ 0.5,  0.5,  1])
        no_pts = 5
        desired = array([[-0.5   , -0.5   , -0.75  ], 
                         [-0.25  , -0.25  , -0.3125], 
                         [ 0.    ,  0.    ,  0.125 ], 
                         [ 0.25  ,  0.25  ,  0.5625], 
                         [ 0.5   ,  0.5   ,  1.    ]])
        assert_array_equal(PostProc.MakeLine(startpt, endpt, no_pts),
                           desired)
        
