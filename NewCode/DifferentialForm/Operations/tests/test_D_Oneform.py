from __future__ import division
import random
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N

from NewCode.tests.TestMeshes import TwoTets, FlatTet
from NewCode.tests.BrickMeshes import TwoBricks, OneBrick
from NewCode.tests import xfail

from NewCode.DifferentialForm import Discretiser, BrickDiscretiser
from NewCode import Mesh
from NewCode.Meshes import BrickMesh

from NewCode.DifferentialForm.Operations import D_Oneform

class test_Oneform_curlmat_CTLN(TestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)
        self.disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        self.disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2)
        desired_mat=N.zeros((7,9), N.int8)
        row = [1, 1, -1]
        desired_mat[0,[0,3,1]] = row
        desired_mat[1,[0,4,2]] = row
        desired_mat[2,[1,5,2]] = row
        desired_mat[3,[3,5,4]] = row
        desired_mat[4,[1,7,6]] = row
        desired_mat[5,[2,8,6]] = row
        desired_mat[6,[5,8,7]] = row
        self.desired_mat = desired_mat


    def test_curlmat(self):
        CM = D_Oneform.D_oneform_tet(self.disc1, self.disc2)
        assert_array_equal(self.desired_mat,
                           CM.toarray())

class test_BrickOneform_curlmat_CTLN(TestCase):
    TestMesh = OneBrick
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)
        self.disc1 = BrickDiscretiser.setup_PformDiscretiser(self.mesh, 1)
        self.disc2 = BrickDiscretiser.setup_PformDiscretiser(self.mesh, 2)
        self.desired_mat=N.array([[0,0,0,0,1,-1,0,0,-1,1,0,0],
                                  [0,0,0,0,0,0,1,-1,0,0,-1,1],
                                  [-1,1,0,0,0,0,0,0,1,0,-1,0],
                                  [0,0,-1,1,0,0,0,0,0,1,0,-1],
                                  [1,0,-1,0,-1,0,1,0,0,0,0,0],
                                  [0,1,0,-1,0,-1,0,1,0,0,0,0]],
                                 N.float64)
        

    def test_curlmat(self):
        CM = D_Oneform.D_oneform_brick(self.disc1, self.disc2)
        assert_almost_equal(self.desired_mat, CM.toarray(), decimal=15)

class test_BrickOneform_curlmats(TestCase):
    """
    The approach here is a little different. Calculating all the matrices by
    hand is tedious and unnecesary when we know that the curl of 1-forms should
    be exactly representable by the 2-forms. So we take a 1-form, give it some
    random DOFs, reconstruct it at the several points and then check that the
    values are exactly the same.
    """
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)
        self.test_lams = self.testMesh.test_local_coords

    def test_curl_CTLN1f_LTCN2f(self):
        disc1 = BrickDiscretiser.setup_PformDiscretiser(self.mesh, 1, mixed=True)
        disc2 = BrickDiscretiser.setup_PformDiscretiser(self.mesh, 2, mixed=True)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        tl, l_tl = self.test_lams, len(self.test_lams)
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl, decimal=14)


    def test_curl_mixed_2nd(self):
        disc1 = BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 1, order=2, mixed=True)
        disc2 = BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 2, order=2, mixed=True)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        #tl, l_tl = self.test_lams, len(self.test_lams)
        tl, l_tl = N.array([[1/2, 1/2, 1/2, 1/2, 1/2, 1/2]]), 1 
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl, decimal=13)
        
    def test_curl_mixed_4th(self):
        disc1 = BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 1, order=4, mixed=True)
        disc2 = BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 2, order=4, mixed=True)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        #tl, l_tl = self.test_lams, len(self.test_lams)
        tl, l_tl = N.array([[1/2, 1/2, 1/2, 1/2, 1/2, 1/2]]), 1 
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl, decimal=13)

class test_Oneform_curlmats_base(TestCase):
    """
    The approach here is a little different. Calculating all the matrices by
    hand is tedious and unnecesary when we know that the curl of 1-forms should
    be exactly representable by the 2-forms. So we take a 1-form, give it some
    random DOFs, reconstruct it at the several points and then check that the
    values are exactly the same.
    """
    TestMesh = TwoTets
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = Mesh.Mesh(self.testMesh.listmesh)
        self.test_lams = N.array([(i,j,k,1-i-j-k)
                                  for i in N.arange(0,1.0000001, 1/4.)
                                  for j in N.arange(0,1.0000001-i, 1/4.)
                                  for k in N.arange(0,1.0000001-i-j, 1/4.)
                                  ])

class test_Oneform_curlmats(test_Oneform_curlmats_base):
    def test_curl_CTLN1f_LTCN2f(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1, mixed=True)
        disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2, mixed=True)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        tl, l_tl = self.test_lams, len(self.test_lams)
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl)


    def test_curl_LTQN1f_LTLN2f(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1, order=2)
        disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2, mixed=False)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        tl, l_tl = self.test_lams, len(self.test_lams)
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,1)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl)

class test_Oneform_curlmats_FlatTet(test_Oneform_curlmats_base):
    TestMesh = FlatTet    
    def test_curl_LTQN1f_LTLN2f(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1, order=3)
        disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2, order=2,mixed=False)
        dof1 = disc1.newDOFs()
        dof1_d = disc1.D().newDOFs()
        dof1.dofArray[:] = [random.random() for i in range(len(dof1.dofArray))]
        dof1_d.dofArray[:] = dof1.dofArray
        dof2 = dof1.D(disc2)
        tl, l_tl = self.test_lams, len(self.test_lams)
        onef_curl = N.array([dof1_d.reconstruct([i]*l_tl, tl) for i in (0,)], N.float64)
        twof_curl_rep = N.array([dof2.reconstruct([i]*l_tl, tl) for i in (0,)], N.float64)
        assert_almost_equal(twof_curl_rep, onef_curl)
