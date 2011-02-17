from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.tests.TestMeshes import FlatTet, TwoTets, InscribedTetMesh
from NewCode.tests import xfail
from NewCode.Utilities import Struct

from NewCode import Mesh
from NewCode.DifferentialForm import Discretiser
from NewCode.DifferentialForm import DiscretiserDOFs, constrained_on_boundary

def test_PformDiscretiserDOFs_init():
    disc = Struct(newDOFsArray=lambda x: x, matrix=None)
    inst = DiscretiserDOFs.PformDiscretiserDOFs(disc, N.float64)
    assert(inst.disc is disc)
    assert_equal(inst.dofArray, N.float64)
    assert_equal(inst.dtype, N.float64)

class test_PformDiscretiserRHS(NumpyTestCase):
    def setUp(self):
        mesh = Mesh.Mesh(TwoTets.listmesh)
        self.disc = Discretiser.setup_PformDiscretiser(mesh, form=1)
        self.inst = DiscretiserDOFs.PformDiscretiserRHS(self.disc, N.float64)
        
    def test_calcProjPointfunRHS(self):
        desired_dofs = N.array([0,0,1/3,0,0,1/6,-1/3,-1/6,-1/2])
        r0 = N.array([-1/3, -1/6, -1/3])
        matchfun = lambda r: N.array([0,0,1.])
        assert_almost_equal(self.inst.calcProjPointfunRHS(matchfun, r0),
                            desired_dofs, decimal=15)
        
class test_PformDiscretiserDOFs(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)

    def test_D_Oneform(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2)
        disc1_dofs = disc1.newDOFs()
        disc1_dofs.dofArray[:] = [-41533/500000,-354693/1000000,34451/500000,-15471/500000,-444307/1000000,11483/50000,-213497/500000,51111/500000,176393/500000]
        disc2_dofs = disc1_dofs.D(disc2) # Curl of one-form represented by the 2-form
        assert_almost_equal(disc2_dofs.dofArray,
                            [48137/200000, -23851/40000, -38787/200000, 25721/40000, 174523/1000000, 424341/500000, 15007/31250], 
                            decimal=15)

    def test_D_Oneform_disc(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        disc1_dofs = disc1.newDOFs()
        disc1_D = disc1.D()
        disc1_D_dofs = disc1_D.newDOFs()
        assert_equal(disc1_dofs.matrix.stiffness().todense(),
                     disc1_D_dofs.matrix.mass().todense())
        self.assertRaises(ValueError, disc1_D_dofs.matrix.stiffness)
        
class test_PformDiscretiserDOFs_hodge(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        self.disc1 = Discretiser.setup_PformDiscretiser(
            self.mesh, 1, freeFun=constrained_on_boundary)
        self.disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2)
        self.disc2_dofs = self.disc2.newDOFs()
        # Calculated on dark and stormy night. For a constant [300, 2, -1] field.
        # Based on some non-cannonical numeric results, not to be 100% trusted.
        self.disc2_dofs.dofArray[:] = N.array([-49.833333333332668,         
                                               50.166666666665982,          
                                               50.166666666666003,          
                                               -49.833333333332668,         
                                               -50.499999999999325,         
                                               -50.166666666665989,         
                                               16.722222222221994,          
                                               16.388888888888673,          
                                               -50.166666666665996,         
                                               -49.499999999999332,         
                                               17.05555555555534,           
                                               16.388888888888655,          
                                               50.499999999999325,          
                                               49.833333333332661,          
                                               -16.944444444444219,         
                                               -16.277777777777562,         
                                               50.499999999999318,          
                                               -49.499999999999318,         
                                               -49.944444444443775,         
                                               50.055555555554882,          
                                               49.833333333332675,          
                                               49.499999999999346,          
                                               -16.944444444444219,         
                                               -16.611111111110887,         
                                               -16.833333333333105,         
                                               16.499999999999776,          
                                               -16.61111111111089,          
                                               16.722222222222001], N.float64)
    def test_hodgeTwoform_value(self):
        #self.disc2_dofs.useLU = True
        actual_disc1_dofs = self.disc2_dofs.hodgeStar(self.disc1)
        # See dark stormy night comment
        desired_1fdofs = N.array([-0.41994750656168195,
                                  -126.82414698162638, 
                                  125.56430446194133,  
                                  -126.40419947506467, 
                                  -1.2598425196850382, 
                                  125.14435695537968], N.float64)
        
        assert_almost_equal(actual_disc1_dofs.dofArray, desired_1fdofs, decimal=12)

    def test_hodgeChecking(self):
        self.assertRaises(TypeError, self.disc2_dofs.hodgeStar, self.disc2)
        disc1_dofs = self.disc1.newDOFs()
        self.assertRaises(TypeError, disc1_dofs.hodgeStar, self.disc1)
        disc1_wrongmesh = Discretiser.setup_PformDiscretiser(
            Mesh.Mesh(TwoTets.listmesh), 1)
        self.assertRaises(
            AssertionError, self.disc2_dofs.hodgeStar, disc1_wrongmesh)

class test_Reconstruct(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(FlatTet.listmesh)

    def test_OneformReconstruct(self):
        disc1 = Discretiser.setup_PformDiscretiser(self.mesh, 1)
        disc1_dofs = disc1.newDOFs()
        disc1_dofs.dofArray[:] = N.array(
            [-3.,   98.66666666666666666667,
             98.33333333333333333,  101.66666666666666666666667,
             101.3333333333333333,   -0.333333333333333333333],
            N.float64)
        assert_almost_equal(disc1_dofs.reconstruct(
                    [0], [N.array([1/4, 1/4, 1/4, 1/4], N.float64)]),
                            [[300., 2., -1.]],
                            decimal=13)

    def test_TwoformReconstruct(self):
        disc2 = Discretiser.setup_PformDiscretiser(self.mesh, 2)
        disc2_dofs = disc2.newDOFs()
        disc2_dofs.dofArray[:] = N.array(
            [-49.83333333333333333,  50.1666666666666666667,
              50.16666666666666667, -49.8333333333333333333], N.float64)
        assert_almost_equal(disc2_dofs.reconstruct(
                    [0], [N.array([1/4, 1/4, 1/4, 1/4], N.float64)]),
                            [[300., 2., -1.]],
                            decimal=13)

        
class test_match_Oneform(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        (a,b,d) = (1.5, 1, 1)
        abd = N.array((a,b,d), N.float64)
        self.matchfun = lambda x: N.array(
            [N.sin(N.pi*(x[0]-a/2)/a)*N.sin(N.pi*(x[2]-d/2)/d), 0, 0], N.float64) 
        self.disc = Discretiser.setup_PformDiscretiser(
            self.mesh, 1, freeFun=constrained_on_boundary)
        self.disc.setIntegrationRule(2) # Desired values calc'd with 2nd order intg
        self.inst = self.disc.newDOFs()
        
    def test_matchRHS(self):
        # Numerical results from development code, take with grain of salt.
        desired_RHS = N.array([-3.7947076036992655e-19,
                               -0.053799608632419949,	  
                               0.064286392425964831,	  
                               -0.064286392425964831,	  
                               -3.7947076036992655e-19, 
                               0.053799608632419949], N.float64)
        assert_almost_equal(self.inst.calcProjRHS(self.matchfun), desired_RHS,
                            decimal=14)
    def test_matchFunction(self):
        # Numerical results from development code, take with grain of salt.
        desired_dofs = N.array([-2.1277467635028025e-18,
                                -0.29413796972261208,
                                0.37532597328553974,
                                -0.37532597328553974,
                              -1.5077186461126546e-19,
                                0.29413796972261208], N.float64)
        self.inst.matchFunction(self.matchfun)
        assert_almost_equal(self.inst.dofArray, desired_dofs,
                            decimal=15)

class test_match_Oneform_complex(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(TwoTets.listmesh)
        self.disc = Discretiser.setup_PformDiscretiser(
            self.mesh, 1)
        self.inst = self.disc.newDOFs(N.complex128)
        self.inst.matrix.dtype = N.complex128
        k = 1.
        self.matchfun = lambda r: N.array([N.e**(-1j*k*r[2]), 0, 0], N.complex128)

    def test_matchRHS(self):
        desired_RHS = 1j*N.array( [0.0e0,2.0366046973002802662113874428144e-2,-8.1855533039967063084777917483089e-3,-8.1855533039967063084777917483089e-3,8.1855533039967063084777917483089e-3,-4.0927766519983531542388958741544e-3,-1.2180493669006096353636082679833e-2,-2.4360987338012192707272165359666e-2,-4.0927766519983531542388958741544e-3]) + N.array( [8.1268515318033284430294287827313e-2,1.2053851409305047559402846644959e-1,0.0e0,0.0e0,-8.1268515318033284430294287827313e-2,-1.2190277297704992664544143174097e-1,-3.9269998775017191163734178622272e-2,-7.8539997550034382327468357244543e-2,-4.0634257659016642215147143913656e-2])
        assert_almost_equal(self.inst.calcProjRHS(self.matchfun), desired_RHS,
                            decimal=11) # Accuracy seems to be limited by integration order

    def test_matchFunction(self):
        desired_dofs = N.array([9.7886020750706461596729936137146e-1-3.2481316450682923609696220479557e-2*1j,3.0430400045409663304121921195249e-1*1j+9.5723417779033257707823457823395e-1,4.6485858269610924714812150669344e-3-1.8520584013492591313899973686088e-1*1j,3.6380236906652028037679074436942e-3-2.6167680896259070024707438943204e-1*1j,2.291954925119077766373781689528e-1*1j-9.7522218381639941316353145392776e-1,-3.3336808899942985316354276753774e-4*1j-9.798707696433605056350126689947e-1,8.1855533039967063084777917483089e-3-4.8526302102046459196383120834071e-2*1j,-3.4085815015819277168364910515011e-1*1j-9.3429441729641588162114249518622e-1,1.4557890630613937758914936250221e-1*1j-2.4556659911990118925433375244914e-2], N.complex128)
        self.inst.matchFunction(self.matchfun)
        assert_almost_equal(self.inst.dofArray, desired_dofs,
                            decimal=8)
    def test_matchErrRMS(self):
        # Note, we need to make a test-case for the new RMS error calculation
        self.inst.dofArray[:] =  [8j+1,7j+2,6j+3,5j+4,4j+5,3j+6,2j+7,1j+8,9]
        # Calculated analytically in Maxima
        desired =  9.040726374072066
        # Note the accuracy is heavily impacted by the order of integration
        assert_almost_equal(self.inst.matchErrRMS(self.matchfun), desired,
                            decimal=8)

class test_match_Twoform(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(InscribedTetMesh.listmesh)
        (a,b,d) = (1.5, 1, 1)
        abd = N.array((a,b,d), N.float64)
        self.matchfun = lambda x: N.array(
            [N.sin(N.pi*(x[0]-a/2)/a)*N.sin(N.pi*(x[2]-d/2)/d), 0, 0], N.float64) 
        self.disc = Discretiser.setup_PformDiscretiser(
            self.mesh, 2, freeFun=constrained_on_boundary)
        self.disc.setIntegrationRule(2) # Desired values calc'd with 2nd order intg
        self.inst = self.disc.newDOFs()
        
    def test_matchRHS(self):
        # Numerical results from development code, take with grain of salt.
        desired_RHS = N.array([0.12051002293748433,
                               -0.12051002293748435,
                               0.083189098900329661,
                               0.083189098900329661,
                               0.1422159079312362,
                               0.14221590793123617,
                               -0.1422159079312362,
                               -0.14221590793123617,
                               -0.12051002293748435,
                               0.12051002293748432,
                               -0.083189098900329661,
                               -0.083189098900329661,
                               -0.14347119036132294,
                               0.14347119036132289,
                               -0.14347119036132291,
                               0.14347119036132294], N.float64)
        assert_almost_equal(self.inst.calcProjRHS(self.matchfun), desired_RHS,
                            decimal=13)
    def test_matchFunction(self):
        # Numerical results from development code, take with grain of salt.
        desired_dofs = N.array([0.096288362037501476,
                                -0.096288362037501476,
                                0.029771756494969139,
                                0.029771756494969142,
                                0.046636559075227707,
                                0.046636559075227693,
                                -0.0466365590752277,
                                -0.046636559075227693,
                                -0.096288362037501476,
                                0.096288362037501476,
                                -0.029771756494969142,
                                -0.029771756494969145,
                                -0.051514405319910458,
                                0.051514405319910431,
                                -0.051514405319910458,
                                0.051514405319910465], N.float64)
        self.inst.matchFunction(self.matchfun)
        assert_almost_equal(self.inst.dofArray, desired_dofs,
                            decimal=14)
