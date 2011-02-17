from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode import Integration

class test_quad_tet(NumpyTestCase):
    def test_quadratic(self):
        def quad_fun(l):
            return l[0] + l[0]**2 + l[1] + l[1]**2 + l[2] + l[2]**2 + 1
        assert_almost_equal(Integration.quad_tet(quad_fun, 2), 41/20,
                            decimal=13)

    def test_ordertest(self):
        self.assertRaises(KeyError, Integration.quad_tet, lambda x: x, 3)
        
class test_TetIntegrator(NumpyTestCase):
    quad_fun = lambda self, l: 1*l[0] + 2*l[0]**2 + 3*l[1] + 4*l[1]**2 + 1 \
               + 5*l[2] + 6*l[2]**2 + 7*l[3] + 8*l[3]**2 \
               + 1*l[0]*l[1] + 2*l[0]*l[2] + 3*l[0]*l[3] \
               - 4*l[1]*l[2] - 5*l[1]*l[3] - 6*l[2]*l[3]

    def setUp(self):
        self.inst2 = Integration.TetIntegrator(2)
        
    def test_eval_pts(self):
        assert_equal(
            N.array(
            [[0.5854101966249685, 0.138196601125015,  0.138196601125015,  0.138196601125015],
             [0.138196601125015,  0.5854101966249685, 0.138196601125015,  0.138196601125015],
             [0.138196601125015,  0.138196601125015,  0.5854101966249685, 0.138196601125015],
             [0.138196601125015,  0.138196601125015,  0.138196601125015,  0.5854101966249685]],
            N.float64),
            self.inst2.evalPoints()
            )

    def test_intg(self):
        funvals = N.array([self.quad_fun(pt) for pt in self.inst2.evalPoints()])
        assert_almost_equal(self.inst2.integrateFun(funvals), 131/20)
        
class test_TetProdIntegrator(NumpyTestCase):
    quad_fun = test_TetIntegrator.quad_fun.im_func

    def test_intg_2(self):
        intg = Integration.TetProdIntegrator(2)
        funvals = N.array([self.quad_fun(pt) for pt in intg.evalPoints()])
        assert_almost_equal(intg.integrateFun(funvals), 131/20)


class _test_PyramIntegrator(NumpyTestCase):
    test_fun = lambda coord: None
    desired_test_intg_val = None
    integrator = Integration.PyramIntegrator
    decimal_tol = 15
    no_pts = None
    
    def test_intg(self):
        inst = self.integrator(self.no_pts)
        fun = self.test_fun
        funvals = N.array([fun(pt) for pt in inst.evalPoints()], N.float64)
        assert_almost_equal(inst.integrateFun(funvals), self.desired_test_intg_val,
                            decimal=self.decimal_tol)

class test_PyramIntegrator_2_2(_test_PyramIntegrator):
    test_fun = staticmethod(lambda lam: (
        5*lam[2]**2*lam[3]**2-6*lam[1]*lam[2]**2*lam[3]+5*lam[1]**2*lam[2]**2)/(4*lam[4]**2-8*lam[4]+4))
    desired_test_intg_val = 21/180
    no_pts = (2,2)
    
class _test_HexaIntegrator(NumpyTestCase):
    third_order_fun = staticmethod(
        lambda x: -3*x[0]**3*x[1]**3*x[2]**3+2*x[0]**2*x[1]**3*x[2]**3-2*x[1]**3*x[2]**3-2*x[0]**2*x[1]**2*x[2]**3+3*x[0]*x[1]**2*x[2]**3+x[1]**2*x[2]**3+3*x[0]**3*x[1]*x[2]**3+x[0]**2*x[1]*x[2]**3-x[0]*x[1]*x[2]**3-3*x[1]*x[2]**3-x[0]**3*x[2]**3-3*x[0]**2*x[2]**3+2*x[0]*x[2]**3+3*x[0]**3*x[1]**3*x[2]**2+x[0]**2*x[1]**3*x[2]**2-x[0]*x[1]**3*x[2]**2-3*x[1]**3*x[2]**2-x[0]**3*x[1]**2*x[2]**2-3*x[0]**2*x[1]**2*x[2]**2+2*x[0]*x[1]**2*x[2]**2+2*x[0]**3*x[1]*x[2]**2-2*x[0]*x[1]*x[2]**2+3*x[1]*x[2]**2-2*x[0]**3*x[2]**2+3*x[0]**2*x[2]**2+x[0]*x[2]**2-x[2]**2+2*x[0]**3*x[1]**3*x[2]-2*x[0]*x[1]**3*x[2]+3*x[1]**3*x[2]-2*x[0]**3*x[1]**2*x[2]+3*x[0]**2*x[1]**2*x[2]+x[0]*x[1]**2*x[2]-x[1]**2*x[2]+x[0]**3*x[1]*x[2]-x[0]**2*x[1]*x[2]-3*x[0]*x[1]*x[2]+2*x[1]*x[2]-3*x[0]**3*x[2]+2*x[0]**2*x[2]-2*x[2]+x[0]**3*x[1]**3-x[0]**2*x[1]**3-3*x[0]*x[1]**3+2*x[1]**3-3*x[0]**3*x[1]**2+2*x[0]**2*x[1]**2-2*x[1]**2-2*x[0]**2*x[1]+3*x[0]*x[1]+x[1]+3*x[0]**3+x[0]**2-x[0]-3)

class test_HexaIntegrator(_test_HexaIntegrator):                            

    def test_rule_2(self):
        inst = Integration.HexaIntegrator(2)
        indep_pts = N.array([[2.11324865405187117745426e-1,2.11324865405187117745426e-1,2.11324865405187117745426e-1],[2.11324865405187117745426e-1,2.11324865405187117745426e-1,7.88675134594812882254574e-1],[2.11324865405187117745426e-1,7.88675134594812882254574e-1,2.11324865405187117745426e-1],[2.11324865405187117745426e-1,7.88675134594812882254574e-1,7.88675134594812882254574e-1],[7.88675134594812882254574e-1,2.11324865405187117745426e-1,2.11324865405187117745426e-1],[7.88675134594812882254574e-1,2.11324865405187117745426e-1,7.88675134594812882254574e-1],[7.88675134594812882254574e-1,7.88675134594812882254574e-1,2.11324865405187117745426e-1],[7.88675134594812882254574e-1,7.88675134594812882254574e-1,7.88675134594812882254574e-1]], N.float64)
        pts = indep_pts[:,[0,1,2,0,1,2]]
        pts[:,3:6] *= -1
        pts[:,3:6] += 1
        assert_almost_equal(pts, inst.evalPoints(), decimal=16)
        assert_almost_equal(N.array([1.25e-1,1.25e-1,1.25e-1,1.25e-1,
                                   1.25e-1,1.25e-1,1.25e-1,1.25e-1], N.float64),
                            inst.rule.weights, decimal=16)

    def test_intg(self):
        inst = Integration.HexaIntegrator(2)
        fun = self.third_order_fun
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])        
        assert_almost_equal(inst.integrateFun(funvals), -551/192, decimal=15)

class test_HexaLobattoIntegrator(_test_HexaIntegrator):
    def test_intg(self):
        inst = Integration.HexaLobattoIntegrator(3)
        fun = self.third_order_fun
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        assert_almost_equal(inst.integrateFun(funvals), -551/192, decimal=15)
        inst = Integration.HexaLobattoIntegrator(4)
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        assert_almost_equal(inst.integrateFun(funvals), -551/192, decimal=14)
        inst = Integration.HexaLobattoIntegrator(5)
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        assert_almost_equal(inst.integrateFun(funvals), -551/192, decimal=14)
        inst = Integration.HexaLobattoIntegrator(6)
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        assert_almost_equal(inst.integrateFun(funvals), -551/192, decimal=14)
        inst = Integration.HexaLobattoIntegrator(2)
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        # 2 point rule has order 1, should not be able to integrate properly
        self.assert_(N.abs(inst.integrateFun(funvals) - -551/192) > .2)
        
class test_QuadIntegrator(NumpyTestCase):
    third_order_fun = staticmethod(
        lambda x: -2*x[0]**3*x[1]**3+x[0]**2*x[1]**3-3*x[0]*x[1]**3-3*x[0]**3*x[1]**2+3*x[0]*x[1]**2-x[1]**2+3*x[0]**3*x[1]-x[0]**2*x[1]+2*x[0]*x[1]-2*x[1]+2*x[0]**3-2*x[0]**2+x[0]-3)
    def test_intg(self):
        inst = Integration.QuadIntegrator(2)
        fun = self.third_order_fun
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])        
        assert_almost_equal(inst.integrateFun(funvals), -83/24, decimal=15)
    

class test_TriIntegrator(NumpyTestCase):
    quad_fun = staticmethod(lambda l: 1*l[0] + 2*l[0]**2 + 3*l[1] + 4*l[1]**2 + 1 \
               + 5*l[2] + 6*l[2]**2 + 1*l[0]*l[1] + 2*l[0]*l[2] \
               - 4*l[1]*l[2]) 

    fourth_order_fun = staticmethod(lambda l:l[2]**4+l[1]*l[2]**3+l[0]*l[2]**3-l[2]**3+l[1]**2*l[2]**2+l[0]*l[1]*l[2]**2-l[1]*l[2]**2+l[0]**2*l[2]**2-l[0]*l[2]**2+l[2]**2+l[1]**3*l[2]+l[0]*l[1]**2*l[2]-l[1]**2*l[2]+l[0]**2*l[1]*l[2]-l[0]*l[1]*l[2]+l[1]*l[2]+l[0]**3*l[2]-l[0]**2*l[2]+l[0]*l[2]-l[2]+l[1]**4+l[0]*l[1]**3-l[1]**3+l[0]**2*l[1]**2-l[0]*l[1]**2+l[1]**2+l[0]**3*l[1]-l[0]**2*l[1]+l[0]*l[1]-l[1]+l[0]**4-l[0]**3+l[0]**2-l[0]+1)

    sixth_order_fun = staticmethod(lambda l: l[2]**6+l[1]*l[2]**5+l[0]*l[2]**5-l[2]**5+l[1]**2*l[2]**4+l[0]*l[1]*l[2]**4-l[1]*l[2]**4+l[0]**2*l[2]**4-l[0]*l[2]**4+l[2]**4+l[1]**3*l[2]**3+l[0]*l[1]**2*l[2]**3-l[1]**2*l[2]**3+l[0]**2*l[1]*l[2]**3-l[0]*l[1]*l[2]**3+l[1]*l[2]**3+l[0]**3*l[2]**3-l[0]**2*l[2]**3+l[0]*l[2]**3-l[2]**3+l[1]**4*l[2]**2+l[0]*l[1]**3*l[2]**2-l[1]**3*l[2]**2+l[0]**2*l[1]**2*l[2]**2-l[0]*l[1]**2*l[2]**2+l[1]**2*l[2]**2+l[0]**3*l[1]*l[2]**2-l[0]**2*l[1]*l[2]**2+l[0]*l[1]*l[2]**2-l[1]*l[2]**2+l[0]**4*l[2]**2-l[0]**3*l[2]**2+l[0]**2*l[2]**2-l[0]*l[2]**2+l[2]**2+l[1]**5*l[2]+l[0]*l[1]**4*l[2]-l[1]**4*l[2]+l[0]**2*l[1]**3*l[2]-l[0]*l[1]**3*l[2]+l[1]**3*l[2]+l[0]**3*l[1]**2*l[2]-l[0]**2*l[1]**2*l[2]+l[0]*l[1]**2*l[2]-l[1]**2*l[2]+l[0]**4*l[1]*l[2]-l[0]**3*l[1]*l[2]+l[0]**2*l[1]*l[2]-l[0]*l[1]*l[2]+l[1]*l[2]+l[0]**5*l[2]-l[0]**4*l[2]+l[0]**3*l[2]-l[0]**2*l[2]+l[0]*l[2]-l[2]+l[1]**6+l[0]*l[1]**5-l[1]**5+l[0]**2*l[1]**4-l[0]*l[1]**4+l[1]**4+l[0]**3*l[1]**3-l[0]**2*l[1]**3+l[0]*l[1]**3-l[1]**3+l[0]**4*l[1]**2-l[0]**3*l[1]**2+l[0]**2*l[1]**2-l[0]*l[1]**2+l[1]**2+l[0]**5*l[1]-l[0]**4*l[1]+l[0]**3*l[1]-l[0]**2*l[1]+l[0]*l[1]-l[1]+l[0]**6-l[0]**5+l[0]**4-l[0]**3+l[0]**2-l[0]+1)

    eighth_order_fun = staticmethod(lambda l:l[2]**8+l[1]*l[2]**7+l[0]*l[2]**7-l[2]**7+l[1]**2*l[2]**6+l[0]*l[1]*l[2]**6-l[1]*l[2]**6+l[0]**2*l[2]**6-l[0]*l[2]**6+l[2]**6+l[1]**3*l[2]**5+l[0]*l[1]**2*l[2]**5-l[1]**2*l[2]**5+l[0]**2*l[1]*l[2]**5-l[0]*l[1]*l[2]**5+l[1]*l[2]**5+l[0]**3*l[2]**5-l[0]**2*l[2]**5+l[0]*l[2]**5-l[2]**5+l[1]**4*l[2]**4+l[0]*l[1]**3*l[2]**4-l[1]**3*l[2]**4+l[0]**2*l[1]**2*l[2]**4-l[0]*l[1]**2*l[2]**4+l[1]**2*l[2]**4+l[0]**3*l[1]*l[2]**4-l[0]**2*l[1]*l[2]**4+l[0]*l[1]*l[2]**4-l[1]*l[2]**4+l[0]**4*l[2]**4-l[0]**3*l[2]**4+l[0]**2*l[2]**4-l[0]*l[2]**4+l[2]**4+l[1]**5*l[2]**3+l[0]*l[1]**4*l[2]**3-l[1]**4*l[2]**3+l[0]**2*l[1]**3*l[2]**3-l[0]*l[1]**3*l[2]**3+l[1]**3*l[2]**3+l[0]**3*l[1]**2*l[2]**3-l[0]**2*l[1]**2*l[2]**3+l[0]*l[1]**2*l[2]**3-l[1]**2*l[2]**3+l[0]**4*l[1]*l[2]**3-l[0]**3*l[1]*l[2]**3+l[0]**2*l[1]*l[2]**3-l[0]*l[1]*l[2]**3+l[1]*l[2]**3+l[0]**5*l[2]**3-l[0]**4*l[2]**3+l[0]**3*l[2]**3-l[0]**2*l[2]**3+l[0]*l[2]**3-l[2]**3+l[1]**6*l[2]**2+l[0]*l[1]**5*l[2]**2-l[1]**5*l[2]**2+l[0]**2*l[1]**4*l[2]**2-l[0]*l[1]**4*l[2]**2+l[1]**4*l[2]**2+l[0]**3*l[1]**3*l[2]**2-l[0]**2*l[1]**3*l[2]**2+l[0]*l[1]**3*l[2]**2-l[1]**3*l[2]**2+l[0]**4*l[1]**2*l[2]**2-l[0]**3*l[1]**2*l[2]**2+l[0]**2*l[1]**2*l[2]**2-l[0]*l[1]**2*l[2]**2+l[1]**2*l[2]**2+l[0]**5*l[1]*l[2]**2-l[0]**4*l[1]*l[2]**2+l[0]**3*l[1]*l[2]**2-l[0]**2*l[1]*l[2]**2+l[0]*l[1]*l[2]**2-l[1]*l[2]**2+l[0]**6*l[2]**2-l[0]**5*l[2]**2+l[0]**4*l[2]**2-l[0]**3*l[2]**2+l[0]**2*l[2]**2-l[0]*l[2]**2+l[2]**2+l[1]**7*l[2]+l[0]*l[1]**6*l[2]-l[1]**6*l[2]+l[0]**2*l[1]**5*l[2]-l[0]*l[1]**5*l[2]+l[1]**5*l[2]+l[0]**3*l[1]**4*l[2]-l[0]**2*l[1]**4*l[2]+l[0]*l[1]**4*l[2]-l[1]**4*l[2]+l[0]**4*l[1]**3*l[2]-l[0]**3*l[1]**3*l[2]+l[0]**2*l[1]**3*l[2]-l[0]*l[1]**3*l[2]+l[1]**3*l[2]+l[0]**5*l[1]**2*l[2]-l[0]**4*l[1]**2*l[2]+l[0]**3*l[1]**2*l[2]-l[0]**2*l[1]**2*l[2]+l[0]*l[1]**2*l[2]-l[1]**2*l[2]+l[0]**6*l[1]*l[2]-l[0]**5*l[1]*l[2]+l[0]**4*l[1]*l[2]-l[0]**3*l[1]*l[2]+l[0]**2*l[1]*l[2]-l[0]*l[1]*l[2]+l[1]*l[2]+l[0]**7*l[2]-l[0]**6*l[2]+l[0]**5*l[2]-l[0]**4*l[2]+l[0]**3*l[2]-l[0]**2*l[2]+l[0]*l[2]-l[2]+l[1]**8+l[0]*l[1]**7-l[1]**7+l[0]**2*l[1]**6-l[0]*l[1]**6+l[1]**6+l[0]**3*l[1]**5-l[0]**2*l[1]**5+l[0]*l[1]**5-l[1]**5+l[0]**4*l[1]**4-l[0]**3*l[1]**4+l[0]**2*l[1]**4-l[0]*l[1]**4+l[1]**4+l[0]**5*l[1]**3-l[0]**4*l[1]**3+l[0]**3*l[1]**3-l[0]**2*l[1]**3+l[0]*l[1]**3-l[1]**3+l[0]**6*l[1]**2-l[0]**5*l[1]**2+l[0]**4*l[1]**2-l[0]**3*l[1]**2+l[0]**2*l[1]**2-l[0]*l[1]**2+l[1]**2+l[0]**7*l[1]-l[0]**6*l[1]+l[0]**5*l[1]-l[0]**4*l[1]+l[0]**3*l[1]-l[0]**2*l[1]+l[0]*l[1]-l[1]+l[0]**8-l[0]**7+l[0]**6-l[0]**5+l[0]**4-l[0]**3+l[0]**2-l[0]+1)

    tenth_order_fun = staticmethod(lambda l: l[2]**10+l[1]*l[2]**9+l[0]*l[2]**9-l[2]**9+l[1]**2*l[2]**8+l[0]*l[1]*l[2]**8-l[1]*l[2]**8+l[0]**2*l[2]**8-l[0]*l[2]**8+l[2]**8+l[1]**3*l[2]**7+l[0]*l[1]**2*l[2]**7-l[1]**2*l[2]**7+l[0]**2*l[1]*l[2]**7-l[0]*l[1]*l[2]**7+l[1]*l[2]**7+l[0]**3*l[2]**7-l[0]**2*l[2]**7+l[0]*l[2]**7-l[2]**7+l[1]**4*l[2]**6+l[0]*l[1]**3*l[2]**6-l[1]**3*l[2]**6+l[0]**2*l[1]**2*l[2]**6-l[0]*l[1]**2*l[2]**6+l[1]**2*l[2]**6+l[0]**3*l[1]*l[2]**6-l[0]**2*l[1]*l[2]**6+l[0]*l[1]*l[2]**6-l[1]*l[2]**6+l[0]**4*l[2]**6-l[0]**3*l[2]**6+l[0]**2*l[2]**6-l[0]*l[2]**6+l[2]**6+l[1]**5*l[2]**5+l[0]*l[1]**4*l[2]**5-l[1]**4*l[2]**5+l[0]**2*l[1]**3*l[2]**5-l[0]*l[1]**3*l[2]**5+l[1]**3*l[2]**5+l[0]**3*l[1]**2*l[2]**5-l[0]**2*l[1]**2*l[2]**5+l[0]*l[1]**2*l[2]**5-l[1]**2*l[2]**5+l[0]**4*l[1]*l[2]**5-l[0]**3*l[1]*l[2]**5+l[0]**2*l[1]*l[2]**5-l[0]*l[1]*l[2]**5+l[1]*l[2]**5+l[0]**5*l[2]**5-l[0]**4*l[2]**5+l[0]**3*l[2]**5-l[0]**2*l[2]**5+l[0]*l[2]**5-l[2]**5+l[1]**6*l[2]**4+l[0]*l[1]**5*l[2]**4-l[1]**5*l[2]**4+l[0]**2*l[1]**4*l[2]**4-l[0]*l[1]**4*l[2]**4+l[1]**4*l[2]**4+l[0]**3*l[1]**3*l[2]**4-l[0]**2*l[1]**3*l[2]**4+l[0]*l[1]**3*l[2]**4-l[1]**3*l[2]**4+l[0]**4*l[1]**2*l[2]**4-l[0]**3*l[1]**2*l[2]**4+l[0]**2*l[1]**2*l[2]**4-l[0]*l[1]**2*l[2]**4+l[1]**2*l[2]**4+l[0]**5*l[1]*l[2]**4-l[0]**4*l[1]*l[2]**4+l[0]**3*l[1]*l[2]**4-l[0]**2*l[1]*l[2]**4+l[0]*l[1]*l[2]**4-l[1]*l[2]**4+l[0]**6*l[2]**4-l[0]**5*l[2]**4+l[0]**4*l[2]**4-l[0]**3*l[2]**4+l[0]**2*l[2]**4-l[0]*l[2]**4+l[2]**4+l[1]**7*l[2]**3+l[0]*l[1]**6*l[2]**3-l[1]**6*l[2]**3+l[0]**2*l[1]**5*l[2]**3-l[0]*l[1]**5*l[2]**3+l[1]**5*l[2]**3+l[0]**3*l[1]**4*l[2]**3-l[0]**2*l[1]**4*l[2]**3+l[0]*l[1]**4*l[2]**3-l[1]**4*l[2]**3+l[0]**4*l[1]**3*l[2]**3-l[0]**3*l[1]**3*l[2]**3+l[0]**2*l[1]**3*l[2]**3-l[0]*l[1]**3*l[2]**3+l[1]**3*l[2]**3+l[0]**5*l[1]**2*l[2]**3-l[0]**4*l[1]**2*l[2]**3+l[0]**3*l[1]**2*l[2]**3-l[0]**2*l[1]**2*l[2]**3+l[0]*l[1]**2*l[2]**3-l[1]**2*l[2]**3+l[0]**6*l[1]*l[2]**3-l[0]**5*l[1]*l[2]**3+l[0]**4*l[1]*l[2]**3-l[0]**3*l[1]*l[2]**3+l[0]**2*l[1]*l[2]**3-l[0]*l[1]*l[2]**3+l[1]*l[2]**3+l[0]**7*l[2]**3-l[0]**6*l[2]**3+l[0]**5*l[2]**3-l[0]**4*l[2]**3+l[0]**3*l[2]**3-l[0]**2*l[2]**3+l[0]*l[2]**3-l[2]**3+l[1]**8*l[2]**2+l[0]*l[1]**7*l[2]**2-l[1]**7*l[2]**2+l[0]**2*l[1]**6*l[2]**2-l[0]*l[1]**6*l[2]**2+l[1]**6*l[2]**2+l[0]**3*l[1]**5*l[2]**2-l[0]**2*l[1]**5*l[2]**2+l[0]*l[1]**5*l[2]**2-l[1]**5*l[2]**2+l[0]**4*l[1]**4*l[2]**2-l[0]**3*l[1]**4*l[2]**2+l[0]**2*l[1]**4*l[2]**2-l[0]*l[1]**4*l[2]**2+l[1]**4*l[2]**2+l[0]**5*l[1]**3*l[2]**2-l[0]**4*l[1]**3*l[2]**2+l[0]**3*l[1]**3*l[2]**2-l[0]**2*l[1]**3*l[2]**2+l[0]*l[1]**3*l[2]**2-l[1]**3*l[2]**2+l[0]**6*l[1]**2*l[2]**2-l[0]**5*l[1]**2*l[2]**2+l[0]**4*l[1]**2*l[2]**2-l[0]**3*l[1]**2*l[2]**2+l[0]**2*l[1]**2*l[2]**2-l[0]*l[1]**2*l[2]**2+l[1]**2*l[2]**2+l[0]**7*l[1]*l[2]**2-l[0]**6*l[1]*l[2]**2+l[0]**5*l[1]*l[2]**2-l[0]**4*l[1]*l[2]**2+l[0]**3*l[1]*l[2]**2-l[0]**2*l[1]*l[2]**2+l[0]*l[1]*l[2]**2-l[1]*l[2]**2+l[0]**8*l[2]**2-l[0]**7*l[2]**2+l[0]**6*l[2]**2-l[0]**5*l[2]**2+l[0]**4*l[2]**2-l[0]**3*l[2]**2+l[0]**2*l[2]**2-l[0]*l[2]**2+l[2]**2+l[1]**9*l[2]+l[0]*l[1]**8*l[2]-l[1]**8*l[2]+l[0]**2*l[1]**7*l[2]-l[0]*l[1]**7*l[2]+l[1]**7*l[2]+l[0]**3*l[1]**6*l[2]-l[0]**2*l[1]**6*l[2]+l[0]*l[1]**6*l[2]-l[1]**6*l[2]+l[0]**4*l[1]**5*l[2]-l[0]**3*l[1]**5*l[2]+l[0]**2*l[1]**5*l[2]-l[0]*l[1]**5*l[2]+l[1]**5*l[2]+l[0]**5*l[1]**4*l[2]-l[0]**4*l[1]**4*l[2]+l[0]**3*l[1]**4*l[2]-l[0]**2*l[1]**4*l[2]+l[0]*l[1]**4*l[2]-l[1]**4*l[2]+l[0]**6*l[1]**3*l[2]-l[0]**5*l[1]**3*l[2]+l[0]**4*l[1]**3*l[2]-l[0]**3*l[1]**3*l[2]+l[0]**2*l[1]**3*l[2]-l[0]*l[1]**3*l[2]+l[1]**3*l[2]+l[0]**7*l[1]**2*l[2]-l[0]**6*l[1]**2*l[2]+l[0]**5*l[1]**2*l[2]-l[0]**4*l[1]**2*l[2]+l[0]**3*l[1]**2*l[2]-l[0]**2*l[1]**2*l[2]+l[0]*l[1]**2*l[2]-l[1]**2*l[2]+l[0]**8*l[1]*l[2]-l[0]**7*l[1]*l[2]+l[0]**6*l[1]*l[2]-l[0]**5*l[1]*l[2]+l[0]**4*l[1]*l[2]-l[0]**3*l[1]*l[2]+l[0]**2*l[1]*l[2]-l[0]*l[1]*l[2]+l[1]*l[2]+l[0]**9*l[2]-l[0]**8*l[2]+l[0]**7*l[2]-l[0]**6*l[2]+l[0]**5*l[2]-l[0]**4*l[2]+l[0]**3*l[2]-l[0]**2*l[2]+l[0]*l[2]-l[2]+l[1]**10+l[0]*l[1]**9-l[1]**9+l[0]**2*l[1]**8-l[0]*l[1]**8+l[1]**8+l[0]**3*l[1]**7-l[0]**2*l[1]**7+l[0]*l[1]**7-l[1]**7+l[0]**4*l[1]**6-l[0]**3*l[1]**6+l[0]**2*l[1]**6-l[0]*l[1]**6+l[1]**6+l[0]**5*l[1]**5-l[0]**4*l[1]**5+l[0]**3*l[1]**5-l[0]**2*l[1]**5+l[0]*l[1]**5-l[1]**5+l[0]**6*l[1]**4-l[0]**5*l[1]**4+l[0]**4*l[1]**4-l[0]**3*l[1]**4+l[0]**2*l[1]**4-l[0]*l[1]**4+l[1]**4+l[0]**7*l[1]**3-l[0]**6*l[1]**3+l[0]**5*l[1]**3-l[0]**4*l[1]**3+l[0]**3*l[1]**3-l[0]**2*l[1]**3+l[0]*l[1]**3-l[1]**3+l[0]**8*l[1]**2-l[0]**7*l[1]**2+l[0]**6*l[1]**2-l[0]**5*l[1]**2+l[0]**4*l[1]**2-l[0]**3*l[1]**2+l[0]**2*l[1]**2-l[0]*l[1]**2+l[1]**2+l[0]**9*l[1]-l[0]**8*l[1]+l[0]**7*l[1]-l[0]**6*l[1]+l[0]**5*l[1]-l[0]**4*l[1]+l[0]**3*l[1]-l[0]**2*l[1]+l[0]*l[1]-l[1]+l[0]**10-l[0]**9+l[0]**8-l[0]**7+l[0]**6-l[0]**5+l[0]**4-l[0]**3+l[0]**2-l[0]+1)
    
    def _test_integ(self, fun, inst, desired, *names, **kwargs):
        funvals = N.array([fun(pt) for pt in inst.evalPoints()])
        assert_almost_equal(inst.integrateFun(funvals), desired, *names, **kwargs)

    def test_quadratic(self):
        self._test_integ(self.quad_fun, Integration.TriIntegrator(2), 71/12, decimal=14)

    def test_fourth_order(self):
        self._test_integ(self.fourth_order_fun, Integration.TriIntegrator(4),
                         7/12, decimal=15)

    def test_sixth_order(self):
        self._test_integ(self.sixth_order_fun, Integration.TriIntegrator(6),
                         517/1008, decimal=14)

    def test_eighth_order(self):
        self._test_integ(self.eighth_order_fun, Integration.TriIntegrator(8),
                         12163/25200, decimal=15)

    def test_tenth_order(self):
        self._test_integ(self.tenth_order_fun, Integration.TriIntegrator(10),
                         6179/13200, decimal=16)

class test_gauss_lobatto_coeffs(NumpyTestCase):
    x_3 = N.array([ -1.0, -0.447213595499958, 0.447213595499958, 1.0], N.float64)
    w_3 = N.array([1/6, 5/6, 5/6, 1/6], N.float64)
    x_5 = N.array([-1.000000000000000,
                   -0.765055323929465,
                   -0.285231516480645,
                   0.285231516480645,
                   0.765055323929465,
                   1.000000000000000], N.float64)
    w_5 = N.array([0.0666666666666667,
                   0.3784749562978469,
                   0.5548583770354865,
                   0.5548583770354865,
                   0.3784749562978469,
                   0.0666666666666667], N.float64)

    def test_3(self):
        x,w,P = Integration.gauss_lobatto_coeffs(3)
        assert_almost_equal(x, self.x_3, decimal=16)
        assert_almost_equal(w, self.w_3, decimal=16)

    def test_5(self):
        x,w,P = Integration.gauss_lobatto_coeffs(5)
        assert_almost_equal(x, self.x_5, decimal=15)
        assert_almost_equal(w, self.w_5, decimal=16)
