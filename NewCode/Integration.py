from __future__ import division

import numpy as N
from scipy.special import orthogonal

from NewCode.Utilities import Struct

import IntegrationRules 



def get_unit_lobatto_rule(no_points):
    x,w,P = gauss_lobatto_coeffs(no_points-1)
    return (x/2 + 0.5, w/2)    

def quad_tet(fun, order):
    """
    Integrate a function of the 4 tetrahedral simplex coordinates over a unit tet

    Arguments
    =========

    fun
      Function that takes an array([l1,l2,l3,l4]), where l1-4 are the simplex
      coordinates and returns a scalar value

    order
      Desired polynomial integration order. Will exactly integrate polynomial
      of this order or less

    Currently only order 2 integration is implemented
    """

    try:
        rules = IntegrationRules.tet_rules[order]
    except KeyError:
        raise KeyError, "Unknown integration order"

    return (N.array([fun(coord) for coord in rules.coords]) * rules.weights).sum()
    
class Integrator(object):
    def __init__(self, order):
        """
        Initialise TetIntegrator with rule of order order
        """
        self.order = order
        try:
            self.rule = self.all_rules[order]
        except KeyError:
            raise KeyError, "Unknown integration order"
        self.noPts = len(self.rule.weights)


    def evalPoints(self):
        """
        Return the evalutation coordinates required for integration
        """
        return self.rule.coords

    def integrateFun(self, funvals):
        """
        Evalute integral for a function with values funvals at evalPoints()

        The user is responsible for retrieving the evalution points using the
        evalPoints() method and evaluating the integrand at these points.

        """
        
        assert len(funvals) == self.noPts
        return N.add.reduce((funvals.T * self.rule.weights).T)

class LineIntegrator(Integrator):
    @staticmethod
    def _get_pts_weights(no_points):
        return tuple(N.real(x) for x in orthogonal.ps_roots(no_points))

    @staticmethod
    def _calc_order(no_points):
        return 2*no_points - 1

    def __init__(self, no_points):
        """
        Initialise LineIntegrator with rule of no_points points

        Accurate to polynomial order 2*no_points -1
        """
        self.no_points = no_points
        pts, weights = self._get_pts_weights(no_points)
        self.rule = IntegrationRules.IntgRule(weights, pts)
        self.noPts = len(self.rule.weights)
        self.order = self._calc_order(no_points)

class TriIntegrator(Integrator):
    all_rules = IntegrationRules.tri_rules

class TetIntegrator(Integrator):
    all_rules = IntegrationRules.tet_rules

class TetProdIntegrator(LineIntegrator):
    def __init__(self, order):
        """
        Initialise tetrahedral product rule integrator with rule of order order
        """
        
        self.order = order 
        no_lin_pts = order // 2 + 1
        lin_pts, lin_wts = self._get_pts_weights(no_lin_pts)
        tri_rule = IntegrationRules.tri_rules[order+1]
        tri_pts, tri_wts = tri_rule.coords, tri_rule.weights
        tet_pts = []; tet_wts = []
        for tr_p, tr_w in zip(tri_pts, tri_wts):
            for l_p, l_w in zip(lin_pts, lin_wts):
                t_xy = tr_p[0:2]
                fac_xy = 1. - N.sum(t_xy)
                t_z = fac_xy*l_p
                t_u1_u3 = N.hstack((t_xy, t_z))
                t_u4 = 1. - N.sum(t_u1_u3)
                tet_pts.append(N.hstack((t_u1_u3, t_u4)))
                tet_wts.append(3.*fac_xy*tr_w*l_w)
                
        self.rule = IntegrationRules.IntgRule(N.array(tet_wts), N.array(tet_pts))
        self.noPts = len(tet_wts)
                
class PyramIntegrator(Integrator):
    all_rules = IntegrationRules.pyr_rules


class HexaIntegrator(LineIntegrator):
    def __init__(self, no_points):
        """
        Initialise HexaIntegrator with rule of no_points^3 points

        Accurate to polynomial order 2*no_points -1  in each axis
        """
        self.no_points = no_points
        pts1d, w1d = self._get_pts_weights(no_points)
        weights = N.array([wx*wy*wz for wx in w1d for wy in w1d for wz in w1d],
                          N.float64)
        pts = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for px in pts1d for py in pts1d for pz in pts1d],
                      N.float64)
        self.rule = IntegrationRules.IntgRule(weights, pts)
        self.noPts = len(self.rule.weights)
        self.order = self._calc_order(no_points)

class HexaLobattoIntegrator(HexaIntegrator):
    def __init__(self, no_points):
        """
        Initialise HexaIntegrator with Gauss-Lobatto rule of no_points^3 points

        Accurate to polynomial order 2*(no_points-1) - 1  in each axis.

        Equivalent to BrickTrapzIntegrator for no_points=2
        """
        HexaIntegrator.__init__(self, no_points)
        
    @staticmethod
    def _calc_order(no_points):
        return 2*no_points - 3

    _get_pts_weights = staticmethod(get_unit_lobatto_rule)

class BrickTrapzIntegrator(Integrator):
    def __init__(self, no_points):
        """
        Use trapezoidal rule, i.e. evaluate functions at the vertices of the brick

        no_points is blithely ignored.
        """
        self.order = None
        pts1d, w1d = N.array([0,1], N.float64), N.array([1/2., 1/2.], N.float64)
        weights = N.array([wx*wy*wz for wx in w1d for wy in w1d for wz in w1d],
                          N.float64)
        pts = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for px in pts1d for py in pts1d for pz in pts1d],
                      N.float64)
        self.rule = IntegrationRules.IntgRule(weights, pts)
        self.noPts = len(self.rule.weights)
        

class HexaGaussLobattoLobattoIntegrator(Integrator):
    """See fisher05 eqn (15) or cohen98 eqn (3.4)"""
    def __init__(self, no_gauss_pts):
        no_lobatto_pts = no_gauss_pts + 1
        g_pts, g_wts = (N.real(x) for x in orthogonal.ps_roots(no_gauss_pts))
        l_pts, l_wts = get_unit_lobatto_rule(no_lobatto_pts)
        weights = N.array([wx*wy*wz for wx in g_wts for wy in l_wts for wz in l_wts],
                          N.float64)
        pts_xx = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for px in g_pts for py in l_pts for pz in l_pts],
                      N.float64)
        pts_yy = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for py in g_pts for pz in l_pts for px in l_pts],
                      N.float64)
        pts_zz = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for pz in g_pts for px in l_pts for py in l_pts],
                      N.float64)
        self.noPts = len(weights)
        self.rule = Struct(weights=weights,coords=Struct(
            xx=pts_xx, yy=pts_yy, zz=pts_zz))
        self.order = None               # Not quite sure how to handle this
        
class HexaLobattoGaussGaussIntegrator(Integrator):
    """See fisher05 eqn (15) or cohen98 eqn (3.4)"""
    def __init__(self, no_gauss_pts):
        no_lobatto_pts = no_gauss_pts + 1
        g_pts, g_wts = (N.real(x) for x in orthogonal.ps_roots(no_gauss_pts))
        l_pts, l_wts = get_unit_lobatto_rule(no_lobatto_pts)
        weights = N.array([wx*wy*wz for wx in l_wts for wy in g_wts for wz in g_wts],
                          N.float64)
        pts_xx = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for px in l_pts for py in g_pts for pz in g_pts],
                      N.float64)
        pts_yy = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for py in l_pts for pz in g_pts for px in g_pts],
                      N.float64)
        pts_zz = N.array([(px,py,pz,1-px,1-py,1-pz)
                       for pz in l_pts for px in g_pts for py in g_pts],
                      N.float64)
        self.noPts = len(weights)
        self.rule = Struct(weights=weights,coords=Struct(
            xx=pts_xx, yy=pts_yy, zz=pts_zz))
        self.order = None               # Not quite sure how to handle this

class QuadIntegrator(Integrator):
    def __init__(self, no_points):
        """
        Initialise QuadIntegrator with rule of no_points^2 points

        Accurate to polynomial order 2*no_points -1  in each axis
        """
        self.no_points = no_points
        pts1d, w1d = tuple(N.real(x) for x in orthogonal.ps_roots(no_points))
        weights = N.array([wx*wy for wx in w1d for wy in w1d],
                          N.float64)
        pts = N.array([(px,py) for px in pts1d for py in pts1d],
                      N.float64)
        self.rule = IntegrationRules.IntgRule(weights, pts)
        self.noPts = len(self.rule.weights)
    
def gauss_lobatto_coeffs(n):
    """
    Computes the Legendre-Gauss-Lobatto nodes, weights and Vandermonde matrix.

    Input
    =====
    n -- Computes an n+1 pt rule of order (2*n - 1)

    Output
    ======
    tuple (x,w,P) with:

    x -- Integration points
    w -- Integration Weights
    P -- Legendre-Gauss-Lobatto Vandermonde matrix

    The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
    integration and spectral methods. 

    based on lglnodes.m written by Greg von Winckel - 04/17/2004
    gregvw@chtm.unm.edu

    Reference used by Greg:
    C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, \"Spectral Methods
    in Fluid Dynamics,\" Section 2.3. Springer-Verlag 1987
    """
    
    eps = N.finfo(float).eps
    # Truncation + 1
    n1 = n+1;
    #Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = N.cos(N.pi*N.arange(0,n+1)/n)

    # The Legendre Vandermonde Matrix
    P = N.zeros((n1,n1));

    # Compute P_(n) using the recursion relation
    # Compute its first and second derivatives and 
    # update x using the Newton-Raphson method.

    xold = 2

    while N.max(N.abs(x-xold)) > eps:
        xold = x
        P[:,0] = 1;    P[:,1] = x;
    
        for k in range(1,n):
            P[:,k+1]=( (2*k+1)*x*P[:,k] - k*P[:,k-1] )/(k+1)
     
        x = xold - ( x*P[:,n1-1] - P[:,n-1] ) / ( n1*P[:,n1-1] )
             

    w=2/(n*n1*P[:,n1-1]**2)

    return (x[::-1],w,P)

