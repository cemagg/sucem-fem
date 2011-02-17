from __future__ import division

import numpy as N
from scipy.special import j1, jv
from scipy.integrate import quad, quadrature, romberg

class WaveguideTransientResponse(object):
    epsrel = 2*N.finfo(float).eps
    epsabs = 2*N.finfo(float).eps
    v = 1.                              # Speed of light

    # Impulse response of WG has one part that is a time-delayed dirac-delta and
    # another part involving a finite function. Convolution with the delta part is
    # trivial (i.e. just delayed version of the input), hence only the finite part
    # requires numerical convolution. Implemented from: "Exact, closed-form
    # expressions for transient fields in homogeneously filled waveguides",
    # S. L. Dvorak, IEEE Transactions on Microwave Theory and Techniques vol42, no
    # 11, nov 94, equation (42), assuming no losses (i.e. tau -> inf)

    def imp_finite(self,t,z):
        v = self.v ; XX = N.sqrt(t**2 - (z/v)**2) ; omega_c = self.omega_c
        return -z*omega_c/v/XX*jv(1,omega_c*XX) if (t-z/v) > 0 else 0    

    def finite_response(self, z, t):
#         return (romberg(lambda tau: self.drv_fun(tau)*self.imp_finite(t-tau, z),
#                         0, self.drv_time, tol=self.epsrel, divmax=20, vec_func=False,
#                         show=True), 'blah')
#         return quadrature(lambda tau: self.drv_fun(tau)*self.imp_finite(t-tau, z),
#                           0, self.drv_time, tol=self.epsrel, maxiter=200, vec_func=False)
        return quad(lambda tau: self.drv_fun(tau)*self.imp_finite(t-tau, z),
                    0, self.drv_time, epsabs=self.epsabs, epsrel=self.epsrel,
                    limit=200)


    def full_response(self, z, t):
        v = self.v
        drv_val = self.drv_fun(t-z/v) if t < self.drv_time+z/v else 0
        return  drv_val + self.finite_response(z, t)[0]


    def __init__(self, omega_c, drv_fun, drv_time):
        """
        omega_c -- Cutoff frequency in rad/s,
        drv_fun(t) -- Driving function of time
        drv_time -- drv_fun(t) is set to 0 for t > drv_time
        """
        self.omega_c = omega_c ; self.drv_fun = drv_fun ; self.drv_time = drv_time
        

