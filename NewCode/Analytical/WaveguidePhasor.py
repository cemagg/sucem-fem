from __future__ import division

import numpy as N

class TE01(object):
    def __init__(self, a, c=1.):
        self.a = a
        self.c = c
        self.k_c = N.pi/a
        self.omega_c = self.k_c/c

    def gamma(self, omega):
        return 1j*self.beta(omega) if omega > self.omega_c \
               else self.alpha(omega)

    def alpha(self, omega):
        return N.pi/self.a*N.sqrt(1.-(omega/self.omega_c)**2)

    def beta(self, omega):
        return omega/self.c*N.sqrt(1-(self.omega_c/omega)**2)

    def TF(self, omega, z):
        gam = self.gamma(omega)
        if z >= 0:
            return N.exp(-gam*z)
        else:
            return N.exp((-1j*N.imag(gam) + N.real(gam))*z)
