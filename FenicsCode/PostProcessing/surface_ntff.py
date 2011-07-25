## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as np
import dolfin
from dolfin import curl, cross, dx, ds, Constant, dot
from FenicsCode.Consts import Z0, c0
import FenicsCode.PostProcessing.ntff_expressions as ntff_expressions

class SurfaceNTFFForms(object):
    def __init__(self, function_space):
        self.function_space = V = function_space
        self.n = V.cell().n
        # \vec{r'}, i.e. rprime is simply the position vector
        self.rprime = V.cell().x
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)
        self.r_hat = ntff_expressions.get_r_hat()
        self.k0 = ntff_expressions.get_k0()
        self.theta_hat = ntff_expressions.get_theta_hat()
        self.phi_hat = ntff_expressions.get_phi_hat()
        # phase term to be used in sin/cos
        self.phase = ntff_expressions.get_phase(self.k0, self.rprime, self.r_hat)

    def set_dofs(self, dofs):
        x_r = np.real(dofs).copy()
        x_i = np.imag(dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i

    def get_N_form(self):
        try:
            return self.N_form
        except AttributeError:
            pass
        # Set up magnetic field and equivalent electric current forms
        H_r = -curl(self.E_i)/(self.k0*Z0)
        H_i = curl(self.E_r)/(self.k0*Z0)
        J_r = cross(self.n, H_r)
        J_i = cross(self.n, H_i)
        #------------------------------
        # Set up form for far field potential N
        theta_hat = self.theta_hat
        phi_hat = self.phi_hat
        phase = self.phase
        N_r = J_r*dolfin.cos(phase) - J_i*dolfin.sin(phase)
        N_i = J_r*dolfin.sin(phase) + J_i*dolfin.cos(phase)
        # UFL does not seem to like vector valued functionals, so we split the
        # final functional from into theta and phi components
        self.N_form = dict(
            r_theta=dot(theta_hat, N_r)*ds,
            r_phi=dot(phi_hat, N_r)*ds,
            i_theta=dot(theta_hat, N_i)*ds,
            i_phi=dot(phi_hat, N_i)*ds)
        return self.N_form

    def get_L_form(self):
        try:
            return self.L_form
        except AttributeError:
            pass
        # Set up eqivalent magnetic current forms
        M_r = -cross(self.n, self.E_r)
        M_i = -cross(self.n, self.E_i)
        #------------------------------
        # Set up form for far field potential L
        theta_hat = self.theta_hat
        phi_hat = self.phi_hat
        phase = self.phase
        L_r = M_r*dolfin.cos(phase) - M_i*dolfin.sin(phase)
        L_i = M_r*dolfin.sin(phase) + M_i*dolfin.cos(phase)
        self.L_form = dict(
            r_theta=dot(theta_hat, L_r)*ds,
            r_phi=dot(phi_hat, L_r)*ds,
            i_theta=dot(theta_hat, L_i)*ds,
            i_phi=dot(phi_hat, L_i)*ds)
        return self.L_form


    def assemble_N(self):
        #------------------------------
        # evaluate numerically, adding the real and imaginary parts
        N = self.get_N_form()
        N_theta = dolfin.assemble(N['r_theta']) + 1j*dolfin.assemble(N['i_theta'])
        N_phi = dolfin.assemble(N['r_phi']) + 1j*dolfin.assemble(N['i_phi'])
        return (N_theta, N_phi)

    def assemble_L(self):
        #------------------------------
        # evaluate numerically, adding the real and imaginary parts
        L = self.get_L_form()
        L_theta = dolfin.assemble(L['r_theta']) + 1j*dolfin.assemble(L['i_theta'])
        L_phi = dolfin.assemble(L['r_phi']) + 1j*dolfin.assemble(L['i_phi'])
        return (L_theta, L_phi)

    def set_parms(self, theta, phi, k0):
        self.k0.k0 = k0
        self.theta_hat.theta = theta
        self.theta_hat.phi = phi
        self.phi_hat.theta = theta
        self.phi_hat.phi = phi
        self.r_hat.theta = theta
        self.r_hat.phi = phi


class NTFF(object):
    def __init__(self, function_space):
        self.function_space = V = function_space
        # Outward pointing normal
        self.forms = SurfaceNTFFForms(V)
        self._L = []
        self._N = []
        self._r_fac = []

    def set_dofs(self, dofs):
        self.forms.set_dofs(dofs)

    def set_frequency(self, frequency):
        self.frequency = frequency
        self.k0 = self.frequency*2*np.pi/c0
            
    def calc_pt(self, theta_deg, phi_deg):
        theta = np.deg2rad(theta_deg)
        phi = np.deg2rad(phi_deg)
        self.forms.set_parms(theta, phi, self.k0)
        L_theta, L_phi = self.forms.assemble_L()
        N_theta, N_phi = self.forms.assemble_N()
        #------------------------------
        # Calculate the far fields normalised to radius 1.
        r_fac = 1j*self.k0*np.exp(-1j*self.k0)/(4*np.pi)
        E_theta = -r_fac*(L_phi + Z0*N_theta)
        E_phi = r_fac*(L_theta - Z0*N_phi)
        self._L.append([L_theta, L_phi])
        self._N.append([N_theta, N_phi])
        self._r_fac.append(r_fac)
        return (E_theta, E_phi)
