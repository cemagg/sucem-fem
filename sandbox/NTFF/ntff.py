from __future__ import division

import numpy as N
import dolfin
from dolfin import curl, cross, dx, ds, Constant, dot
from FenicsCode.Consts import Z0, c0

class NTFF(object):
    def __init__(self, function_space):
        self.function_space = V = function_space
        # Outward pointing normal
        self.n = V.cell().n
        # \vec{r'}, i.e. rprime is simply the position vector
        self.rprime = V.cell().x
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)
        self._L = []
        self._N = []

    def set_dofs(self, dofs):
        x_r = N.real(dofs).copy()
        x_i = N.imag(dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i
        #------------------------------
        # Set up eqivalent magnetic current forms
        self.M_r = -cross(self.n, self.E_r)
        self.M_i = -cross(self.n, self.E_i)

    def set_frequency(self, frequency):
        self.frequency = frequency
        self.init_calc()

    def init_calc(self):
        self.k0 = k0 = self.frequency*2*N.pi/c0
        #------------------------------
        # Set up magnetic field and equivalent electric current forms
        self.H_r = -curl(self.E_i)/(self.k0*Z0)
        self.H_i = curl(self.E_r)/(self.k0*Z0)
        self.J_r = cross(self.n, self.H_r)
        self.J_i = cross(self.n, self.H_i)
        
        
    def calc_pt(self, theta_deg, phi_deg):
        theta = N.deg2rad(theta_deg)
        phi = N.deg2rad(phi_deg)
        rhat_ = N.array([N.sin(theta)*N.cos(phi),
                         N.sin(theta)*N.sin(phi),
                         N.cos(theta)], dtype=N.float64)
        theta_hat_ = N.array([N.cos(theta)*N.cos(phi),
                              N.cos(theta)*N.sin(phi),
                              -N.sin(theta)], dtype=N.float64)
        phi_hat_ = N.array([-N.sin(phi), N.cos(phi), 0.], dtype=N.float64)
        rhat = Constant(list(rhat_))
        theta_hat = Constant(list(theta_hat_))
        phi_hat = Constant(list(phi_hat_))
        M_r = self.M_r
        M_i = self.M_i
        J_r = self.J_r
        J_i = self.J_i
        k0 = self.k0 
        # phase term to be used in sin/cos
        phase = k0*dot(self.rprime, rhat)
        #------------------------------
        # Set up form for far field potential N
        N_r = J_r*dolfin.cos(phase) - J_i*dolfin.sin(phase)
        N_i = J_r*dolfin.sin(phase) + J_i*dolfin.cos(phase)
        # UFL does not seem to like vector valued functionals, so we split the
        # final functional from into theta and phi components
        N_r_theta = dot(theta_hat, N_r)*ds
        N_r_phi = dot(phi_hat, N_r)*ds
        N_i_theta = dot(theta_hat, N_i)*ds
        N_i_phi = dot(phi_hat, N_i)*ds

        #------------------------------
        # Set up form for far field potential L
        L_r = M_r*dolfin.cos(phase) - M_i*dolfin.sin(phase)
        L_i = M_r*dolfin.sin(phase) + M_i*dolfin.cos(phase)
        L_r_theta = dot(theta_hat, L_r)*ds
        L_r_phi = dot(phi_hat, L_r)*ds
        L_i_theta = dot(theta_hat, L_i)*ds
        L_i_phi = dot(phi_hat, L_i)*ds

        # from ufl.algorithms import estimate_total_polynomial_degree
        # print "L degree: ", estimate_total_polynomial_degree(L_r_theta)

        #------------------------------
        # Now evaluate numerically, adding the real and imaginary parts
        # ffc_opt = {"quadrature_degree": 2, "representation": "quadrature"}
        # L_theta = dolfin.assemble(L_r_theta, form_compiler_parameters=ffc_opt) \
        #           + 1j*dolfin.assemble(L_i_theta, form_compiler_parameters=ffc_opt)
        # L_phi = dolfin.assemble(L_r_phi, form_compiler_parameters=ffc_opt) \
        #         + 1j*dolfin.assemble(L_i_phi, form_compiler_parameters=ffc_opt)
        # N_theta = dolfin.assemble(N_r_theta, form_compiler_parameters=ffc_opt) \
        #           + 1j*dolfin.assemble(N_i_theta, form_compiler_parameters=ffc_opt)
        # N_phi = dolfin.assemble(N_r_phi, form_compiler_parameters=ffc_opt) \
        #         + 1j*dolfin.assemble(N_i_phi, form_compiler_parameters=ffc_opt)    

        L_theta = dolfin.assemble(L_r_theta) + 1j*dolfin.assemble(L_i_theta)
        L_phi = dolfin.assemble(L_r_phi) + 1j*dolfin.assemble(L_i_phi)
        N_theta = dolfin.assemble(N_r_theta) + 1j*dolfin.assemble(N_i_theta)
        N_phi = dolfin.assemble(N_r_phi) + 1j*dolfin.assemble(N_i_phi)    

        self._L.append([L_theta, L_phi])
        self._N.append([N_theta, N_phi])

        #------------------------------
        # Calculate the far fields normalised to radius 1.
        r_fac = 1j*k0*N.exp(-1j*k0)/(4*N.pi)
        E_theta = -r_fac*(L_phi + Z0*N_theta)
        E_phi = r_fac*(L_theta - Z0*N_phi)
        return (E_theta, E_phi)
