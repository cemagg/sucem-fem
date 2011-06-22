from __future__ import division

import numpy as np
import dolfin
from dolfin import curl, cross, dx, ds, Constant, dot
from FenicsCode.Consts import Z0, c0
from surface_ntff import SurfaceNTFFForms

class SurfaceInterpolant(object):
    """Calculate a complex surface interpolant on a 3D domain"""
    def __init__(self, function_space):
        V = self.function_space = function_space
        self.u_r = dolfin.Function(V)
        self.u_i = dolfin.Function(V)
        
    def set_interpolant(self, interpolant):
        """Set interpolant f([x,y,z]) -> [f_x, f_y, f_z]"""
        self.interpolant = interpolant

        class int_Re(dolfin.Expression):
            def __init__(self, interpolant):
                self.interpolant = interpolant
                
            def value_shape(self):
                return (3,)

            def eval(self, value, x):
                value[:] = self.interpolant(x).real

        class int_Im(dolfin.Expression):
            def __init__(self, interpolant):
                self.interpolant = interpolant

            def value_shape(self):
                return (3,)

            def eval(self, value, x):
                value[:] = self.interpolant(x).imag

        self.interpolant_expression_Re = int_Re(interpolant)
        self.interpolant_expression_Im = int_Im(interpolant)
        
    def calculate_interpolation(self):
        def boundary(x, on_boundary):
            return on_boundary
        V = self.function_space
        bc_Re = dolfin.DirichletBC(
            V, self.interpolant_expression_Re, boundary)
        bc_Im = dolfin.DirichletBC(
            V, self.interpolant_expression_Im, boundary)
        bc_Re.apply(self.u_r.vector())
        bc_Im.apply(self.u_i.vector())
        x = np.zeros(self.u_r.vector().size(), np.complex128)
        x[:] = self.u_r.vector().array() + 1j*self.u_i.vector().array()
        return x

class TransformTestingFunction(object):
    def __init__(self, rhat, ahat, k0):
        self.rhat, self.ahat, self.k0 = rhat, ahat, k0
        self.constfac = 1j/k0*np.cross(rhat, np.cross(ahat, rhat))
        
    def __call__(self, rprime):
        return self.constfac*np.exp(1j*self.k0*np.dot(self.rhat, rprime))

class CalcEMFunctional(object):
    """Evaluate EM functional, assuming freespace

    Evaluate <curl(E), curl(g)> - k0^2 <E, g>
    where E is the E-field solution and g is an H(curl) testing function
    
    """

    def __init__(self, function_space, testing_space=None):
        V = self.function_space = function_space
        if testing_space is None:
            self.testing_space = V
        else:
            self.testing_space = testing_space
        Vt = self.testing_space
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)
        self.g_r = dolfin.Function(Vt)
        self.g_i = dolfin.Function(Vt)

    def set_k0(self, k0):
        self.k0 = k0
        self.form_r, self.form_i = self._get_forms()

    def _get_forms(self):
        E_r, E_i, g_r, g_i = self.E_r, self.E_i, self.g_r, self.g_i
        k0 = self.k0
        form_r = (dot(curl(E_r), curl(g_r)) - dot(curl(E_i), curl(g_i)) \
                 + k0**2*dot(E_r, g_r) - dot(E_i, g_i))*dx
        form_i = (dot(curl(E_r), curl(g_i)) + dot(curl(E_i), curl(g_r)) \
                 + k0**2*dot(E_r, g_i) - dot(E_i, g_r))*dx
        return form_r, form_i

    def set_E_dofs(self, E_dofs):
        x_r = np.real(E_dofs).copy()
        x_i = np.imag(E_dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i
        
    def set_g_dofs(self, g_dofs):
        x_r = np.real(g_dofs).copy()
        x_i = np.imag(g_dofs).copy()
        self.g_r.vector()[:] = x_r
        self.g_i.vector()[:] = x_i

    def calc_functional(self):
        """Calculate functional using given E, g and k0 values"""
        I_r = dolfin.assemble(self.form_r)
        I_i = dolfin.assemble(self.form_i)
        return I_r + 1j*I_i

class NTFF(object):
    def __init__(self, function_space, testing_space=None):
        self.function_space = function_space
        if testing_space is None: 
            testing_space = function_space 
        self.testing_space = testing_space 
        self.functional = CalcEMFunctional(function_space, testing_space)
        self.testing_interpolator = SurfaceInterpolant(testing_space)
        self.forms = SurfaceNTFFForms(self.function_space)

    def set_k0(self, k0):
        self.k0 = k0
        self.functional.set_k0(k0)

    def set_frequency(self, frequency):
        self.frequency = frequency
        self.set_k0(self.frequency*2*np.pi/c0)

    def set_dofs(self, dofs):
        self.functional.set_E_dofs(dofs)
        self.forms.set_dofs(dofs)

    def calc_pt(self, theta_deg, phi_deg):
        # H-field contribution using variational calculation
        E_H_theta, E_H_phi = self.calc_pt_E_H(theta_deg, phi_deg)
        theta = np.deg2rad(theta_deg)
        phi = np.deg2rad(phi_deg)
        self.forms.set_parms(theta, phi, self.k0)
        L_theta, L_phi = self.forms.assemble_L()
        #------------------------------
        # Calculate the far fields normalised to radius 1.
        r_fac = 1j*self.k0*np.exp(-1j*self.k0)/(4*np.pi)
        E_theta = r_fac*(-L_phi + E_H_theta)
        E_phi = r_fac*(L_theta + E_H_phi)
        return (E_theta, E_phi)
    
    def calc_pt_E_H(self, theta_deg, phi_deg):
        print theta_deg
        theta = np.deg2rad(theta_deg)
        phi = np.deg2rad(phi_deg)
        rhat = np.array([np.sin(theta)*np.cos(phi),
                        np.sin(theta)*np.sin(phi),
                        np.cos(theta)], dtype=np.float64)
        theta_hat = np.array([np.cos(theta)*np.cos(phi),
                             np.cos(theta)*np.sin(phi),
                             -np.sin(theta)], dtype=np.float64)
        phi_hat = np.array([-np.sin(phi), np.cos(phi), 0.], dtype=np.float64)
        E_H_theta = self.calc_ff_func(rhat, theta_hat)
        E_H_phi = self.calc_ff_func(rhat, phi_hat)
        return E_H_theta, E_H_phi
    
    def calc_ff_func(self, rhat, ahat):
        fit_fn = TransformTestingFunction(rhat, ahat, self.k0)
        self.testing_interpolator.set_interpolant(fit_fn)
        g_dofs = self.testing_interpolator.calculate_interpolation()
        self.functional.set_g_dofs(g_dofs)
        return self.functional.calc_functional()
