from __future__ import division

import numpy as np
import dolfin
from FenicsCode.Consts import Z0, c0
from FenicsCode import Geometry 
from FenicsCode.Interpolation import SurfaceInterpolant
from FenicsCode.PostProcessing import CalcEMFunctional
from FenicsCode.PostProcessing.surface_ntff import SurfaceNTFFForms

class TransformTestingExpression(object):
    """Dolfin expression for

    1j/k0 rhat x (ahat x rhat)*e^(rhat . rprime)
    """
    def __init__(self):
        self.expr_r = dolfin.Expression([
            '-cf_x*sin(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)',
            '-cf_y*sin(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)',
            '-cf_z*sin(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)'
            ],
                                        cf_x=0, cf_y=0, cf_z=0, k0=0,
                                        rhat_z=0, rhat_y=0, rhat_x=0)
        self.expr_i = dolfin.Expression([
            'cf_x*cos(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)',
            'cf_y*cos(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)',
            'cf_z*cos(x[2]*k0*rhat_z+x[1]*k0*rhat_y+x[0]*k0*rhat_x)',
            ],
                                        cf_x=0, cf_y=0, cf_z=0, k0=0,
                                        rhat_z=0, rhat_y=0, rhat_x=0)
        
    def set_parms(self, rhat, ahat, k0):
        crossy = np.cross(rhat, np.cross(ahat, rhat))
        constfac = crossy/k0
        self._set_expr_parms(self.expr_r, rhat, k0, constfac)
        self._set_expr_parms(self.expr_i, rhat, k0, constfac)

    def _set_expr_parms(self, expr, rhat, k0, constfac):
        expr.cf_x, expr.cf_y, expr.cf_z = constfac
        expr.k0 = k0
        expr.rhat_x, expr.rhat_y, expr.rhat_z = rhat

    def get_expression(self):
        return self.expr_r, self.expr_i
    

class NTFF(object):
    def __init__(self, function_space, testing_space=None):
        self.function_space = function_space
        if testing_space is None: 
            testing_space = function_space 
        self.testing_space = testing_space 
        self.functional = CalcEMFunctional(function_space, testing_space)
        self.testing_interpolator = SurfaceInterpolant(testing_space)
        self.cell_domains = dolfin.CellFunction('uint', self.function_space.mesh())
        self.cell_domains.set_all(0)
        self.cell_region = 1
        boundary_cells = Geometry.BoundaryEdgeCells(self.function_space.mesh())
        boundary_cells.mark(self.cell_domains, self.cell_region)
        self.surface_forms = SurfaceNTFFForms(self.function_space)
        self.testing_expression_gen = TransformTestingExpression()
        self.functional.set_cell_domains(self.cell_domains, self.cell_region)

    def set_k0(self, k0):
        self.k0 = k0
        self.functional.set_k0(k0)

    def set_frequency(self, frequency):
        self.frequency = frequency
        self.set_k0(self.frequency*2*np.pi/c0)

    def set_dofs(self, dofs):
        self.functional.set_E_dofs(dofs)
        self.surface_forms.set_dofs(dofs)

    def calc_pt(self, theta_deg, phi_deg):
        # H-field contribution using variational calculation
        E_H_theta, E_H_phi = self.calc_pt_E_H(theta_deg, phi_deg)
        theta = np.deg2rad(theta_deg)
        phi = np.deg2rad(phi_deg)
        self.surface_forms.set_parms(theta, phi, self.k0)
        L_theta, L_phi = self.surface_forms.assemble_L()
        #------------------------------
        # Calculate the far fields normalised to radius 1.
        r_fac = 1j*self.k0*np.exp(-1j*self.k0)/(4*np.pi)
        E_theta = r_fac*(-L_phi + E_H_theta)
        E_phi = r_fac*(L_theta + E_H_phi)
        return (E_theta, E_phi)
    
    def calc_pt_E_H(self, theta_deg, phi_deg):
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
        self.testing_expression_gen.set_parms(rhat, ahat, self.k0)
        fit_expr_r, fit_expr_i = self.testing_expression_gen.get_expression()
        self.testing_interpolator.set_interpolant_expression(fit_expr_r, fit_expr_i)
        g_dofs = self.testing_interpolator.calculate_interpolation()
        self.functional.set_g_dofs(g_dofs)
        return self.functional.calc_functional()

