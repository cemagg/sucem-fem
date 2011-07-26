from __future__ import division

import dolfin

from sucemfem.Consts import eps0, mu0, c0, Z0
from sucemfem.PostProcessing import CalcEMFunctional
from sucemfem import Geometry 

class SurfaceFlux(object):
    """Caclulate the real part of the power-flux through a surface using a surface integral

    I.e. the dot product of the Poynting (E x complex_conjugate(H))
    vector with the surface normal is integrated to yield the power
    flux. Only the real part of this is calucated. RMS powerflow is
    given by half this quantity, since we are dealing with
    peak-values.

    Currently only the problem bounding surface can be used.

    """
    def __init__(self, function_space):
        self.function_space = V = function_space
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)

    def set_dofs(self, dofs):
        x_r = np.real(dofs).copy()
        x_i = np.imag(dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i

    def set_k0(self,k0):
        self.k0 = k0

    def _get_mur_function(self):
        if self.epsr_function is not None:
            mur_func = self.mur_function
        else:
            mur_func = dolfin.Constant(1)

        return mur_func

    def set_mur_function(self, mur_function):
        self.mur_function = mur_function

    def _get_form(self):
        n = self.function_space.cell().x
        k0 = self.k0
        E_i = self.E_i
        E_r = self.E_r
        mu_r = self._get_mur_function()
        return (1/k0/Z0)*dolfin.dot(n, (dolfin.cross(E_r, -dolfin.curl(E_i)/mu_r) +
                                        dolfin.cross(E_i, dolfin.curl(E_r)/mu_r)))*dolfin.ds

    def calc_flux(self):
        """Calculate the power flux"""
        return dolfin.assemble(self._get_form)

class VariationalSurfaceFlux(object):
    def __init__(self, function_space):
        self.function_space = V = function_space
        self.E_r = dolfin.Function(V)
        self.E_i = dolfin.Function(V)
        self.mur_function = None
        ## Set CalcEMCalcEMFunctional to only integrate along a skin
        ## of elements connected to the boundary
        boundary_cells = Geometry.BoundaryEdgeCells(V.mesh())
        cell_domains = dolfin.CellFunction('uint', V.mesh())
        cell_domains.set_all(0)
        cell_region = 1
        boundary_cells.mark(cell_domains, cell_region)
        self.functional = CalcEMFunctional(V)
        self.functional.set_cell_domains(cell_domains, cell_region)

    def set_dofs(self, dofs):
        x_r = np.real(dofs).copy()
        x_i = np.imag(dofs).copy()
        self.E_r.vector()[:] = x_r
        self.E_i.vector()[:] = x_i
        self.dofs = dofs
        
        boundary = dolfin.DomainBoundary()
        E_r_dirich = dolfin.DirichletBC(
            self.function_space, self.E_r, boundary)
        E_i_dirich = dolfin.DirichletBC(
            self.function_space, self.E_i, boundary)
        x_r_dirich = np.zeros(len(x_r))
        x_i_dirich = np.zeros(len(x_r))
        E_r_dirich.apply(x_r_dirich)
        E_i_dirich.apply(x_i_dirich)
        self.dirich_dofs = x_r_dirich.array() + 1j*x_i_dirich.array()
        self.functional.set_E_dofs(dofs)
        
    def set_k0(self,k0):
        self.k0 = k0
        self.functional.set_k0(k0)
        self.functional.set_g_dofs(1j*self.dirich_dofs.conjugate()/k0/Z0)
        
    def _get_mur_function(self):
        if self.epsr_function is not None:
            mur_func = self.mur_function
        else:
            mur_func = dolfin.Constant(1)

        return mur_func

    def set_mur_function(self, mur_function):
        self.mur_function = mur_function
        self.functional.set_mur_function(mur_function)
