from __future__ import division

from itertools import izip
import numpy as N
from NewCode.Utilities import CacheLast, Struct
from NewCode.DifferentialForm import DiscretiserMatrices
from NewCode.MatrixUtils import extract_diag_mat

spaceDims = 3

class PMLMatrices(DiscretiserMatrices.Matrices):
    no_dt_dep = set(['sigma_x', 'sigma_y', 'sigma_z',
                     'inv_sigma_x_inv', 'inv_sigma_y_inv', 'inv_sigma_z_inv'])
    eps0 = 1                       # Fake!
    
    def __init__(self, discs, sigma_fns):
        """
        sigma_fns is a Struct such that sigma.x returns the material property
        sigma_x, etc.
        """
        self.discs = discs
        self.sigma_tensor_fns = Struct(
            x=[sigma_fns.x, sigma_fns.y, sigma_fns.z],
            y=[sigma_fns.y, sigma_fns.z, sigma_fns.x],
            z=[sigma_fns.z, sigma_fns.x, sigma_fns.y])
        
    def set_dt(self, dt):
        for name, ca in self.cachedAttrs.iteritems():
            if name not in self.no_dt_dep: ca.clearCache()
        self.dt = dt

    def calc_material_values(self, disc, material_fn):
        print "Calculating material array"
        m = N.zeros(disc.totalDOFs, N.float64)
        int_points, int_coords = disc.interpolationPointsCoords()
        for i in range(spaceDims):
            ind = int_coords==i
            m[ind] = material_fn[i](int_points[ind])
        return m

    @CacheLast.CachedMethod
    def sigma_x(self):
        return self.calc_material_values(self.discs.E, self.sigma_tensor_fns.x)

    @CacheLast.CachedMethod
    def sigma_y(self):
        return self.calc_material_values(self.discs.E, self.sigma_tensor_fns.y)

    @CacheLast.CachedMethod
    def sigma_z(self):
        return self.calc_material_values(self.discs.E, self.sigma_tensor_fns.z)

    @CacheLast.CachedMethod
    def inv_sigma_x_inv(self):
        return self.calc_material_values(self.discs.B, self.sigma_tensor_fns.x)

    @CacheLast.CachedMethod
    def inv_sigma_y_inv(self):
        return self.calc_material_values(self.discs.B, self.sigma_tensor_fns.y)

    @CacheLast.CachedMethod
    def inv_sigma_z_inv(self):
        return self.calc_material_values(self.discs.B, self.sigma_tensor_fns.z)

    @CacheLast.CachedMethod
    def A_dy(self):
        L1 = self.discs.E.matrix.mass()
        return L1*(1/self.dt + 1/2/self.eps0*self.sigma_y())

    @CacheLast.CachedMethod
    def B_dy(self):
        L1 = self.discs.E.matrix.mass()
        return L1*(1/self.dt - 1/2/self.eps0*self.sigma_y())

    @CacheLast.CachedMethod
    def A_dx(self):
        return 1/self.dt + 1/2/self.eps0*self.sigma_x()

    @CacheLast.CachedMethod
    def B_dx(self):
        return 1/self.dt - 1/2/self.eps0*self.sigma_x()

    @CacheLast.CachedMethod
    def A_ez(self):
        return 1/self.dt + 1/2/self.eps0*self.sigma_z()

    @CacheLast.CachedMethod
    def B_ez(self):
        return 1/self.dt - 1/2/self.eps0*self.sigma_z()

    @CacheLast.CachedMethod
    def A_by(self):
        return 1/self.dt + 1/2/self.eps0*self.inv_sigma_y_inv()

    @CacheLast.CachedMethod
    def B_by(self):
        return 1/self.dt - 1/2/self.eps0*self.inv_sigma_y_inv()

    @CacheLast.CachedMethod
    def A_hz(self):
        return 1/self.dt + 1/2/self.eps0*self.inv_sigma_z_inv()

    @CacheLast.CachedMethod
    def B_hz(self):
        return 1/self.dt - 1/2/self.eps0*self.inv_sigma_z_inv()

    @CacheLast.CachedMethod
    def A_bx(self):
        return 1/self.dt + 1/2/self.eps0*self.inv_sigma_x_inv()

    @CacheLast.CachedMethod
    def B_bx(self):
        return 1/self.dt - 1/2/self.eps0*self.inv_sigma_x_inv()

