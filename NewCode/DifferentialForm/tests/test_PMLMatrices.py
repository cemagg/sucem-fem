from __future__ import division

import pickle, os
import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.Meshes import BrickMesh, BrickMeshGen
from NewCode.DifferentialForm import PMLMatrices, BrickDiscretiser, allfree

class test_PMLMatrices(NumpyTestCase):
    h = 1.
    a,b,c = 1.,2.,3.
    no_PML_cells = 2
    PML_m = 3.
    order = 2 
    def setUp(self):
        ga,gb,gc = self.a, self.b, self.c
        h = self.h
        PML_pos = (ga/2, gb/2, gc/2)
        PML_m = self.PML_m
        PML_sigma_max = 0.8*(PML_m+1)/h
        PML_len = h*self.no_PML_cells

        a,b,c = N.array([ga,gb,gc]) + PML_len*2
        self.mesh = BrickMesh.Mesh(BrickMeshGen.make_rect_cavity_brick_listmesh(
            a,b,c, h, grid_offset=[-a/2,-b/2,-c/2]))

        sigma_x_fn = lambda r: N.where(
            N.abs(r.T[0]) >= PML_pos[0],
            PML_sigma_max*((N.abs(r.T[0])-PML_pos[0])/PML_len)**PML_m,
            0.)
        sigma_y_fn = lambda r: N.where(
            N.abs(r.T[1]) >= PML_pos[1],
            PML_sigma_max*((N.abs(r.T[1])-PML_pos[1])/PML_len)**PML_m,
            0.)
        sigma_z_fn = lambda r: N.where(
            N.abs(r.T[2]) >= PML_pos[2],
            PML_sigma_max*((N.abs(r.T[2])-PML_pos[2])/PML_len)**PML_m,
            0.)

        self.sigma_fns = Struct(x=sigma_x_fn, y=sigma_y_fn, z=sigma_z_fn)

        order = self.order
        self.discs = Struct(E=BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 1, order, mixed=True, freeFun=allfree, btype='cohen98'),
                            B=BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, 2, order, mixed=True, freeFun=allfree, btype='cohen98'))
        self.discs.E.diagonalise()
        self.discs.B.diagonalise()
        self.inst = PMLMatrices.PMLMatrices(self.discs, self.sigma_fns)

    def test_sigma_xyz(self):
        sigma_x = self.inst.sigma_x()
        sigma_y = self.inst.sigma_y()
        sigma_z = self.inst.sigma_z()
        assert_almost_equal(sigma_x, sigma_mats.sigma_x, decimal=15)
        assert_almost_equal(sigma_y, sigma_mats.sigma_y, decimal=15)
        assert_almost_equal(sigma_z, sigma_mats.sigma_z, decimal=15)
        
    def test_inv_sigma_xyz_inv(self):
        inv_sigma_x_inv = self.inst.inv_sigma_x_inv()
        inv_sigma_y_inv = self.inst.inv_sigma_y_inv()
        inv_sigma_z_inv = self.inst.inv_sigma_z_inv()
        assert_almost_equal(inv_sigma_x_inv, sigma_mats.inv_sigma_x_inv, decimal=15)
        assert_almost_equal(inv_sigma_y_inv, sigma_mats.inv_sigma_y_inv, decimal=15)
        assert_almost_equal(inv_sigma_z_inv, sigma_mats.inv_sigma_z_inv, decimal=15)

# Not neccesarily correct, as generated using
# nmarais@sun.ac.za--femcode/newcode--main--0.2--patch-142
sigma_mats = pickle.load(file(os.path.join(
    os.path.dirname(__file__), 'sigma_mats.pickle')))
