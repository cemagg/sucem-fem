# Makes 1/2 (etc.) return floating point rather than 0.
from __future__ import division

from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N
from NewCode.tests import xfail
from NewCode.Utilities import Struct
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm.BasisFunction import BrickTwoform
import brick_twoform_vals

class test_basisfunctions(TestCase):
    def setUp(self):
        self.facedual_edge_numbering = BrickMesh.Element.FACEDUAL_EDGE_NUMBERING
        no_steps = 3
        self.lams = N.array([(i,j,k,no_steps-i,no_steps-j,no_steps-k)
                             for i in range(no_steps+1)
                             for j in range(no_steps+1)
                             for k in range(no_steps+1)], N.float64)/no_steps

    def gen_values_face(self, bf_genlist):
        return [[facefun(None)(coord) for coord in self.lams]
                for facefun in bf_genlist]

    def gen_values_vol(self, bf_genlist):
        return [[fn(coord) for coord in self.lams]
                for fn in bf_genlist]

    def test_face0(self):
        assert_almost_equal(brick_twoform_vals.face_vals0,
                            self.gen_values_face(BrickTwoform.facefuns0),
                            decimal=15)

    def test_face0_D(self):
        # lam is actually not used by the lowest order face_D funcs so just
        # pass 0 as a dummy value
        assert_almost_equal(
            [facefun_D(None)(0) for facefun_D in BrickTwoform.facefuns0_D],
            N.array([[-1], [-1], [-1], [1], [1], [1]]), decimal=15)

    def test_face_rieben04_2m(self):
        facefuns = BrickTwoform.basis_set(2).fns['face']
        assert_almost_equal(brick_twoform_vals.rieben04_2m_facefns,
                            self.gen_values_face(facefuns),
                            decimal=15)

    def test_face_rieben04_2m_D(self):
        facefuns = BrickTwoform.basis_set(2).fns_D['face']
        assert_almost_equal(brick_twoform_vals.rieben04_2m_facefns_D,
                            self.gen_values_face(facefuns),
                            decimal=15)

    def test_vol_rieben04_2m(self):
        volfuns = BrickTwoform.basis_set(2).fns['vol']
        assert_almost_equal(brick_twoform_vals.rieben04_2m_volfns,
                            self.gen_values_vol(volfuns),
                            decimal=15)

    def test_vol_rieben04_2m_D(self):
        volfuns = BrickTwoform.basis_set(2).fns_D['vol']
        assert_almost_equal(brick_twoform_vals.rieben04_2m_volfns_D,
                            self.gen_values_vol(volfuns),
                            decimal=15)


class test_basis_set(TestCase):
    def test_1_mixed(self):
        bset = BrickTwoform.basis_set(1)
        desired_info = Struct(form=2, type='rieben04', order=0.5)
        assert_equal(bset.info, desired_info)
        self.assertRaises(AssertionError, BrickTwoform.basis_set, 0)
        self.assertRaises(NotImplementedError, BrickTwoform.basis_set, 1000000)
       


