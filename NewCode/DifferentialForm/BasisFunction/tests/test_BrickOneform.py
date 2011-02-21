# Makes 1/2 (etc.) return floating point rather than 0.
from __future__ import division

from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N
from NewCode.tests import xfail
from NewCode.Utilities import Struct, partial
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm.BasisFunction import BrickOneform
import brick_oneform_vals

class _test_basisfunctions(TestCase):
    def setUp(self):
        self.facedual_edge_numbering = BrickMesh.Element.FACEDUAL_EDGE_NUMBERING
        no_steps = 3
        self.lams = N.array([(i,j,k,no_steps-i,no_steps-j,no_steps-k)
                             for i in range(no_steps+1)
                             for j in range(no_steps+1)
                             for k in range(no_steps+1)], N.float64)/no_steps

    def gen_values_edgeface(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(self.facedual_edge_numbering)
                 for bf_gen in bf_genlist]]

    def gen_values_vol(self, bf_genlist):
        return [[fn(coord) for coord in self.lams]
                for fn in bf_genlist]

class test_cohen98_basisfunctions(_test_basisfunctions):
    btype = 'cohen98'
    def test_edge_gen0(self):
        edgefuns0 = BrickOneform.basis_set(1, btype=self.btype).fns['edge']
        assert_almost_equal(self.gen_values_edgeface(edgefuns0), brick_oneform_vals.edge_vals0, 
                            decimal=15)

class test_rieben04_basisfunctions(_test_basisfunctions):
    btype = 'rieben04'
    def test_edge0(self):
        assert_almost_equal(brick_oneform_vals.edge_vals0,
                            self.gen_values_edgeface(BrickOneform.edgefuns0),
                            decimal=15)

    def test_edge_gen0(self):
        edgefuns0 = BrickOneform.basis_set(1, btype=self.btype).fns['edge']
        assert_almost_equal(self.gen_values_edgeface(edgefuns0), brick_oneform_vals.edge_vals0, 
                            decimal=15)

    def test_edge_rieben04_2m(self):
        edgefuns = BrickOneform.basis_set(2, btype=self.btype).fns['edge']
        assert_almost_equal(self.gen_values_edgeface(edgefuns),
                            brick_oneform_vals.rieben04_2m_edgefns,
                            decimal=15)

    def test_face_rieben04_2m(self):
        facefuns = BrickOneform.basis_set(2, btype=self.btype).fns['face']
        assert_almost_equal(self.gen_values_edgeface(facefuns),
                            brick_oneform_vals.rieben04_2m_facefns,
                            decimal=15)

    def test_vol_rieben04_2m(self):
        volfuns = BrickOneform.basis_set(2, btype=self.btype).fns['vol']
        assert_almost_equal(self.gen_values_vol(volfuns),
                            brick_oneform_vals.rieben04_2m_volfns,
                            decimal=15)

    def test_edge_rieben04_2m_D(self):
        edgefuns = BrickOneform.basis_set(2, btype=self.btype).fns_D['edge']
        assert_almost_equal(self.gen_values_edgeface(edgefuns),
                            brick_oneform_vals.rieben04_2m_edgefns_D,
                            decimal=15)

    def test_face_rieben04_2m_D(self):
        facefuns = BrickOneform.basis_set(2, btype=self.btype).fns_D['face']
        assert_almost_equal(self.gen_values_edgeface(facefuns),
                            brick_oneform_vals.rieben04_2m_facefns_D,
                            decimal=15)

    def test_vol_rieben04_2m_D(self):
        volfuns = BrickOneform.basis_set(2, btype=self.btype).fns_D['vol']
        assert_almost_equal(self.gen_values_vol(volfuns),
                            brick_oneform_vals.rieben04_2m_volfns_D,
                            decimal=15)

class test_basis_set(TestCase):
    def test_1_mixed(self):
        bset = BrickOneform.basis_set(1, mixed=True)
        # Check default value of mixed
        assert_equal(bset.info, Struct(form=1, type='rieben04', order=0.5))
        self.assertRaises(AssertionError, BrickOneform.basis_set, 0)
        self.assertRaises(NotImplementedError, BrickOneform.basis_set, 1000000)


