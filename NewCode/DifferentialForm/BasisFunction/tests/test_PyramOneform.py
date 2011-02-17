from __future__ import division

from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N
from NewCode.Utilities import Struct, partial
from NewCode.tests import xfail
from NewCode.Meshes import PyramMesh
from NewCode.tests.PyramMeshes import SixPyram
from NewCode.DifferentialForm.BasisFunction import PyramOneform
import pyram_oneform_vals, pyram_oneform_interp_vals

class _test_basisfunctions(NumpyTestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.lams = self.TestMesh.test_local_coords
        self.facedual_edge_numbering = PyramMesh.Element.FACEDUAL_EDGE_NUMBERING
        self.facedual_edge_sense = PyramMesh.Element.FACEDUAL_EDGE_SENSE
        self.apexface_edges = PyramMesh.Element.LOCAL_APEXFACEEDGES

    def gen_values_edge(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(self.facedual_edge_numbering, self.facedual_edge_sense )
                 for bf_gen in bf_genlist]]

    def gen_values_apexface(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(self.apexface_edges,
                        self.facedual_edge_numbering,
                        self.facedual_edge_sense)
                 for bf_gen in bf_genlist]]

    def gen_values_basefacevol(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(self.facedual_edge_numbering, self.facedual_edge_sense)
                 for bf_gen in bf_genlist]]

class test_graglia99_basisfunctions(_test_basisfunctions):
    btype = 'graglia99'
    def test_edge0(self):
        edgefuns0 = PyramOneform.basis_set(1,btype=self.btype).fns.edge
        assert_almost_equal(self.gen_values_edge(edgefuns0), 
                            pyram_oneform_vals.edgevals_graglia99_0, decimal=15)

    def test_D_edge0(self):
        D_edgefuns0 = PyramOneform.basis_set(1,btype=self.btype).fns_D.edge
        assert_almost_equal(self.gen_values_edge(D_edgefuns0), 
                            pyram_oneform_vals.D_edgevals_graglia99_0, decimal=15)

    def test_D_edge1(self):
        D_edgefuns1 = PyramOneform.basis_set(2,btype=self.btype).fns_D.edge
        assert_almost_equal(self.gen_values_edge(D_edgefuns1), 
                            pyram_oneform_vals.D_edgevals_graglia99_1, decimal=15)

    def test_edge1(self):
        edgefuns = PyramOneform.basis_set(2,btype=self.btype).fns.edge
        assert_almost_equal(self.gen_values_edge(edgefuns), 
                            pyram_oneform_vals.edgevals_graglia99_1, decimal=15)
        
    def test_apexface1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns.apexface
        assert_almost_equal(self.gen_values_apexface(funs), 
                            pyram_oneform_vals.apexfacevals_graglia99_1, decimal=15)

    def test_D_apexface1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns_D.apexface
        assert_almost_equal(self.gen_values_apexface(funs), 
                            pyram_oneform_vals.D_apexfacevals_graglia99_1, decimal=15)

    def test_baseface1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns.baseface
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.basefacevals_graglia99_1, decimal=15)

    def test_D_baseface1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns_D.baseface
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.D_basefacevals_graglia99_1, decimal=15)

    def test_basevol1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns.vol
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.basevolvals_graglia99_1, decimal=15)

    def test_D_basevol1(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns_D.vol
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.D_basevolvals_graglia99_1, decimal=15)


class test_graglia99_interpfns(_test_basisfunctions):
    btype = 'graglia99'

    def gen_apex_values_edge(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i+4, self.facedual_edge_numbering)
                 for bf_gen in bf_genlist for i in range(4)]]

    def gen_D_apex_values_edge(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i+4, self.facedual_edge_numbering).D
                 for bf_gen in bf_genlist for i in range(4)]]

    def gen_base_values_edge(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i, self.facedual_edge_numbering)
                 for bf_gen in bf_genlist for i in range(4)]]

    def gen_D_base_values_edge(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i, self.facedual_edge_numbering).D
                 for bf_gen in bf_genlist for i in range(4)]]

    def gen_apex_values_face(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i) for bf_gen in bf_genlist for i in range(4)]]

    def gen_D_apex_values_face(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i).D for bf_gen in bf_genlist for i in range(4)]]

    def gen_base_values_face(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i) for i in (1,2) for bf_gen in bf_genlist]]

    def gen_D_base_values_face(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i).D for i in (1,2) for bf_gen in bf_genlist]]

    def gen_base_values_vol(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i) for i in (1,2) for bf_gen in bf_genlist]]

    def gen_D_base_values_vol(self, bf_genlist):
        return [[fun(coord) for coord in self.lams]
                for fun in
                [bf_gen(i).D for i in (1,2) for bf_gen in bf_genlist]]

    def test_apexedge_efs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        apexedge_ifns = ifn.get_apexedgeFuncs()
        assert_almost_equal(self.gen_apex_values_edge(apexedge_ifns), 
                            pyram_oneform_interp_vals.apexedge_efs_graglia99_p1,
                            decimal=15)

    def test_D_apexedge_efs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        apexedge_ifns = ifn.get_apexedgeFuncs()
        assert_almost_equal(self.gen_D_apex_values_edge(apexedge_ifns), 
                            pyram_oneform_interp_vals.D_apexedge_efs_graglia99_p1,
                            decimal=15)

    def test_baseedge_efs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        baseedge_ifns = ifn.get_baseedgeFuncs()
        assert_almost_equal(self.gen_base_values_edge(baseedge_ifns), 
                            pyram_oneform_interp_vals.baseedge_efs_graglia99_p1,
                            decimal=15)

    def test_D_baseedge_efs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        baseedge_ifns = ifn.get_baseedgeFuncs()
        assert_almost_equal(self.gen_D_base_values_edge(baseedge_ifns), 
                            pyram_oneform_interp_vals.D_baseedge_efs_graglia99_p1,
                            decimal=15)

    def test_apex_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        apexface_ifns = ifn.get_apexfaceFuncs()
        assert_almost_equal(self.gen_apex_values_face(apexface_ifns), 
                            pyram_oneform_interp_vals.apexedge_ffs_graglia99_p1,
                            decimal=15)

    def test_D_apex_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        apexface_ifns = ifn.get_apexfaceFuncs()
        assert_almost_equal(self.gen_D_apex_values_face(apexface_ifns), 
                            pyram_oneform_interp_vals.D_apexedge_ffs_graglia99_p1,
                            decimal=15)

    def test_base_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        baseface_ifns = ifn.get_basefaceFuncs()
        assert_almost_equal(self.gen_base_values_face(baseface_ifns), 
                            pyram_oneform_interp_vals.baseedge_ffs_graglia99_p1,
                            decimal=15)

    def test_D_base_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        baseface_ifns = ifn.get_basefaceFuncs()
        assert_almost_equal(self.gen_D_base_values_face(baseface_ifns), 
                            pyram_oneform_interp_vals.D_baseedge_ffs_graglia99_p1,
                            decimal=15)

    def test_base_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        basevol_ifns = ifn.get_basevolFuncs()
        assert_almost_equal(self.gen_base_values_vol(basevol_ifns), 
                            pyram_oneform_interp_vals.baseedge_vfs_graglia99_p1,
                            decimal=15)

    def test_D_base_ffs_p1(self):
        ifn = PyramOneform.InterpFuncs(order=1)
        basevol_ifns = ifn.get_basevolFuncs()
        assert_almost_equal(self.gen_D_base_values_vol(basevol_ifns), 
                            pyram_oneform_interp_vals.D_baseedge_vfs_graglia99_p1,
                            decimal=15)

class test_coulomb97_basisfunctions(_test_basisfunctions):
    btype = 'coulomb97'
    def test_edge0(self):
        edgefuns0 = PyramOneform.basis_set(1,btype=self.btype).fns.edge
        assert_almost_equal(self.gen_values_edge(edgefuns0), 
                            pyram_oneform_vals.edgevals_graglia99_0, decimal=15)

    def test_D_edge0(self):
        D_edgefuns0 = PyramOneform.basis_set(1,btype=self.btype).fns_D.edge
        assert_almost_equal(self.gen_values_edge(D_edgefuns0), 
                            pyram_oneform_vals.D_edgevals_graglia99_0, decimal=15)

    def test_edge1(self):
        edgefuns1 = PyramOneform.basis_set(2,btype=self.btype).fns.edge
        assert_almost_equal(self.gen_values_edge(edgefuns1), 
                            pyram_oneform_vals.edgevals_coulomb97_1, decimal=15)

    def test_D_edge1(self):
        D_edgefuns1 = PyramOneform.basis_set(2,btype=self.btype).fns_D.edge
        assert_almost_equal(self.gen_values_edge(D_edgefuns1), 
                            pyram_oneform_vals.D_edgevals_coulomb97_1, decimal=15)

    def test_base_ffs_2m(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns.baseface
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.basefacevals_coulomb97_2m, decimal=15)

    def D_test_base_ffs_2m(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns_D.baseface
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.D_basefacevals_coulomb97_2m, decimal=15)

        
    def test_basevol_2m(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns.vol
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.volfunvals_coulomb97_2m, decimal=15)

    def test_D_basevol_2m(self):
        funs = PyramOneform.basis_set(2,btype=self.btype).fns_D.vol
        assert_almost_equal(self.gen_values_basefacevol(funs), 
                            pyram_oneform_vals.D_volfunvals_coulomb97_2m, decimal=15)
