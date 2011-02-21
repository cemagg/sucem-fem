from __future__ import division

from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from scipy import sparse

from NewCode.Utilities import Struct
from NewCode.DifferentialForm import BrickDiscretiser
from NewCode.Integration import BrickTrapzIntegrator
from NewCode.Meshes import BrickMesh
import NewCode.ImplicitExplicit as IEH
import NewCode.ImplicitExplicit.NewmarkHybridMatrices as NHM

from test_ImplicitExplicit import _FourBrickLinearHybridBlockMatrix
import HybridMatrices

class test_HybridBlockMats(_FourBrickLinearHybridBlockMatrix):
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)
        self.elgroups = IEH.gen_elgroups(self.mesh, self.implicitElements)
        self.group_freefuns = IEH.gen_group_freefuns(
            self.mesh, self.elgroups, self.global_freefun)
        self.discs = Struct(
            E=Struct((k,self.setup_disc(k)) for k in self.group_freefuns))
        self.inst = NHM.HybridBlockMats(
            self.discs, self.implicit_beta, self.implicitElements)
        self.inst.set_dt(self.dt)
        
    def setup_disc(self, name):
        ff, ff_vol = self.group_freefuns[name]
        disc = BrickDiscretiser.setup_PformDiscretiser(
            self.mesh, form=1,order=1,mixed=True,freeFun=ff, vol_freeFun=ff_vol)
        disc.set_integrator(BrickTrapzIntegrator)
        disc.setIntegrationRule(0)
        disc.D().set_integrator(BrickTrapzIntegrator)
        disc.D().setIntegrationRule(0)
        return disc

    def test_blockmats_A(self):
        self.verify_block_mat('A', self.block_mats['A'], self.inst, HybridMatrices)

    def test_blockmats_B(self):
        self.verify_block_mat('B', self.block_mats['B'], self.inst, HybridMatrices)
        

class test_HybridMergedMats(TestCase):
    def setUp(self):
        class MockHybridBlockMats(object):
            discs = HybridMatrices.discs
            def __getattr__(self, attr):
                return lambda : sparse.coo_matrix(getattr(HybridMatrices, attr))
        
        self.inst = NHM.HybridMergedMats(MockHybridBlockMats())

    def test_A_imp(self):
        assert_almost_equal(self.inst.A_imp().todense(), HybridMatrices.A_imp,
                            decimal=16)

    def test_A_exp(self):
        assert_almost_equal(self.inst.A_exp().todense(), HybridMatrices.A_exp,
                            decimal=16)
    def test_A(self):
        assert_almost_equal(self.inst.A().todense(), HybridMatrices.A,
                            decimal=16)

    def test_B(self):
        assert_almost_equal(self.inst.B().todense(), HybridMatrices.B,
                            decimal=16)

