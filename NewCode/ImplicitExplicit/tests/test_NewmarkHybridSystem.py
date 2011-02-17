from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Meshes import BrickMesh
from NewCode.ImplicitExplicit.NewmarkHybridSystem import NewmarkHybridSystem
import HybridMatrices
from test_ImplicitExplicit import _FourBrickLinearHybridBlockMatrix

class test_NewmarkHybridSystem(_FourBrickLinearHybridBlockMatrix):
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)
        self.inst = NewmarkHybridSystem(self.mesh, self.implicit_beta)

    def test_init(self):
        self.inst.init_elgroups(self.implicitElements)
        assert_equal(self.inst.elgroups, self.elgroups)
        self.inst.init_group_freefuns(self.global_freefun)
        self.verify_group_freefuns(self.inst.group_freefuns)

    def test_blockmats(self):
        self.inst.init_elgroups(self.implicitElements)
        self.inst.init_group_freefuns(self.global_freefun)
        self.inst.init_discs()
        self.inst.init_block_matrices()
        self.inst.set_dt(self.dt)
        self.verify_block_mat('A', self.block_mats['A'],
                              self.inst.block_matrices, HybridMatrices)
        self.verify_block_mat('B', self.block_mats['B'],
                              self.inst.block_matrices, HybridMatrices)
