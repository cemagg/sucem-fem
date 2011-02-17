from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests.BrickMeshes import TestBrick
from NewCode import ImplicitExplicit as IE
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import BrickDiscretiser, allfree
from NewCode.Integration import BrickTrapzIntegrator

import HybridMatrices

class FourBricksLinear(TestBrick):
    listmesh = {
        'GridDimension' : N.array([2,5,2], N.int32),
        'GridStepSize' : N.array([1,2,3], N.float64),
        'GridOffset' : N.array([0,0,0], N.float64),
        }
    nodes = N.array([[0,0,0], [0,0,3],
                     [0,2,0], [0,2,3],
                     [0,4,0], [0,4,3],
                     [0,6,0], [0,6,3],
                     [0,8,0], [0,8,3],
                     [1,0,0], [1,0,3],
                     [1,2,0], [1,2,3],
                     [1,4,0], [1,4,3],
                     [1,6,0], [1,6,3],
                     [1,8,0], [1,8,3]], N.float64)

class _FourBrickLinearTestInfo(NumpyTestCase):
    TestMesh = FourBricksLinear
    implicitElements = set([0,1])
    dt = 1/33
    DiscretiserModule = BrickDiscretiser
    implicit_beta = 0.25
    elgroups = Struct(ii=set([0]), ie=set([1]), ei=set([2]), ee=set([3]))
    global_freefun = staticmethod(allfree)
    edgesets = dict(a=set([0,1,10,11,18,19,26,31]),
                    b=set([2,3,12,13,20,21,27,32]),
                    c=set([4,5,28,33]),
                    d=set([6,7,14,15,22,23,29,34]),
                    e=set([8,9,16,17,24,25,30,35]))
    facesets = dict(a=set([0,4,8,13,14]),
                    b=set([1,5,9,15,16]),
                    c=set([10]),
                    d=set([2,6,11,17,18]),
                    e=set([3,7,12,19,20]))
    volsets = dict(a=set([0]),
                   b=set([1]),
                   c=set([]),
                   d=set([2]),
                   e=set([3]))

    def verify_group_freefuns(self, group_freefuns):
        edgesets = self.edgesets ; facesets = self.facesets
        volsets = self.volsets
        for k,(freefun, freefun_vol) in group_freefuns.iteritems():
            t_edgeset = set([i for i,e in enumerate(self.mesh.edges) if freefun(e)])
            assert_equal(t_edgeset, edgesets[k])
            t_faceset = set([i for i,f in enumerate(self.mesh.faces) if freefun(f)])
            assert_equal(t_faceset, facesets[k])
            t_volset = set([i for i,f in enumerate(self.mesh.elements)
                            if freefun_vol(f)])
            assert_equal(t_volset, volsets[k])

    
class test_DOF_GroupNumbering(_FourBrickLinearTestInfo):
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)

    def test_groupFreefuns(self):
        group_freefuns = IE.gen_group_freefuns(
            self.mesh, self.elgroups, self.global_freefun)
        self.verify_group_freefuns(group_freefuns)
        
    def test_elGroups(self):
        assert_equal(IE.gen_elgroups(self.mesh, self.implicitElements),
                     self.elgroups)

class _FourBrickLinearHybridBlockMatrix(_FourBrickLinearTestInfo):
    block_mats = {'A':(('aa', 12), ('bb', 12), ('ab', 16), ('bc', 16),
                       ('cc', 12), ('dd', 12), ('ee', 12),),
                  'B':(('aa', 12), ('ab', 16), ('bb', 11), ('bc', 16),
                       ('cc', 12), ('cd', 16), ('dd', 11),
                       ('de', 16), ('ee', 12))
                  }

    @staticmethod
    def verify_block_mat(matname, subblocks, submats_calc, submats_ref):
        for submat, decim_tol in subblocks:
            mat = matname+'_'+submat
            print mat
            assert_almost_equal(getattr(submats_calc, mat)().todense(),
                                getattr(submats_ref, mat), decim_tol)


