from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.Utilities import Struct
from NewCode.tests.TestMeshes import FlatTet, TwoTets
from NewCode.ProxyList import ProxyList, ItemClassFactory
from NewCode import Mesh
from NewCode.DifferentialForm import Discretiser, DiscretiserEntities 
from NewCode.DifferentialForm import Permuter

class test_Permuter(TestCase):
    PermuterClass = Permuter.Permuter
    def setUp(self):
        self.elm = Struct(edgenos=N.arange(6), facenos=N.arange(4),
                          noDOFs=Struct(edge=Struct(perEntity=2, perElement=12),
                                        face=Struct(perEntity=1, perElement=4),
                                        element=16),
                          entityNames=('edge', 'face', 'vol'),
                          index=0, dead_elements=set()
                          )
        EntClass=ItemClassFactory(('freeNo',), 'DiscretiserEntityMock')
        self.geomEntities = Struct(edge=ProxyList(EntClass(
            dict(freeNo=N.array([-1, 10, 100, 1000, 10000, 100000])))
                                                  ),
                                   face=ProxyList(EntClass(
            dict(freeNo=N.array([1, 10, -1, 1000])))
                                                  ))
        self.geomEntities.edge.noFree = 5
        self.geomEntities.face.noFree = 3
        self.inst = self.PermuterClass(self.elm, self.geomEntities)
        
    def test_init(self):
        self.assertRaises(AssertionError, self.PermuterClass,
                          self.elm, ('node', 'garbage'))
        
        assert(self.geomEntities is self.inst.geomEntities)
        assert_equal(self.inst.noDOFs, dict(edge=10, face=3))
        assert_equal(self.inst.offsetDOFs, dict(edge=0,face=10))
        assert_equal(self.inst.totalDOFs, 13)

    def test_permuteElement(self):
        # Make the mesh consistent
        self.geomEntities.edge[:]._freeNo=N.array([-1,0,1,2,3,4])
        self.geomEntities.face[:]._freeNo=N.array([0,1,-1,2])
        inst = self.PermuterClass(self.elm, self.geomEntities)
        elperms = ((1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 15),
                   (0, 2, 4, 6, 8, 1, 3, 5,  7,  9, 10, 11, 12))
        self.inst.permuteElement(self.elm)
        assert_equal(self.inst.permuteElement(self.elm), elperms)
        
    def test_permuteElementPerEntity(self):
        elm, geomEntities, inst = self.elm, self.geomEntities, self.inst
        
        perms=dict(edge=N.array([-1,20,200,2000,20000,200000,-1,21,201,2001,20001,200001]),
                   face=N.array([1, 10, -1, 1000]))
        actual_perm=inst.permuteElementPerEntity(elm)
        # Fix all negative DOFs (i.e. constrained) to be -1
        for prm in actual_perm.values(): prm[prm < 0] = -1
        assert_equal(actual_perm, perms)

    def test_permuteEntities_LQTN(self):
        # I guess you could call this one a bit of an integration test
        from NewCode.DifferentialForm.BasisFunction import Oneform
        basisSet = Oneform.basis_set(2,mixed=True)
        mesh = Mesh.Mesh(TwoTets.listmesh)
        geomEntities = DiscretiserEntities.make_geomEntities(mesh, basisSet,
                                                             lambda x: True)
        disc = Discretiser.PformDiscretiser(
            1, mesh, geomEntities, self.PermuterClass, basisSet)
        p = disc.permuter
        # The permutations were verified by hand so is almost certainly correct
        (l, g) = p.permuteElement(disc.elements[0])
        assert_equal((l,g), (range(20), [0,2,4,6,8,10,1,3,5,7,9,11,
                                         18,20,22,24,19,21,23,25]))
        (l, g) = p.permuteElement(disc.elements[1])
        assert_equal((l,g), (range(20), [2,4,12,10,14,16,3,5,13,11,15,17,
                                         22,26,28,30,23,27,29,31]))

    def test_permuteEntities_QTCuN(self):
        # I guess you could call this one a bit of an integration test
        from NewCode.DifferentialForm.BasisFunction import Oneform
        basisSet = Oneform.basis_set(3 ,mixed=True)
        mesh = Mesh.Mesh(TwoTets.listmesh)
        geomEntities = DiscretiserEntities.make_geomEntities(mesh, basisSet,
                                                             lambda x: True)
        disc = Discretiser.PformDiscretiser(
            1, mesh, geomEntities, self.PermuterClass, basisSet)
        p = disc.permuter
        # These permutations weren't manually verified, but seem correct.
        (l, g) = p.permuteElement(disc.elements[0])
        assert_equal((l,g), (range(45), [
            0,  3,  6,  9, 12, 15,  1,  4,  7, 10, 13, 16,  2,  5,  8, 11, 14,
            17, 27, 33, 39, 45, 28, 34, 40, 46, 29, 35, 41, 47, 30, 36, 42, 48,
            31, 37, 43, 49, 32, 38, 44, 50, 69, 70, 71]))
        (l, g) = p.permuteElement(disc.elements[1])
        assert_equal((l,g), (range(45), [
            3,  6, 18, 15, 21, 24,  4,  7, 19, 16, 22, 25,  5,  8, 20, 17, 23,
            26, 39, 51, 57, 63, 40, 52, 58, 64, 41, 53, 59, 65, 42, 54, 60, 66,
            43, 55, 61, 67, 44, 56, 62, 68, 72, 73, 74])) 
    def test_newDOFs(self):
        dofs = self.inst.newDOFs()
        assert_equal(dofs, N.zeros(13))
        assert_equal(dofs.dtype, N.float64)
                     
class test_PermuterWithEntities(TestCase):
    PermuterClass = Permuter.Permuter
    freefun = staticmethod(
        lambda ent: not (len(ent.nodes) == 2 and ent.index == 8
                         or (len(ent.nodes)==3 and ent.index==5)
                         ))
    TestMesh = TwoTets
    
    def setUp(self):
        self.testMesh = self.TestMesh()

    def setup_disc(self, order, mixed=True):
        # I guess you could call this one a bit of an integration test
        from NewCode.DifferentialForm.BasisFunction import Oneform
        basisSet = Oneform.basis_set(order ,mixed=mixed)
        mesh = Mesh.Mesh(self.testMesh.listmesh)
        freefun = self.freefun
        geomEntities = DiscretiserEntities.make_geomEntities(mesh, basisSet,
                                                             freefun)
        disc = Discretiser.PformDiscretiser(
            1, mesh, geomEntities, self.PermuterClass, basisSet)
        p = disc.permuter
        return disc, p

    def test_permuteEntities_QTCuN(self):
        disc, p = self.setup_disc(3)
        ent_dofs = p.permuteElementEntities(disc.elements[0])
        # These permutations weren't manually verified, but seem correct.
        per_edge, per_face, per_vol = 3,6,3
        no_edgefuns, no_facefuns, no_volfuns = 8*per_edge, 6*per_face, 2*per_vol
        face_offset, vol_offset = no_edgefuns, no_edgefuns + no_facefuns
        assert_equal(ent_dofs['edge'], [N.arange(6), N.arange(18).reshape((6,per_edge))])
        assert_equal(ent_dofs['vol'], [0, [N.arange(per_vol)+vol_offset]])
        assert_equal(ent_dofs['face'], [N.arange(4),
                                        N.arange(24).reshape((4,per_face)) + face_offset])

        ent_dofs = p.permuteElementEntities(disc.elements[1])
        re = N.arange(per_edge)
        assert_equal(ent_dofs['edge'],
                     [N.arange(5), N.array([3+re, 6+re, 18+re, 15+re, 21+re])])
        assert_equal(ent_dofs['vol'], [0, [N.arange(3)+vol_offset+3]])
        rf = N.arange(per_face)
        assert_equal(ent_dofs['face'], [N.array([0,1,3]), N.array(
            [2*6+rf, 4*6+rf, 5*6+rf])+face_offset])

    def test_permuteGlobalEntities_QTCuN(self):
        disc, p = self.setup_disc(3)
        per_edge, per_face, per_vol = 3,6,3
        no_edgefuns, no_facefuns, no_volfuns = 8*per_edge, 6*per_face, 2*per_vol
        face_offset, vol_offset = no_edgefuns, no_edgefuns + no_facefuns
        re = N.arange(per_edge)
        edge_dofmap = N.array([re, re+per_edge, re+2*per_edge, re+3*per_edge, re+4*per_edge,
                              re+5*per_edge, re+6*per_edge, re+7*per_edge, re*0-1])
        rf = N.arange(per_face)
        face_dofmap = N.array([rf, rf+per_face, rf+2*per_face, rf+3*per_face, 
                              rf+4*per_face, rf*0-1, rf+5*per_face])
        face_dofmap[face_dofmap >= 0] += face_offset
        rv = N.arange(per_vol)
        vol_dofmap = N.array([rv, rv+per_vol])
        vol_dofmap[vol_dofmap >= 0] += vol_offset
        assert_equal(p.globalEntityPermutationTable('edge'), edge_dofmap)
        assert_equal(p.globalEntityPermutationTable('face'), face_dofmap)
        # Won't work until the volume entity type supports slicing
        #assert_equal(p.globalEntityPermutationTable('vol'), vol_dofmap)
