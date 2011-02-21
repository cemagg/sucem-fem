from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.DifferentialForm import Discretiser, DiscretiserEntities 
from NewCode.DifferentialForm import BasisFunction
from NewCode.DifferentialForm import DiscretiserElement
from NewCode import Mesh
from NewCode.tests.TestMeshes import FlatTet, TwoTets
from NewCode.Utilities import Struct
from NewCode.tests import xfail
from NewCode.ProxyList import ProxyList, ItemClassFactory

class test_PformDiscretiser(TestCase):
    Permuter = Discretiser.Permuter
    def setUp(self):
        self.mesh=Mesh.Mesh(FlatTet.listmesh)
        self.geomEdge = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Edge(attrs=self.mesh.edges.list_repr(),
                                             mesh=self.mesh, freefun=lambda x: True))
        self.geomFace = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Face(attrs=self.mesh.faces.list_repr(),
                                             mesh=self.mesh, freefun=lambda x: True))
        self.geomEntities=dict(edge=self.geomEdge, face=self.geomFace,
                               vol=DiscretiserEntities.FakeVolEntity())
        self.basisSet = BasisFunction.Oneform.basis_set(3)
        
    def test_init(self):
        PfD=Discretiser.PformDiscretiser
        Permuter = Discretiser.Permuter
        self._test_init(PfD, Permuter)

    def _test_init(self, PfD, PermuterClassOrInstance):
        mesh=self.mesh
        geomEntities = self.geomEntities
        basisSet = self.basisSet
        # Test that p != (1,2) raises an error
        self.assertRaises(NotImplementedError, PfD, 0, mesh, 
                          geomEntities, PermuterClassOrInstance, basisSet)
        inst = PfD(1, mesh, geomEntities, PermuterClassOrInstance, basisSet)
        self.inst = inst                # Store incase it is to be re-used
        self.assert_(isinstance(inst.elements[0], DiscretiserElement.OneformElement))
        self.assert_(inst.geomEntities is geomEntities)
        self.assert_(isinstance(inst.permuter, self.Permuter))
        from NewCode.DifferentialForm.DiscretiserMatrices import DiscretiserMatrices
        self.assert_(isinstance(inst.matrix, DiscretiserMatrices))
        self.assert_(inst.matrix.disc is self.inst)
        # Is PfD checking geometry entities are defined on the correct same mesh
        geomEntities['edge'].mesh = Mesh.Mesh(FlatTet.listmesh)
        self.assertRaises(
            AssertionError, PfD, 1, mesh, geomEntities, PermuterClassOrInstance, basisSet)

    def _test_permutation(self,inst):
        # QT/Cun element, all edges/faces unconstrained 
        assert_equal(inst.elements[0].permutation(), [
            range(45), [ 0,  3,  6,  9, 12, 15,  1,  4,  7, 10, 13, 16,  2,  5,  8, 11, 14,
                         17, 18, 24, 30, 36, 19, 25, 31, 37, 20, 26, 32, 38, 21, 27, 33, 39,
                         22, 28, 34, 40, 23, 29, 35, 41, 42, 43, 44]])
        dofs = inst.newDOFsArray()
        assert_equal(dofs, N.zeros(45))
        assert_equal(dofs.dtype, N.float64)
        assert_equal(inst.totalDOFs, 45)

    def test_permutation(self):
        inst = Discretiser.PformDiscretiser(
            1, self.mesh, self.geomEntities, Discretiser.Permuter, self.basisSet)
        self._test_permutation(inst)

    def test_basisType_D(self):
        inst = Discretiser.PformDiscretiser(
            1, self.mesh, self.geomEntities, Discretiser.Permuter, self.basisSet)
        inst_D = inst.D()
        # The basis functions of a 1-form should be i.t.o. the 4 covariant
        # component basis vectors
        assert_equal(len(inst.elements[0].refVals()[0][0]), 4)
        # The basis functions of a 2-form should be i.t.o. the 6 contravariant
        # component basis vectors
        assert_equal(len(inst_D.elements[0].refVals()[0][0]), 6)
        
class test_PformDiscretiser_D(test_PformDiscretiser):
    def _get_inst(self):
        permuter = Discretiser.PformDiscretiser(
            1, self.mesh, self.geomEntities,
            self.Permuter, self.basisSet).permuter
        return Discretiser.PformDiscretiser_D(
            1, self.mesh, self.geomEntities, permuter, self.basisSet)

    def test_init(self):
        PfD=Discretiser.PformDiscretiser_D
        permuter = Discretiser.PformDiscretiser(
            1, self.mesh, self.geomEntities,
            self.Permuter, self.basisSet).permuter
        self._test_init(PfD, permuter)
        
    def test_permutation(self):
        inst = self._get_inst()
        self._test_permutation(inst)
        
    def test_basisType_D(self):
        inst = self._get_inst()
        self.assertRaises(AttributeError, getattr, inst, 'D')
        

# class test_freefuns(TestCase, FlatTet):
#     def setUp(self):
        
#         self.mesh = Mesh.Mesh(self.listmesh)
#         # Make the boundary set more interesting, since the orignal mesh
#         # has only boundary edges/faces
#         self.mesh.edges.entity._onBoundary = set([1,2,4,5])
#         self.mesh.faces.entity._onBoundary = set([0,2,3])
        
#     def test_E_PEC_free_edge(self):
#         edges = self.mesh.edges
#         desired_free = [not edge.onBoundary for edge in edges]
#         # scalar operation
#         assert_array_equal(
#             [Discretiser.E_PECboundary_free_edge(edge) for edge in edges],
#             desired_free)
#         # array operation
#         assert_array_equal(Discretiser.E_PECboundary_free_edge(edges[:]),
#                            desired_free)
        

#     def test_B_PECboundary_free_face(self):
#         faces = self.mesh.faces
#         desired_free = [not face.onBoundary for face in faces]
#         # scalar operation
#         assert_array_equal(
#             [Discretiser.B_PECboundary_free_face(face) for face in faces],
#             desired_free)
#         # array operation
#         assert_array_equal(Discretiser.B_PECboundary_free_face(faces[:]),
#                            desired_free)

#     def test_EB_PEC_with_Discretiser(self):
#         disc1form = Discretiser.OneformDiscretiser(
#             self.mesh,
#             free_edge=Discretiser.E_PECboundary_free_edge)
#         # 2 unconstrained edges
#         assert_equal(disc1form.edges.nodofs, 2)
#         assert_array_equal(disc1form.edges[:].dofno, [0,-1,-1,1,-1,-1])
#         assert_array_equal(disc1form.edges[:].isfree,
#                            [True, False, False, True, False, False])

#         disc2form = Discretiser.TwoformDiscretiser(
#             self.mesh,
#             free_face=Discretiser.B_PECboundary_free_face)
#         # 1 unconstrained face
#         assert_equal(disc2form.faces.nodofs, 1)
#         assert_array_equal(disc2form.faces[:].dofno, [-1,0,-1,-1])
#         assert_array_equal(disc2form.faces[:].isfree, [False, True, False, False])

