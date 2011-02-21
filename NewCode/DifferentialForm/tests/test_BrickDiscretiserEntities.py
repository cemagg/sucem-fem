from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.tests.BrickMeshes import OneBrick, TwoBricks, FourBricks
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import constrained_on_boundary
from NewCode.DifferentialForm import BrickDiscretiserEntities

class test_DiscretiserEntities(TestCase):
    TestMesh = FourBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        self.mesh = BrickMesh.Mesh(self.listmesh)

    def test_Edge(self):
        edges = BrickDiscretiserEntities.DiscretiserEntityList(
            BrickDiscretiserEntities.Edge(mesh=self.mesh, freefun=lambda x: True,
                                             attrs=self.mesh.edges.list_repr()
                                             ))
        noEdges = self.mesh.noEdges
        list_repr = edges.list_repr()
        desired_list_repr = {'nodes': self.testMesh.edgeNodes,
                             'nodeCoords': self.testMesh.nodes,
                             'connect2elem': self.testMesh.edgeConnect2Elem,
                             'freeNo' : N.arange(noEdges, dtype=N.int32),
                             }
        desired_list_repr['onBoundary'] = set(i for i, onBoundary
                                              in enumerate(self.testMesh.boundaryEdges)
                                              if onBoundary)
        assert_equal(desired_list_repr, list_repr)

    def test_Face(self):
        faces = BrickDiscretiserEntities.DiscretiserEntityList(
            BrickDiscretiserEntities.Face(mesh=self.mesh, freefun=lambda x: True,
                                             attrs=self.mesh.faces.list_repr()
                                             ))
        noFaces = self.mesh.noFaces
        list_repr = faces.list_repr()
        desired_list_repr = {'nodes': self.testMesh.faceNodes,
                             'nodeCoords': self.testMesh.nodes,
                             'connect2elem': self.testMesh.faceConnect2Elem,
                             'freeNo' : N.arange(noFaces, dtype=N.int32),
                             'edgenos' : self.testMesh.faceEdges
                             }
        desired_list_repr['onBoundary'] = set(i for i, onBoundary
                                              in enumerate(self.testMesh.boundaryFaces)
                                              if onBoundary)
        assert_equal(desired_list_repr, list_repr)        


class test_freefuns(TestCase):
    TestMesh = OneBrick
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        self.mesh = BrickMesh.Mesh(self.listmesh)
        # Make the boundary set more interesting, since the orignal mesh
        # has only boundary edges/faces
        self.mesh.edges.entity._onBoundary = set([1,2,4,5,6,7,8,9,10,11])
        self.mesh.faces.entity._onBoundary = set([0,2,3,4,5])
        
    def test_constrained_on_boundary_edge(self):
        edges = self.mesh.edges
        desired_free = [not edge.onBoundary for edge in edges]
        # scalar operation
        assert_array_equal(
            [constrained_on_boundary(edge) for edge in edges],
            desired_free)
        # array operation
        assert_array_equal(constrained_on_boundary(edges[:]),
                           desired_free)
 
    def test_constrained_on_boundary_face(self):
        faces = self.mesh.faces
        desired_free = [not face.onBoundary for face in faces]
        # scalar operation
        assert_array_equal(
            [constrained_on_boundary(face) for face in faces],
            desired_free)
        # array operation
        assert_array_equal(constrained_on_boundary(faces[:]),
                           desired_free)

    def test_constrained_on_boundary_edge_freeNo(self):
        edges = BrickDiscretiserEntities.DiscretiserEntityList(
            BrickDiscretiserEntities.Edge(mesh=self.mesh,
                                             freefun=constrained_on_boundary,
                                             attrs=self.mesh.edges.list_repr()
                                             )
            )
        # 2 unconstrained edges
        assert_equal(edges.noFree, 2)
        assert_array_equal(edges[:].freeNo, [0,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1])
        assert_array_equal(edges[:].isFree,
                           [True,  False, False, True,  False, False,
                            False, False, False, False, False, False])

    def test_constrained_on_boundary_face_freeNo(self):
        faces = BrickDiscretiserEntities.DiscretiserEntityList(
            BrickDiscretiserEntities.Face(mesh=self.mesh,
                                             freefun=constrained_on_boundary,
                                             attrs=self.mesh.faces.list_repr()
                                             ))
        # 1 unconstrained face
        assert_equal(faces.noFree, 1)
        assert_array_equal(faces[:].freeNo, [-1,0,-1,-1,-1,-1])
        assert_array_equal(faces[:].isFree, [False, True, False, False, False, False])
