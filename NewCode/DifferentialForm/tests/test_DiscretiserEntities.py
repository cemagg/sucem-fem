from __future__ import division

from numpy import zeros, array, arange, alltrue, float64, int32
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode.tests.TestMeshes import InscribedTetMesh, FlatTet
from NewCode import  Mesh
from NewCode.DifferentialForm import constrained_on_boundary
from NewCode.DifferentialForm import DiscretiserEntities

class test_DiscretiserEntities(TestCase, InscribedTetMesh):
    def setUp(self):
        self.mesh = Mesh.Mesh(self.listmesh)

    def test_Edge(self):
        edges = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Edge(mesh=self.mesh, freefun=lambda x: True,
                                             attrs=self.mesh.edges.list_repr()
                                             ))
        noEdges = self.mesh.noEdges
        list_repr = edges.list_repr()
        desired_list_repr = {'nodes': self.listmesh['EdgeNodes'],
                             'nodeCoords': self.listmesh['Nodes'],
                             'connect2elem': self.listmesh['EdgeConnect2Elem'],
                             'freeNo' : arange(noEdges, dtype=int32),
                             }
        desired_list_repr['onBoundary'] = set(i for i, onBoundary
                                              in enumerate(self.BoundaryEdges)
                                              if onBoundary)
        assert_equal(desired_list_repr, list_repr)

    def test_Face(self):
        faces = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Face(mesh=self.mesh, freefun=lambda x: True,
                                             attrs=self.mesh.faces.list_repr()
                                             ))
        noFaces = self.mesh.noFaces
        list_repr = faces.list_repr()
        desired_list_repr = {'nodes': self.listmesh['FaceNodes'],
                             'nodeCoords': self.listmesh['Nodes'],
                             'connect2elem': self.listmesh['FaceConnect2Elem'],
                             'freeNo' : arange(noFaces, dtype=int32),
                             'edgenos' : self.FaceEdges
                             }
        desired_list_repr['onBoundary'] = set(i for i, onBoundary
                                              in enumerate(self.BoundaryFaces)
                                              if onBoundary)
        assert_equal(desired_list_repr, list_repr)        

class test_freefuns(TestCase, FlatTet):
    def setUp(self):
        self.mesh = Mesh.Mesh(self.listmesh)
        # Make the boundary set more interesting, since the orignal mesh
        # has only boundary edges/faces
        self.mesh.edges.entity._onBoundary = set([1,2,4,5])
        self.mesh.faces.entity._onBoundary = set([0,2,3])
        
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
        edges = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Edge(mesh=self.mesh,
                                             freefun=constrained_on_boundary,
                                             attrs=self.mesh.edges.list_repr()
                                             )
            )
        # 2 unconstrained edges
        assert_equal(edges.noFree, 2)
        assert_array_equal(edges[:].freeNo, [0,-1,-1,1,-1,-1])
        assert_array_equal(edges[:].isFree,
                           [True, False, False, True, False, False])

    def test_constrained_on_boundary_face_freeNo(self):
        faces = DiscretiserEntities.DiscretiserEntityList(
            DiscretiserEntities.Face(mesh=self.mesh,
                                             freefun=constrained_on_boundary,
                                             attrs=self.mesh.faces.list_repr()
                                             ))
        # 1 unconstrained face
        assert_equal(faces.noFree, 1)
        assert_array_equal(faces[:].freeNo, [-1,0,-1,-1])
        assert_array_equal(faces[:].isFree, [False, True, False, False])
