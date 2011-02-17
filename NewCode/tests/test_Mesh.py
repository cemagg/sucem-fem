from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal,\
     assert_array_almost_equal,assert_almost_equal, assert_equal
from numpy import array, float64, int32, sqrt
import sys
#
# Local Imports
#
from NewCode import Mesh
from NewCode.tests import xfail
import NewCode.tests.TestMeshes 
from NewCode.tests.TestMeshes import FlatTet, TwoTets, InscribedTetMesh


class test_Mesh(NumpyTestCase, FlatTet):
    def setUp(self):
        FlatTet.setUp(self)                     
        self.mesh=Mesh.Mesh(self.listmesh)
        pass

    def test_CovBaseVecs(self):
        assert_array_almost_equal(self.mesh.elements[0].covBaseVecs().T,
                                  [[-1.5, 0.5,-0.5], # Calculated with Maxima
                                   [-1.5,-0.5, 0.5],
                                   [ 1.5, 1.5, 1.5],
                                   [ 1.5,-1.5,-1.5]],
                                  decimal=15)
        pass

    def test_FaceBasisVecs(self):
        # Calculated with Maxima 
        facebases = [[[-1.5, 3.0,-1.5],[-1.5,-1.5, 3.0],[ 0.0, 1.5, 1.5]],
                     [[ 1.5,-1.5, 3.0],[ 1.5, 3.0,-1.5],[ 0.0, 1.5, 1.5]],
                     [[ 0.0, 4.5,-4.5],[ 1.5, 3.0,-1.5],[ 1.5, 1.5,-3.0]],
                     [[ 0.0, 4.5,-4.5],[-1.5, 1.5,-3.0],[-1.5, 3.0,-1.5]]]

        assert_array_almost_equal([self.mesh.elements[0].FaceBasisVecs(i)
                                   for i in range(4)],
                                  facebases, decimal=15)

        pass

    def test_outwardNormals(self):
        desired = [[-N.sqrt(3)/3,1/N.sqrt(3),1/N.sqrt(3)],
                   [-N.sqrt(3)/3,-N.sqrt(3)/3,-N.sqrt(3)/3],
                   [3/N.sqrt(11),1/N.sqrt(11),-1/N.sqrt(11)],
                   [3/N.sqrt(11),-1/N.sqrt(11),1/N.sqrt(11)]]
        assert_almost_equal(self.mesh.elements[0].outwardNormals(), desired,
                            decimal=16)

    def test_Element(self):
        elms = self.mesh.elements
        assert_array_equal([elm.facenos for elm in elms], self.ElementFaces)
        assert_array_equal([elm.edgenos for elm in elms], self.ElementEdges)
        assert_array_equal([elm.nodes for elm in elms], self.ElementNodes)
        for elm in elms:
            assert_array_equal(elm.nodeCoords,
                               [self.Nodes[i] for i in elm.nodes])
    def test_Face(self):
        faces = self.mesh.faces
        assert_array_equal([face.nodes for face in faces], self.FaceNodes)
        assert_array_equal([face.connect2elem for face in faces],
                           self.FaceConnect2Elem)
        assert_array_equal(faces[:].edgenos, self.FaceEdges)
        for face in faces:
            assert_array_equal(face.nodeCoords,
                               [self.Nodes[i] for i in face.nodes])
    def test_faces_onBoundary(self):
        faces=self.mesh.faces
        assert_array_equal(faces.onBoundary.nodes,
                           [face.nodes for face in faces if face.onBoundary])
        assert_array_equal(faces.onBoundary.nodes,
                           [face.nodes for face in faces.onBoundary])
        assert_equal(set(faces.onBoundary.index), faces.boundarySet)
        
    def test_Edge(self):
        edges = self.mesh.edges
        assert_array_equal([edge.nodes for edge in edges], self.EdgeNodes)
        for edge in edges:
            assert_array_equal(edge.nodeCoords,
                               [self.Nodes[i] for i in edge.nodes])

    def test_edges_onBoundary(self):
        edges=self.mesh.edges
        assert_array_equal(edges.onBoundary.nodes,
                           [edge.nodes for edge in edges if edge.onBoundary])
        assert_array_equal(edges.onBoundary.nodes,
                           [edge.nodes for edge in edges.onBoundary])
        

class test_MeshKDTree(NumpyTestCase):
    def setUp(self):
        self.mesh = Mesh.MeshWithKDTree(InscribedTetMesh.listmesh)

    def test_kdTree(self):
        kdt1 = self.mesh.kdTree
        kdt2 = self.mesh.kdTree
        self.assert_(kdt1 is kdt2)

    def test_findNodesRadius(self):
        assert_equal(N.unique(self.mesh.findNodesRadius([0,0,0], 10)),
                     N.arange(8))
        assert_equal(self.mesh.findNodesRadius([1./6, 1./6, -1./6], 0.01),
                     N.array([4]))

    def test_findClosestNode(self):
        assert_equal(self.mesh.findClosestNode([1./2, 1./2, 1./2-0.05]), 1)

    def test_maxEdgeLength(self):
        assert_almost_equal(self.mesh.maxEdgeLength, 1.41421356237,
                            decimal=11)
        
class test_Mesh_TwoTets(NumpyTestCase):
    def setUp(self):
        self.listmesh = TwoTets.listmesh
        self.mesh=Mesh.Mesh(self.listmesh)

    def test_Face_area(self):
        desired = array([sqrt(3)/2,sqrt(3)/2,sqrt(3)/2,sqrt(3)/2,1/2,1/2,1/2])
        assert_array_almost_equal([f.area() for f in self.mesh.faces], desired,
                                  decimal=16)
        
    def test_Element_Midpt(self):
        assert_array_equal([elm.midpoint for elm in self.mesh.elements],
                           [self.listmesh['Nodes'][elnodes].sum(axis=0)/4.
                            for elnodes in self.listmesh['ElementNodes']])
        assert_array_equal(self.mesh.elements[:].midpoint,
                           [self.listmesh['Nodes'][elnodes].sum(axis=0)/4.
                            for elnodes in self.listmesh['ElementNodes']])
        
class test_Mesh_Exeptions(NumpyTestCase):
    def setUp(self):
        self.listmesh={}
        self.listmesh['ElementNodes']=array([[1,2,3,4]]) - 1
        self.listmesh['ElementEdges']=array([[1,2,3,4,5,6]]) - 1
        self.listmesh['ElementFaces']=array([[1,2,3,4]]) - 1
        self.listmesh['FaceNodes']=array([[1,2,3],
                                      [1,2,4],
                                      [1,3,4],
                                      [2,3,4]]) - 1
        self.listmesh['EdgeNodes']=array([[1,2],
                                      [1,3],
                                      [1,4],
                                      [2,3],
                                      [2,4],
                                      [3,4]]) - 1
        self.listmesh['Nodes']=array([[-0.5, 0.5,-0.5],
                                  [ 0.5, 0.5, 0.5],
                                  [ 0.5,-0.5,-0.5],
                                  [-0.5,-0.5, 0.5]], float64)
        
    def test_element_node_sort(self):
        self.listmesh['ElementNodes']=array([[1,3,2,4]])-1
        self.assertRaises(AssertionError, Mesh.Mesh, self.listmesh)
    
    def test_face_node_sort(self):
        self.listmesh['FaceNodes'][0]=array([3,2,1]) - 1
        self.assertRaises(AssertionError, Mesh.Mesh, self.listmesh)

    def test_edge_node_sort(self):
        self.listmesh['EdgeNodes'][0]=array([2,1]) - 1
        self.assertRaises(AssertionError, Mesh.Mesh, self.listmesh)

class test_Mesh_Connection(NumpyTestCase, InscribedTetMesh):
    def setUp(self):
        self.mesh = Mesh.Mesh(self.listmesh)
        # Arrays sourced from test_eMAGUSImport.py. Perhaps they should be put
        # somewhere common :)
        self.connect2elem = array([[0, 0, 8, 7], 
                                   [0, 0, 8, 9], 
                                   [0, 0, 9, 7], 
                                   [0, 0, 8, 10],
                                   [0, 0, 10, 9],
                                   [0, 0, 10, 7],
                                   [1, 6, 3, 11],
                                   [4, 2, 1, 11],
                                   [2, 5, 3, 11],
                                   [4, 5, 6, 11],
                                   [8, 10, 9, 7]], int32) - 1
        self.connect2face = array([[0, 0, 3, 1],
                                   [0, 0, 2, 1],
                                   [0, 0, 3, 3],
                                   [0, 0, 1, 1],
                                   [0, 0, 2, 2],
                                   [0, 0, 3, 2],
                                   [4, 4, 4, 4],
                                   [3, 3, 3, 1],
                                   [4, 4, 3, 3],
                                   [4, 3, 3, 2],
                                   [4, 4, 4, 4]], int32) - 1
    def test_element_face_connections(self):
        elements = self.mesh.elements
        assert_array_equal(elements[:].connect2elem,self.connect2elem)
        assert_array_equal(elements[:].connect2face,self.connect2face)

    def test_node_element_connections(self):
        desired_connections = [[1, 2, 4, 8],
                               [4, 5, 6, 10],
                               [2, 3, 5, 9],
                               [1, 3, 6, 7],
                               [2, 4, 5, 8, 9, 10, 11],
                               [1, 4, 6, 7, 8, 10, 11],
                               [1, 2, 3, 7, 8, 9, 11],
                               [3, 5, 6, 7, 9, 10, 11]]
        assert_equal([(self.mesh.nodeElementConnections[i] + 1).tolist()
                      for i in range(len(self.mesh.nodes))],
                     desired_connections)

    def test_edge_nodemap(self):
        em = self.mesh.edges.nodemap
        assert_equal (len(em), len(self.mesh.edges))
        for k,v in em.iteritems():
            assert_equal(N.array(k, N.int32), self.mesh.edges[v].nodes)

    def test_face_nodemap(self):
        em = self.mesh.faces.nodemap
        assert_equal (len(em), len(self.mesh.faces))
        for k,v in em.iteritems():
            assert_equal(N.array(k, N.int32), self.mesh.faces[v].nodes)

    def test_element_nodemap(self):
        em = self.mesh.elements.nodemap
        assert_equal (len(em), len(self.mesh.elements))
        for k,v in em.iteritems():
            assert_equal(N.array(k, N.int32), self.mesh.elements[v].nodes)

    def test_face_onBoundary(self):
        faces=self.mesh.faces
        assert_array_equal(faces[:].onBoundary, self.BoundaryFaces)
        # To test that the iter_over_slice code works
        assert_array_equal([face.onBoundary for face in faces[3:7]],
                           self.BoundaryFaces[3:7])
        assert_equal(faces[10].onBoundary, self.BoundaryFaces[10])
        assert_array_equal(faces[[1,3,5,10]].onBoundary, self.BoundaryFaces[[1,2,5,10]])

    def test_edge_onBoundary(self):
        edges=self.mesh.edges
        assert_array_equal(edges[:].onBoundary, self.BoundaryEdges)

class test_EntityClasses(NumpyTestCase, FlatTet):
    def setUp(self):
        FlatTet.setUp(self)                     
        self.mesh=Mesh.Mesh(self.listmesh)
        mesh = self.mesh
        self.el = mesh.elements[0]
        self.face = mesh.faces[0]
        self.edge = mesh.edges[0]
    def test_proxy_attrs(self):
        assert_equal(self.el.proxy_attrs,
                     ('facenos', 'edgenos', 'nodes', 'nodeCoords',
                      'connect2elem', 'connect2face'))
        assert_equal(self.face.proxy_attrs,
                     ('edgenos', 'nodes', 'nodeCoords', 'connect2elem', 'onBoundary'))
        assert_equal(self.edge.proxy_attrs,
                     ('nodes', 'nodeCoords', 'connect2elem', 'onBoundary'))

    def test_element_list_repr(self):
        listmesh = self.listmesh
        el_list_repr = {'facenos': listmesh['ElementFaces'],
                        'edgenos': listmesh['ElementEdges'],
                        'nodes': listmesh['ElementNodes'],
                        'connect2elem': listmesh['ElementConnect2Elem'],
                        'connect2face': listmesh['ElementConnect2Face'],
                        'nodeCoords': listmesh['Nodes']}
        assert_equal(el_list_repr, self.el.list_repr())

    def test_face_list_repr(self):
        listmesh = self.listmesh
        face_list_repr = {'edgenos': self.FaceEdges,
                          'nodes': listmesh['FaceNodes'],
                          'connect2elem': listmesh['FaceConnect2Elem'],
                          'nodeCoords': listmesh['Nodes'],
                          }
        face_list_repr['onBoundary'] = set(i for i, onBoundary
                                           in enumerate(self.BoundaryFaces)
                                           if onBoundary)
        assert_equal(face_list_repr, self.face.list_repr())

    def test_edge_list_repr(self):
        listmesh = self.listmesh
        edge_list_repr = {'nodes': listmesh['EdgeNodes'],
                          'nodeCoords': listmesh['Nodes'],
                          'connect2elem': listmesh['EdgeConnect2Elem'],
                          'onBoundary' : set(
            i for i, onBoundary in enumerate(self.BoundaryEdges)
            if onBoundary)
                          ,}
        assert_equal(edge_list_repr, self.edge.list_repr())

class test_find_face_edges(NumpyTestCase, FlatTet):
    def setUp(self):
        self.faces = Mesh.ProxyFaces(Mesh.Face({
            'nodes': self.listmesh['FaceNodes'],
            'nodeCoords': self.listmesh['FaceNodes'],
            'connect2elem': self.listmesh['FaceConnect2Elem'],
            })
                                     )
        self.elements = Mesh.ProxyElements(Mesh.Element({
            'facenos': self.listmesh['ElementFaces'],
            'edgenos': self.listmesh['ElementEdges'],
            'nodes': self.listmesh['ElementNodes'],
            'nodeCoords' : self.listmesh['ElementNodes'],
            'connect2elem': self.listmesh['ElementConnect2Elem'],
            'connect2face': self.listmesh['ElementConnect2Face']
            })
                                           )
    def test_find_face_edges(self):
        face_edges = Mesh.find_face_edges(self.faces, self.elements)
        desired_face_edges = array([[1,4,2],
                                    [1,5,3],
                                    [2,6,3],
                                    [4,6,5]]) - 1
        assert_array_equal(face_edges, desired_face_edges)
            
