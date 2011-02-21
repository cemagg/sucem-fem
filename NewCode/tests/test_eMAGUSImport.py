from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, \
     assert_equal, assert_array_almost_equal
from numpy import array, float64, int32
import sys, os
#
# Local Imports
#
from NewCode import Mesh, eMAGUSImport
import NewCode.tests.TestMeshes as TestMeshes

InscribedTetMesh = TestMeshes.InscribedTetMesh


class test_eMAGUSImport(TestCase, InscribedTetMesh):
    def setUp(self):
        InscribedTetMesh.setUp(self)
        self.benchmark_mesh=Mesh.Mesh(self.listmesh)
        eMAGUSImport.init('NewCode/tests/testdata/test_eMAGUSImport')
        self.testlistmesh = eMAGUSImport.get_listmesh()
        self.testmesh = Mesh.Mesh(self.testlistmesh)
        pass
    
    def test_femmesh_filename(self):
        assert_equal(self.testlistmesh['FemmeshFilename'], 'test_tet.femmesh')

    def test_node_elements(self):
        assert_equal(self.testlistmesh['NodeConnect2Element'],
                     array([1, 2, 4, 8, 4, 5,  6,  10, 2, 3, 5, 9,  1, 3, 6, 7,
                            2, 4, 5, 8, 9, 10, 11, 1,  4, 6, 7, 8, 10, 11,
                            1, 2, 3, 7, 8, 9,  11, 3,  5, 6, 7, 9, 10, 11],
                           int32) - 1)
        assert_equal(self.testlistmesh['NodeConnect2ElementPtr'],
                     array([1, 5, 9, 13, 17, 24, 31, 38, 45], int32) - 1)
    

    def test_CovBaseVecs(self):
        grad_lambdas = array(
            [[[-1.5, 0.5, -0.5], [-1.5, -0.5, 0.5], [1.5, 1.5, 1.5], [1.5, -1.5, -1.5]],
             [[-0.5, 0.5, -1.5], [0.5, -0.5, -1.5], [1.5, 1.5, 1.5], [-1.5, -1.5, 1.5]],
             [[0.5, -1.5, -0.5], [-0.5, -1.5, 0.5], [-1.5, 1.5, -1.5], [1.5, 1.5, 1.5]],
             [[-0.5, 1.5, -0.5], [0.5, 1.5, 0.5], [1.5, -1.5, -1.5], [-1.5, -1.5, 1.5]],
             [[1.5, 0.5, 0.5], [1.5, -0.5, -0.5], [-1.5, 1.5, -1.5], [-1.5, -1.5, 1.5]],
             [[0.5, 0.5, 1.5], [-0.5, -0.5, 1.5], [-1.5, 1.5, -1.5], [1.5, -1.5, -1.5]],
             [[-0.75, -0.75, 0.75],[-0.75, 2.25, 0.75],[-0.75, -0.75, -2.25], [2.25, -0.75, 0.75]],
             [[-0.75, 0.75, -0.75],[2.25, 0.75, -0.75],[-0.75, 0.75, 2.25], [-0.75, -2.25, -0.75]],
             [[0.75, -0.75, -0.75],[0.75, 2.25, -0.75],[-2.25, -0.75, -0.75],[0.75, -0.75, 2.25]],
             [[0.75, 0.75, 0.75],[0.75, 0.75, -2.25],[-2.25, 0.75, 0.75],[0.75, -2.25, 0.75]],
             [[1.5, 1.5, -1.5], [-1.5, 1.5, 1.5], [-1.5, -1.5, -1.5], [1.5, -1.5, 1.5]]],
            float64)
        assert_almost_equal(grad_lambdas,
                            [el.covBaseVecs().transpose() for el in self.testmesh.elements],
                            decimal=14)
        
        pass

    def test_Element(self):
        elms = self.testmesh.elements
        assert_array_equal([elm.facenos for elm in elms], self.ElementFaces)
        assert_array_equal([elm.edgenos for elm in elms], self.ElementEdges)
        assert_array_equal([elm.nodes for elm in elms], self.ElementNodes)
        for elm in elms:
            assert_array_almost_equal(elm.nodeCoords,
                                      [self.Nodes[i] for i in elm.nodes],
                                      decimal=16)
          

    def test_element_faceconnections(self):
        # These arrays were obtained from an eMAGUS run using the script
        # read_element_face_conn_data.py
        
        testmesh = self.testlistmesh
        listmesh = self.listmesh
        assert_array_equal(testmesh['ElementConnect2Elem'], listmesh['ElementConnect2Elem'])
        assert_array_equal(testmesh['ElementConnect2Face'], listmesh['ElementConnect2Face'])
        
    def test_faceconnectelem(self):
        faceconnectelem = self.testlistmesh['FaceConnect2Elem']
        assert_array_equal(faceconnectelem, self.listmesh['FaceConnect2Elem'])

    def test_edgeconnectelem(self):
        faceconnectelem = self.testlistmesh['EdgeConnect2Elem']
        assert_array_equal(faceconnectelem, self.listmesh['EdgeConnect2Elem'])

    def test_Face(self):
        faces = self.testmesh.faces
        assert_array_equal([face.nodes for face in faces], self.FaceNodes)
        for face in faces:
            assert_array_almost_equal(face.nodeCoords,
                                       [self.Nodes[i] for i in face.nodes],
                                       decimal=16)

    def test_Edge(self):
        edges = self.testmesh.edges
        assert_array_equal([edge.nodes for edge in edges], self.EdgeNodes)
        for edge in edges:
            assert_array_almost_equal(edge.nodeCoords,
                                      [self.Nodes[i] for i in edge.nodes],
                                      decimal=16)
