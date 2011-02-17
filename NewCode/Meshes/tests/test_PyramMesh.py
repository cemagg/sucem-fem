from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal

from NewCode import ProxyList
from NewCode.Meshes import PyramMesh, BrickMesh
from NewCode.tests.PyramMeshes import SixPyram, TwelvePyram

class test_GridGeneration(NumpyTestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()

    def test_NodeCoordinates(self):
        nodes = self.nodes
        assert_equal(nodes, self.testMesh.nodes)

    def test_nodeIjk(self):
        nodes = self.nodes
        assert_equal([nodes.nodeIjk(i) for i in range(nodes.pyram_offset)],
                      self.testMesh.nodesIjk)
        assert_equal([nodes.nodeNo(*nodes.nodeIjk(i)) for i in range(len(nodes))],
                     range(len(nodes)))
    @property
    def nodes(self):
        lm = self.testMesh.listmesh
        return PyramMesh.PyramGridNodes(lm['GridDimension'],
                                          lm['GridStepSize'],
                                          offset=lm['GridOffset'])
    

class _StructuredPyramEntityNodeGeneration(NumpyTestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()

    @property
    def edge_nodes(self):
        return PyramMesh.StructuredPyramEdges(self.testMesh.listmesh['GridDimension'])
    def test_EdgeNodeNumbers(self):
        edge_nodes = self.edge_nodes
        assert_equal(edge_nodes, self.testMesh.edgeNodes)
        assert_equal(edge_nodes.pyram_offset, self.testMesh.pyram_edge_offset)

    @property
    def element_nodes(self):
        base_faces = PyramMesh.StructuredPyramBaseFaces(
            self.testMesh.listmesh['GridDimension'])
        return PyramMesh.StructuredPyramNodes(self.testMesh.listmesh['GridDimension'],
                                              base_faces)
    def test_ElementNodeNumbers(self):
        elm_nodes = self.element_nodes
        assert_array_equal(elm_nodes, self.testMesh.elementNodes)

    @property
    def apexface_nodes(self):
        brick_els = ProxyList.ProxyList(BrickMesh.Element(dict(
            nodes=self.testMesh.brickElementNodes,
            nodeCoords=self.testMesh.nodes,
            edgenos=self.testMesh.brickElementEdges,
            gridStepSize=self.testMesh.listmesh['GridStepSize'])))
        brick_edges = ProxyList.ProxyList(BrickMesh.Edge(dict(
            nodes=self.testMesh.edgeNodes, nodeCoords=self.testMesh.nodes,
            gridStepSize=self.testMesh.listmesh['GridStepSize'])))
        return PyramMesh.StructuredApexfaceNodes(
            brick_els, brick_edges, self.testMesh.pyram_node_offset)
    
    def test_apexface_nodes(self):
        assert_array_equal(self.apexface_nodes, self.testMesh.apexfaceNodes)
        
class test_StructuredPyramEntityNodeGeneration(_StructuredPyramEntityNodeGeneration):
    pass

class test_PyramElementEntityNumbers(NumpyTestCase):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        base_faces = PyramMesh.StructuredPyramBaseFaces(
            self.listmesh['GridDimension'])
        self.pyram_nodes =  PyramMesh.StructuredPyramNodes(
            self.listmesh['GridDimension'], base_faces)
        self.pyram_elements = PyramMesh.ProxyElements(
            PyramMesh.Element({'nodeCoords': self.testMesh.nodes.copy(),
                               'nodes': self.testMesh.elementNodes.copy(),
                               }))
    

    def test_ElementEdgeNumbers(self):
        edge_nodes = self.testMesh.edgeNodes
        edge_map = dict((tuple(e_n), i) for i, e_n in enumerate(edge_nodes))
        assert_equal(PyramMesh.numberElementEdges(self.pyram_elements, edge_map),
                     self.TestMesh.elementEdges)

    def test_ElementApexfaceNumbers(self):
        apex_facenodes = self.testMesh.apexfaceNodes
        apexface_map = dict((tuple(f_n), i)
                            for i, f_n in enumerate(apex_facenodes))
        assert_equal(
            PyramMesh.numberElementApexfaces(self.pyram_elements, apexface_map),
            self.testMesh.elementApexfaces)

class test_StructuredPyramMesh(test_GridGeneration,
                               _StructuredPyramEntityNodeGeneration):
    TestMesh = SixPyram
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = PyramMesh.Mesh(self.testMesh.listmesh)
    
    def test_ElementEdges(self):
        assert_equal(self.mesh.elements[:].edgenos, self.testMesh.elementEdges)
    
    def test_ElementNodes(self):
        elms = self.mesh.elements
        assert_array_equal([elm.nodes for elm in elms],
                           self.testMesh.elementNodes)
        for elm in elms:
            assert_array_equal(elm.nodeCoords,
                               [self.testMesh.nodes[i] for i in elm.nodes])

    def test_ElementNodeCoords(self):
        assert_array_equal(self.mesh.elements[:].nodeCoords,
                           self.testMesh.elementNodeCoords)

    def test_ElementApexfaces(self):
        assert_equal(self.mesh.elements[:].apexfacenos,
                     self.testMesh.elementApexfaces)

    def test_ElementBasefaces(self):
        assert_equal(self.mesh.elements[:].basefacenos,
                     self.testMesh.elementBasefaces)

    def test_baseface_nodes(self):
        assert_array_equal(self.mesh.basefaceNodes, self.testMesh.basefaceNodes)


class test_StructuredPyramMesh_twelve_pyram(test_StructuredPyramMesh):
    TestMesh = TwelvePyram
