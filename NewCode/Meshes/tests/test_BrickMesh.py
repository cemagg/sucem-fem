from __future__ import division

import numpy as N
from numpy.testing import NumpyTestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests import xfail
from NewCode.tests.BrickMeshes import OneBrick, TwoBricks, TwoBricksX, FourBricks
from NewCode import ProxyList
from NewCode.Meshes import BrickMesh

class test_BrickElementLocalNumbering(NumpyTestCase):
    def setUp(self):
        self.Element = BrickMesh.Element
    def test_DualLocalEdgeNodeNumbering(self):
        # Test that the edges defined by their nodes are consitent with the
        # edges as defined by the two faces they intersect with by using the
        # dual node numbering scheme
        Element = self.Element
        # Twelve edges intersecting 2 faces each
        assert_equal(Element.LOCAL_EDGEFACES.shape, (Element.COUNT.edge, 2))
        assert_equal(Element.FACEDUAL_EDGE_NUMBERING, Element.LOCAL_EDGEFACES)
        d_nodes = Element.FACEDUAL_NODE_NUMBERING
        for edge_nodes, dual_edge_nodes in zip(
            Element.LOCAL_EDGENODES, Element.FACEDUAL_EDGE_NUMBERING):
            d_nodeset = [set(d_nodes[i]) for i in edge_nodes]
            # Note that the edge facenos are sorted() in ascending order
            assert_equal(sorted(d_nodeset[0].intersection(d_nodeset[1])),
                         dual_edge_nodes)

    def test_DualLocalFaceNodeNumbering(self):
        # Test that the face node defn actually forms the face of a brick using
        # the dual node numbering scheme
        Element = self.Element
        d_nodes = Element.FACEDUAL_NODE_NUMBERING
        for faceno, face_nodes in enumerate(Element.LOCAL_FACENODES):
            d_nodeset = [set(d_nodes[i]) for i in face_nodes]
            assert_equal(set([faceno]), reduce(lambda a,b: a&b, d_nodeset))

    def test_LocalEdgeOrdering(self):
        # Test that the edges are defined in ascending order of node number
        for n in self.Element.LOCAL_EDGENODES: self.assert_(n[0] < n[1])
        # Twelve edges should be defined by 2 nodes each
        assert_equal(self.Element.LOCAL_EDGENODES.shape,
                     (self.Element.COUNT.edge,2))

    def test_LocalFacenodeOrdering(self):
        for fns in self.Element.LOCAL_FACENODES:
            assert_equal(N.argsort(fns), [0,1,3,2])

    def test_FaceEdges(self):
        Element = self.Element
        ellocal_edgenodes = Element.LOCAL_EDGENODES
        ellocal_facenodes = Element.LOCAL_FACENODES
        ellocal_faceedges = Element.LOCAL_FACEEDGES
        # Six faces defined by 4 nodes each
        assert_equal(ellocal_faceedges.shape, (Element.COUNT.face,4))
        # The face-local edges are cyclicaly defined using the face nodes. The
        # element local edges are always defined from the lowest to the highest
        # node number, hence this reordering is needed.
        facelocal_edgenodes = N.array([[1,2], [2,3], [4,3], [1,4]], N.int32) - 1
        for edgenos, facenodes in zip(ellocal_faceedges, ellocal_facenodes):
            edgenodes_via_facedges = ellocal_edgenodes[edgenos]
            edgenodes_via_facenodes = [facenodes[n] for n in facelocal_edgenodes]
            assert_equal(edgenodes_via_facedges, edgenodes_via_facenodes)

    def test_Count(self):
        c = self.Element.COUNT
        assert_equal([c.node, c.edge, c.face, c.vol],
                     [8, 12, 6, 1])


class test_BrickElements(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.elements = BrickMesh.ProxyElements(
            BrickMesh.Element({'nodeCoords': self.testMesh.nodes.copy(),
                               'nodes': self.testMesh.elementNodes.copy(),
                               'gridStepSize' : self.testMesh.listmesh['GridStepSize']}))

    def test_ElementNodes(self):
        elms = self.elements
        assert_array_equal([elm.nodes for elm in elms],
                           self.testMesh.elementNodes)
        for elm in elms:
            assert_array_equal(elm.nodeCoords,
                               [self.testMesh.nodes[i] for i in elm.nodes])

class test_BrickElementsX(test_BrickElements):
    TestMesh = TwoBricksX


class test_GridGeneration(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh

    def test_NodeCoordinates(self):
        nodes = self.nodes
        assert_equal(nodes, self.testMesh.nodes)

    def test_nodeIjk(self):
        nodes = self.nodes
        assert_equal([nodes.nodeIjk(i) for i in range(len(nodes))],
                      self.testMesh.nodesIjk)
        assert_equal([nodes.nodeNo(*nodes.nodeIjk(i)) for i in range(len(nodes))],
                     range(len(nodes)))
    @property
    def nodes(self):
        lm = self.listmesh
        return BrickMesh.CartesianGridNodes(lm['GridDimension'],
                                            lm['GridStepSize'],
                                            offset=lm['GridOffset'])

class test_GridGenerationX(test_GridGeneration):
    TestMesh = TwoBricksX

class _StructuredBrickEntityNodeGeneration(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh

    @property
    def brick_nodes(self):
        return BrickMesh.StructuredBrickNodes(self.listmesh['GridDimension'])
    def test_BrickNodeNumbers(self):
        brick_nodes = self.brick_nodes
        assert_equal(brick_nodes, self.testMesh.elementNodes)
        assert_equal(brick_nodes.gridDimension, self.testMesh.gridDimension)
        assert_equal(brick_nodes.brickGridDimension, self.testMesh.elementGridDimension)
        
    @property
    def edge_nodes(self):
        return BrickMesh.StructuredBrickEdges(self.listmesh['GridDimension'])
    def test_EdgeNodeNumbers(self):
        edge_nodes = self.edge_nodes
        assert_equal(edge_nodes, self.testMesh.edgeNodes)
        assert_equal(edge_nodes.edgeGridDimension, self.testMesh.edgeGridDim)
        
        
    @property
    def face_nodes(self):
        return BrickMesh.StructuredBrickFaces(self.listmesh['GridDimension'])
    def test_FaceNodeNumbers(self):
        face_nodes = self.face_nodes
        assert_equal(face_nodes, self.testMesh.faceNodes)
        assert_equal(face_nodes.faceGridDimension, self.testMesh.faceGridDim)

class test_StructuredBrickEntityNodeGeneration(_StructuredBrickEntityNodeGeneration):
    TestMesh = TwoBricksX

class test_BrickElementEntityNumbers(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        self.brick_nodes =  BrickMesh.StructuredBrickNodes(self.listmesh['GridDimension'])
        self.brick_elements = BrickMesh.ProxyElements(
            BrickMesh.Element({'nodeCoords': self.testMesh.nodes.copy(),
                               'nodes': self.testMesh.elementNodes.copy(),
                               'gridStepSize' : self.testMesh.listmesh['GridStepSize']}))
                               

    def test_ElementEdgeNumbers(self):
        edge_nodes = self.testMesh.edgeNodes
        edge_map = dict((tuple(e_n), i) for i, e_n in enumerate(edge_nodes))
        assert_equal(BrickMesh.numberElementEdges(self.brick_elements, edge_map),
                     self.TestMesh.elementEdges)

    def test_ElementFaceNumbers(self):
        face_nodes = self.testMesh.faceNodes
        face_map = dict((tuple(f_n), i) for i, f_n in enumerate(face_nodes))
        assert_equal(BrickMesh.numberElementFaces(self.brick_elements, face_map),
                     self.TestMesh.elementFaces)

    def test_FaceEdgeNumbers(self):
        edge_nodes = self.testMesh.edgeNodes
        edge_map = dict((tuple(e_n), i) for i, e_n in enumerate(edge_nodes))
        faces = ProxyList.ProxyList(
            BrickMesh.Face({'nodes': self.testMesh.faceNodes,
                            'nodeCoords': self.testMesh.nodes,
                            'connect2elem': self.testMesh.faceConnect2Elem}))
        face_edges = BrickMesh.numberFaceEdges(faces, edge_map)
        assert_equal(face_edges, self.testMesh.faceEdges)
        
class test_BrickElementEntityNumbersX(test_BrickElementEntityNumbers):
    TestMesh = TwoBricksX

class test_BrickElementConnections(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh
        self.brick_nodes =  BrickMesh.StructuredBrickNodes(self.listmesh['GridDimension'])
 
    def test_ElementConnect2Elem(self):
        brick_nodes = self.brick_nodes
        elementConnect2Elem, elementConnect2Face = BrickMesh.calcElementConnectivity(
            brick_nodes)
        assert_equal(elementConnect2Elem, self.testMesh.elementConnect2Elem)
        assert_equal(elementConnect2Face, self.testMesh.elementConnect2Face)
    

    @property
    def brick_elements(self):
        return BrickMesh.ProxyElements(BrickMesh.Element(
            {'nodeCoords': self.testMesh.nodes.copy(),
             'nodes': self.testMesh.elementNodes.copy(),
             'facenos': self.testMesh.elementFaces.copy(),
             'connect2elem' : self.testMesh.elementConnect2Elem.copy(),
             'connect2face' : self.testMesh.elementConnect2Face.copy(),
             'edgenos' : self.testMesh.elementEdges.copy(),
             'gridStepSize' : self.testMesh.listmesh['GridStepSize']}))

    def test_FaceConnect2Elem(self):
        brick_elements = self.brick_elements
        assert_equal(BrickMesh.calcFaceConnectivity(
            brick_elements, len(self.testMesh.faceNodes)),
                     self.testMesh.faceConnect2Elem)
        

    def test_EdgeConnect2Elem(self):
        brick_elements = self.brick_elements
        assert_equal(BrickMesh.calcEdgeConnectivity(
            brick_elements, len(self.testMesh.edgeNodes)),
                     self.testMesh.edgeConnect2Elem)

class test_BrickElementConnectionsX(test_BrickElementConnections):
    TestMesh = TwoBricksX

class test_BoundaryEntities(NumpyTestCase):
    TestMesh = TwoBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.listmesh = self.testMesh.listmesh

    @property
    def brick_elements(self):
        return BrickMesh.ProxyElements(BrickMesh.Element(
            {'nodeCoords': self.testMesh.nodes.copy(),
             'nodes': self.testMesh.elementNodes.copy(),
             'facenos': self.testMesh.elementFaces.copy(),
             'connect2elem' : self.testMesh.elementConnect2Elem.copy(),
             'connect2face' : self.testMesh.elementConnect2Face.copy(),
             }))
    
    def test_boundaryFaces(self):
        faces = BrickMesh.ProxyFaces(
            BrickMesh.Face({'nodes': self.testMesh.faceNodes,
                            'nodeCoords': self.testMesh.nodes,
                            'connect2elem': self.testMesh.faceConnect2Elem}))
        assert_equal(faces.entity._onBoundary, self.testMesh.boundaryFaceSet)

    def test_boundaryEdges(self):
        faces = BrickMesh.ProxyFaces(
            BrickMesh.Face({'nodes': self.testMesh.faceNodes,
                            'nodeCoords': self.testMesh.nodes,
                            'connect2elem': self.testMesh.faceConnect2Elem,
                            'edgenos': self.testMesh.faceEdges}))
        edges = BrickMesh.ProxyEdges(
            BrickMesh.Edge({'nodes': self.testMesh.edgeNodes,
                            'nodeCoords': self.testMesh.nodes}, faces))
        assert_equal(edges.entity._onBoundary, self.testMesh.boundaryEdgeSet)        

class test_BoundaryEntitiesX(test_BoundaryEntities): 
    TestMesh = TwoBricksX
class test_BoundaryEntities4(test_BoundaryEntities): 
    TestMesh = FourBricks

class test_StructuredBrickMesh(test_BrickElements, test_GridGeneration,
                               _StructuredBrickEntityNodeGeneration):
    TestMesh = FourBricks
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)
    
    def test_ElementEdges(self):
        assert_equal(self.mesh.elements[:].edgenos, self.testMesh.elementEdges)

    def test_ElementFaces(self):
        assert_equal(self.mesh.elements[:].facenos, self.testMesh.elementFaces)

    def test_FaceEdges(self):
        assert_equal(self.mesh.faces[:].edgenos, self.testMesh.faceEdges)

    def test_ElementConnect2(self):
        assert_equal(self.mesh.elements[:].connect2elem, self.testMesh.elementConnect2Elem)
        assert_equal(self.mesh.elements[:].connect2face, self.testMesh.elementConnect2Face)

    def test_FaceConnect2(self):
        assert_equal(self.mesh.faces[:].connect2elem, self.testMesh.faceConnect2Elem)

    def test_EdgeConnect2(self):
        assert_equal(self.mesh.edges[:].connect2elem, self.testMesh.edgeConnect2Elem)

    def test_BoundaryFaces(self):
        assert_equal(self.mesh.faces[:].onBoundary, self.testMesh.boundaryFaces)

    def test_BoundaryEdges(self):
        assert_equal(self.mesh.edges[:].onBoundary, self.testMesh.boundaryEdges)

    def test_findClosestNode(self):
        assert_equal(self.mesh.findClosestNode([1,2,3-0.5]), 13)

    def test_locatePoint(self):
        assert_equal(self.mesh.locatePoint([1-0.5,2-0.5,3-0.5]),
                     (0, [1/2, 3/4, 5/6, 1-1/2, 1-3/4, 1-5/6]))
        
        assert_equal(self.mesh.locatePoint([1-0.5,4-0.5,6-0.5]),
                     (3, [1/2, 3/4, 5/6, 1-1/2, 1-3/4, 1-5/6]))
    @property
    def elements(self): return self.mesh.elements
    @property
    def nodes(self): return self.mesh.nodes
    @property
    def brick_nodes(self): return self.mesh.elementNodes
    @property
    def edge_nodes(self): return self.mesh.edgeNodes
    @property
    def face_nodes(self): return self.mesh.faceNodes

class test_BrickFaceArea(NumpyTestCase):
    TestMesh = OneBrick
    def setUp(self):
        self.testMesh = self.TestMesh()
        self.mesh = BrickMesh.Mesh(self.testMesh.listmesh)

    def test_area(self):
        assert_equal([f.area()
                      for f in self.mesh.faces[self.mesh.elements[0].facenos]],
                     [6,3,2,6,3,2])
