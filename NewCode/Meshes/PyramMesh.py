from __future__ import division
import pdb
import numpy as N
from NewCode.Utilities import Struct, ProxyArray
from NewCode import Coordinates, ProxyList
from NewCode.Meshes import BrickMesh as BM
from NewCode.Meshes import TetMesh as TM


pyr_per_brick = 6; nodesPerEdge = 2 ; nodes_per_pyr = 5

class ProxyElement(BM.MeshItemClassFactory(
    ('apexfacenos', 'basefacenos', 'edgenos', 'nodes', 'nodeCoords',
     'connect2elem', 'connect2face'), 'ProxyElement',
    computed_attrs=('nodeCoords',))): pass

def calc_no_bricks(gridDimension):
    nbx, nby, nbz = gridDimension - 1
    return nbx*nby*nbz, (nbx, nby, nbz)

class Mesh(object):
    localCoordLen = 5
    def __init__(self, listmesh):
        lm = listmesh
        self.nodes = PyramGridNodes(lm['GridDimension'],
                                    lm['GridStepSize'],
                                    offset=lm['GridOffset'])        
        self.gridStepSize = lm['GridStepSize']
        self.offset = lm['GridOffset']
        self.basefaceNodes = StructuredPyramBaseFaces(lm['GridDimension'])
        self.basefaces = ProxyList.ProxyList(BaseFace(dict(
            nodes=self.basefaceNodes, nodeCoords=self.nodes)))
        self.elementNodes = StructuredPyramNodes(lm['GridDimension'],
                                                 self.basefaceNodes)
        self.elements = ProxyElements(
            Element({'nodes': self.elementNodes,
                     'nodeCoords': self.nodes,}))
        self.edgeNodes = StructuredPyramEdges(lm['GridDimension'])
        edge_map = dict((tuple(e_n), i)
                        for i, e_n in enumerate(self.edgeNodes))
        self.elements.entity._edgenos = numberElementEdges(self.elements, edge_map)
        self.edges = ProxyEdges(
            Edge({'nodes': self.edgeNodes,
                  'nodeCoords': self.nodes}))
        brick_els = ProxyList.ProxyList(BM.Element(dict(
            nodes=BM.StructuredBrickNodes(lm['GridDimension']),
            nodeCoords=self.nodes, gridStepSize=lm['GridStepSize'])))
        brick_els.entity._edgenos = BM.numberElementEdges(brick_els, edge_map)
        self.apexfaceNodes = StructuredApexfaceNodes(
            brick_els, self.edges, self.nodes.pyram_offset)
        self.apexfaces = ProxyList.ProxyList(ApexFace(dict(
            nodes=self.apexfaceNodes, nodeCoords=self.nodes)))
        apexface_map = dict((tuple(f_n), i)
                            for i, f_n in enumerate(self.apexfaceNodes))
        self.elements.entity._apexfacenos = numberElementApexfaces(
            self.elements, apexface_map)
        baseface_map = dict((tuple(f_n), i)
                            for i, f_n in enumerate(self.basefaceNodes))
        self.elements.entity._basefacenos = numberElementBasefaces(
            self.elements, baseface_map)
        self.noEdges = len(self.edges)
        self.noNodes = len(self.nodes)
        self.noApexfaces = len(self.apexfaces)
        self.noBasefaces = len(self.basefaces)
        
class PyramGridNodes(BM.CartesianGridNodes):
    def __init__(self, gridDimension, gridStepSize, offset=None, dtype=None):
        BM.CartesianGridNodes.__init__(self, gridDimension, gridStepSize,
                                       offset=offset, dtype=dtype)
        pyr_phys_offs = gridStepSize/2.
        if offset is not None: pyr_phys_offs += offset
        pyr_nodes = self.makeGridPoints(gridDimension-1, gridStepSize, dtype=dtype)
        pyr_nodes += pyr_phys_offs
        self.pyram_offset = len(self.subarr)
        self.subarr = N.vstack((self.subarr, pyr_nodes))

class StructuredPyramNodes(ProxyArray):
    def __init__(self, gridDimension, base_facenodes, dtype=N.int32):
        no_bricks, (nbx,nby,nbz) = calc_no_bricks(gridDimension)
        nnx, nny, nnz = gridDimension   # Number of nodes 
        pyr_node_offset = nnx*nny*nnz
        pyr_nodes = N.empty((pyr_per_brick*no_bricks, nodes_per_pyr), dtype=dtype)
        brick_no = BM.node_no_gen(nbx, nby, nbz)
        node_no = BM.node_no_gen(nnx, nny, nnz)
        for i in range(nbx):
            for j in range(nby):
                for k in range(nbz):
                    apex_node = pyr_node_offset + brick_no(i,j,k)
                    el_offs = brick_no(i,j,k)*pyr_per_brick
                    b_fs = [base_facenodes.x_face_no(i,j,k),
                            base_facenodes.x_face_no(i+1,j,k),
                            base_facenodes.y_face_no(i,j,k),
                            base_facenodes.y_face_no(i,j+1,k),
                            base_facenodes.z_face_no(i,j,k),
                            base_facenodes.z_face_no(i,j,k+1)]
                    el_base_nodes = N.sort(base_facenodes[b_fs])
                    pyr_nodes[el_offs:el_offs+pyr_per_brick, 0:4] = el_base_nodes
                    pyr_nodes[el_offs:el_offs+pyr_per_brick, 4] = apex_node
        self.subarr = pyr_nodes
        
class StructuredPyramEdges(BM.StructuredBrickEdges):
    def __init__(self, gridDimension, dtype=N.int32):
        BM.StructuredBrickEdges.__init__(self, gridDimension, dtype=dtype)
        self.pyram_offset = len(self.subarr)
        no_bricks, (nbx, nby, nbz) = calc_no_bricks(gridDimension)
        nnx, nny, nnz = gridDimension   # Number of nodes 
        pyr_node_offset = nnx*nny*nnz
        no_pyr = no_bricks*pyr_per_brick
        extra_edges_per_brick = 8 ;
        no_pyr_internal_edges = extra_edges_per_brick*no_bricks
        edge_node_arr = N.empty((no_pyr_internal_edges, nodesPerEdge), dtype=dtype)
        brick_no = BM.node_no_gen(nbx, nby, nbz)
        node_no = BM.node_no_gen(*gridDimension)
        for i in range(nbx):
            for j in range(nby):
                for k in range(nbz):
                    apex_node = pyr_node_offset + brick_no(i,j,k)
                    edge_offs = brick_no(i,j,k)*extra_edges_per_brick
                    edge_nodes = N.array([(node_no(i+ii, j+jj, k+kk), apex_node)
                                          for ii in range(nodesPerEdge)
                                          for jj in range(nodesPerEdge)
                                          for kk in range(nodesPerEdge)], N.int32)
                    edge_node_arr[
                        edge_offs:edge_offs+extra_edges_per_brick] = edge_nodes

        self.subarr = N.vstack((self.subarr, edge_node_arr))
                    
        
class StructuredPyramBaseFaces(BM.StructuredBrickFaces): pass

def StructuredApexfaceNodes(brick_els, brick_edges, apex_node_offset):
    face_nodes = []
    for bel in brick_els:
        apex_node = bel.index + apex_node_offset
        for edge in  brick_edges[bel.edgenos]:
            face_nodes.append(edge.nodes.tolist() + [apex_node])
    return face_nodes
            
numberElementEdges = BM.numberElementEdges

class Element(ProxyElement, Coordinates.PyramCoord):
    COUNT=Struct(node=5, edge=8, face=5, apexface=4, baseface=1, vol=1)
    # Maps an element local edge number to element local node numbers
    LOCAL_EDGENODES = N.array([[1,2],
                               [1,3],
                               [2,4],
                               [3,4],
                               [1,5],
                               [2,5],
                               [3,5],
                               [4,5]], N.int32) - 1
    # Which local faces is each edge connected to? See graglia99 fig 2a. for
    # the edge direction convention that this one is relative to.
    # loc edgno, edge sense sign relative to graglia99
    LOCAL_EDGEFACES = N.array([[1,5],   # 1, -
                               [2,5],   # 2, +
                               [4,5],   # 3, -
                               [3,5],   # 4, +
                               [1,2],   # 5, +
                               [4,1],   # 6, +
                               [2,3],   # 7, +
                               [3,4]],  # 8, +
                              N.int32) - 1
    LOCAL_APEXFACENODES = N.array([[1,2,5], [1,3,5], [3,4,5], [2,4,5]],
                                  N.int32) - 1
    LOCAL_APEXFACEEDGES = N.array([[1,5,6],
                                   [2,5,7],
                                   [4,7,8],
                                   [3,6,8]], N.int32) - 1
    LOCAL_APEXFACE_OTHERBASENODES = N.array([[3,4], [2,4], [1,2], [1,3]],
                                            N.int32) - 1
    LOCAL_BASEFACENODES = N.array([[1,2,4,3]], N.int32) - 1
    LOCAL_BASEFACENODES_DIAG_OPPOSED = N.array([4,3,2,1], N.int32) - 1
    LOCAL_BASEFACEEDGES = N.array([1,3,4,2], N.int32) - 1
    FACEDUAL_EDGE_NUMBERING = LOCAL_EDGEFACES
    FACEDUAL_EDGE_SENSE = N.array([-1,1,-1,1,1,1,1,1], N.int32)
    FACEDUAL_FACE_NUMBERING = N.arange(COUNT.face, dtype=N.int32)    

    @property
    def volume(self):
        nc = self.nodeCoords
        l1, l2, l3 = l123 = nc[[1,2,4]] - nc[0]
        return N.abs(N.dot(l1, N.cross(l2,l3)))/3

    size = volume

    
class ProxyElements(ProxyList.ProxyList): pass
class ProxyEdges(TM.ProxyEdges): pass
class ApexFace(TM.ProxyFace, TM.FaceArea): pass
class BaseFace(TM.ProxyFace, BM.FaceArea): pass
class Edge(TM.ProxyEdge): pass

class NumberElementApexfaces(BM.NumberEntities):
    entity_name = 'apexface'
numberElementApexfaces = NumberElementApexfaces()
class NumberElementBasefaces(BM.NumberEntities):
    entity_name = 'baseface'
numberElementBasefaces = NumberElementBasefaces()
