from __future__ import division

import numpy as N
import scipy
from NewCode.Utilities import Struct, ProxyArray
from NewCode import ProxyList, Coordinates
from NewCode import Mesh as TetMesh

MeshItemClassFactory = TetMesh.MeshItemClassFactory

class ProxyElement(MeshItemClassFactory(
    ('facenos', 'edgenos', 'nodes', 'nodeCoords',
     'connect2elem', 'connect2face'), 'ProxyElement',
    computed_attrs=('nodeCoords', 'gridStepSize'))):
    @property
    def gridStepSize(self):
        return self._gridStepSize

class Mesh(object):
    localCoordLen = 6
    def __init__(self, listmesh):
        lm = listmesh
        self.nodes = CartesianGridNodes(lm['GridDimension'],
                                        lm['GridStepSize'],
                                        offset=lm['GridOffset'])        
        self.gridStepSize = lm['GridStepSize']
        self.offset = lm['GridOffset']
        self.elementNodes = StructuredBrickNodes(lm['GridDimension'])
        elementConnect2Elem, elementConnect2Face = calcElementConnectivity(
            self.elementNodes)
        self.elements = ProxyElements(
            Element({'nodes': self.elementNodes,
                     'nodeCoords': self.nodes,
                     'connect2face' : elementConnect2Face,
                     'connect2elem' : elementConnect2Elem,
                     'gridStepSize' : self.gridStepSize}))
        self.edgeNodes = StructuredBrickEdges(lm['GridDimension'])
        self.faceNodes = StructuredBrickFaces(lm['GridDimension'])
        edge_map = dict((tuple(e_n), i)
                        for i, e_n in enumerate(self.edgeNodes))
        face_map = dict((tuple(f_n), i)
                        for i, f_n in enumerate(self.faceNodes))
        self.elements.entity._edgenos = numberElementEdges(self.elements, edge_map)
        self.elements.entity._facenos = numberElementFaces(self.elements, face_map)
        faceConnect2Elem = calcFaceConnectivity(self.elements, len(self.faceNodes))
        self.faces = ProxyFaces(
            Face({'nodes': self.faceNodes,
                  'nodeCoords': self.nodes,
                  'connect2elem': faceConnect2Elem}))
        self.faces.entity._edgenos = numberFaceEdges(self.faces, edge_map)
        edgeConnect2Elem = calcEdgeConnectivity(self.elements, len(self.edgeNodes))
        self.edges = ProxyEdges(
            Edge({'nodes': self.edgeNodes,
                  'nodeCoords': self.nodes,
                  'connect2elem': edgeConnect2Elem},
                 self.faces))
        self.noFaces = len(self.faces)
        self.noEdges = len(self.edges)
        self.noNodes = len(self.nodes)
        
    def findClosestNode(self, r):
        dr = self.gridStepSize
        ijk = N.round((r-self.offset)/dr).astype(N.int32)
        return self.nodes.nodeNo(*ijk)
    
    def locatePoint(self, r):
        dr = self.gridStepSize
        ijk = N.round((r-self.offset)/dr-[1/2,1/2,1/2]).astype(N.int32)
        elno = self.elements.elementNo(*ijk)
        return elno, self.elements[elno].global2local(r)

class CartesianGridNodes(ProxyArray):
    def __init__(self, gridDimension, gridStepSize, offset=None, dtype=None):
        subarr = self.makeGridPoints(gridDimension, gridStepSize, dtype=dtype)
        self.dtype = dtype
        if offset is not None: subarr += offset
        self.gridDimension = gridDimension
        self.gridStepSize = gridStepSize
        self.offset = offset if offset is not None else N.array([0,0,0.])
        self.subarr = subarr
        self.nodeNo = node_no_gen(*gridDimension)
        self.nodeIjk = ijk_gen(*gridDimension)
        
    @staticmethod
    def makeGridPoints(gridDimension, gridStepSize, dtype):
        spaceDim=3
        nx,ny,nz = gridDimension
        dx, dy, dz = gridStepSize
        arr = N.empty((nx*ny*nz, spaceDim), dtype=dtype)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    arr[ny*nz*i + nz*j + k] = [i*dx, j*dy, k*dz]
        return arr

def node_no_gen(nnx, nny, nnz, offset=0):
    return lambda i,j,k : nny*nnz*i + nnz*j + k + offset

def ijk_gen(nx,ny,nz):
    nynz = ny*nz
    def nodeno_to_ijk(node_no):
        i,r = divmod(node_no, nynz)
        j,k = divmod(r, nz)
        return i,j,k
    return nodeno_to_ijk

class StructuredBrickNodes(ProxyArray):
    def __init__(self, gridDimension, dtype=N.int32):
        nx,ny,nz = gridDimension - 1    # Number of bricks in each dimension
        nnx, nny, nnz = gridDimension   # Number of nodes 
        nodesPerBrick = 8
        arr = N.empty((nx*ny*nz, nodesPerBrick), dtype=dtype)
        node_no = node_no_gen(nnx, nny, nnz)
        element_no = node_no_gen(nx, ny, nz)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    node_nos = [node_no(i+io,j+jo,k+ko)
                        for io in range(2) for jo in range(2) for ko in range(2)]
                    arr[element_no(i,j,k)] = node_nos
        self.subarr = arr
        self.gridDimension = gridDimension
        self.brickGridDimension = gridDimension - 1
        self.brickNo = node_no_gen(*self.brickGridDimension)

class StructuredBrickEdges(ProxyArray):
    def __init__(self, gridDimension, dtype=N.int32):
        nnx, nny, nnz = gridDimension   # Number of nodes in each dimension
        x_nx, x_ny, x_nz = (nnx-1),nny,nnz
        y_nx, y_ny, y_nz = nnx,(nny-1),nnz
        z_nx, z_ny, z_nz = nnx,nny,(nnz-1)
        x_n = x_nx*x_ny*x_nz            # Number of x-directed edges
        y_n = y_nx*y_ny*y_nz            # Number of y-directed edges
        z_n = z_nx*z_ny*z_nz            # Number of z-directed edges
        nodesPerEdge = 2
        arr = N.empty((x_n+y_n+z_n, nodesPerEdge), dtype=dtype)
        node_no = node_no_gen(nnx, nny, nnz)
        x_edge_no = node_no_gen(x_nx, x_ny, x_nz)
        y_edge_no = node_no_gen(y_nx, y_ny, y_nz)
        z_edge_no = node_no_gen(z_nx, z_ny, z_nz)
        for i in range(x_nx):
            for j in range(x_ny):
                for k in range(x_nz):
                    arr[x_edge_no(i,j,k)] = (
                        node_no(i,j,k), node_no(i+1,j,k))
        for i in range(y_nx):
            for j in range(y_ny):
                for k in range(y_nz):
                    arr[y_edge_no(i,j,k)+x_n] = (
                        node_no(i,j,k), node_no(i,j+1,k))
        for i in range(z_nx):
            for j in range(z_ny):
                for k in range(z_nz):
                    arr[z_edge_no(i,j,k)+x_n+y_n] = (
                        node_no(i,j,k), node_no(i,j,k+1))
        self.subarr = arr
        self.gridDimension = gridDimension
        self.edgeGridDimension = (x_n, y_n, z_n)

class StructuredBrickFaces(ProxyArray):
    def __init__(self, gridDimension, dtype=N.int32):
        nnx, nny, nnz = gridDimension   # Number of nodes in each dimension
        x_nx, x_ny, x_nz = nnx,(nny-1),(nnz-1)
        y_nx, y_ny, y_nz = (nnx-1),nny,(nnz-1)
        z_nx, z_ny, z_nz = (nnx-1),(nny-1),nnz
        x_n = x_nx*x_ny*x_nz            # Number of x-normal faces
        y_n = y_nx*y_ny*y_nz            # Number of y-normal faces
        z_n = z_nx*z_ny*z_nz            # Number of z-normal faces
        nodesPerFace = 4
        arr = N.empty((x_n+y_n+z_n, nodesPerFace), dtype=dtype)
        node_no = node_no_gen(nnx, nny, nnz)
        self.x_face_no = x_face_no = node_no_gen(x_nx, x_ny, x_nz)
        self.y_face_no = y_face_no = node_no_gen(y_nx, y_ny, y_nz, offset=x_n)
        self.z_face_no = z_face_no = node_no_gen(z_nx, z_ny, z_nz, offset=x_n+y_n)
        for i in range(x_nx):
            for j in range(x_ny):
                for k in range(x_nz):
                    arr[x_face_no(i,j,k)] = (
                        node_no(i,j,k), node_no(i,j,k+1),
                        node_no(i,j+1,k+1),node_no(i,j+1,k))
        for i in range(y_nx):
            for j in range(y_ny):
                for k in range(y_nz):
                    arr[y_face_no(i,j,k)] = (
                        node_no(i,j,k), node_no(i,j,k+1),
                        node_no(i+1,j,k+1),node_no(i+1,j,k))
        for i in range(z_nx):
            for j in range(z_ny):
                for k in range(z_nz):
                    arr[z_face_no(i,j,k)] = (
                        node_no(i,j,k), node_no(i,j+1,k),
                        node_no(i+1,j+1,k),node_no(i+1,j,k))

        self.subarr = arr
        self.gridDimension = gridDimension
        self.faceGridDimension = (x_n, y_n, z_n)

class Element(ProxyElement, Coordinates.BrickCoord):
    meshtype = 'hex'
    # Maps an element local edge number to element local node numbers
    LOCAL_EDGENODES = N.array([[1,5],   # 1 x
                               [2,6],   # 2 x
                               [3,7],   # 3 x
                               [4,8],   # 4 x
                               [1,3],   # 5 y
                               [2,4],   # 6 y
                               [5,7],   # 7 y
                               [6,8],   # 8 y
                               [1,2],   # 9 z
                               [3,4],   # 10 z
                               [5,6],   # 11 z
                               [7,8]], N.int32) - 1

    # Maps an element local face number to element local node numbers
    # These faces are numbered in accordance to Graglia97's convention
    LOCAL_FACENODES = N.array([[1,2,4,3], # 1 
                               [1,2,6,5], # 2 
                               [1,3,7,5], # 3 
                               [5,6,8,7], # 4 
                               [3,4,8,7], # 5 
                               [2,4,8,6]], N.int32) - 1
    # Which local faces is each edge connected to?
    LOCAL_EDGEFACES = N.array([[2,3],   # 1
                               [2,6],   # 2
                               [3,5],   # 3
                               [5,6],   # 4
                               [1,3],   # 5
                               [1,6],   # 6
                               [3,4],   # 7
                               [4,6],   # 8
                               [1,2],   # 9
                               [1,5],   # 10
                               [2,4],   # 11
                               [4,5]], N.int32) - 1
    # Maps an element local face number to element local edge number. Note that
    # the edge direction sense of the last two face edges are negative
    LOCAL_FACEEDGES = N.array([[9,6,10,5],  # 1 
                               [9,2,11,1],  # 2 
                               [5,3,7,1],   # 3 
                               [11,8,12,7], # 4 
                               [10,4,12,3], # 5 
                               [6,4,8,2]], N.int32) - 1
    # Using graglia97's dual numbering scheme where everything is defined
    # i.t.o.  face numbers rather than node numbers, nodes are numbered by the
    # three faces that intersect there
    FACEDUAL_NODE_NUMBERING = N.array([[1,2,3], # 1
                                       [1,2,6], # 2
                                       [1,3,5], # 3
                                       [1,5,6], # 4
                                       [2,3,4], # 5
                                       [2,4,6], # 6
                                       [3,4,5], # 7
                                       [4,5,6]], N.int32) - 1
    # Similarly, edges are numbered by the two faces that intersect there
    FACEDUAL_EDGE_NUMBERING = LOCAL_EDGEFACES
    FACEDUAL_FACE_NUMBERING = N.arange(6, dtype=N.int32)
    COUNT=Struct(**{'node':8, 'edge': 12, 'face':6, 'vol':1})

    size = property(lambda self: N.multiply.reduce(self.gridStepSize))

    def normals(self):
        return N.array([[1,0,0], [0,1,0], [0,0,1],
                        [1,0,0], [0,1,0], [0,0,1]], N.float64)

    def outwardNormals(self):
        return N.array([[-1,0,0], [0,-1,0], [0,0,-1],
                        [1,0,0], [0,1,0], [0,0,1]], N.float64)

def make_boundary_edge_set(edge_nodes):
    nnx, nny, nnz = edge_nodes.gridDimension   # Number of nodes in each dimension

class NumberEntities(object):
    entity_name = 'entity'
    def __call__(self, elements, entity_map):
        n_el, n_ent = len(elements), elements.entity.COUNT[self.entity_name]
        element_entities = N.zeros((n_el, n_ent), dtype=N.int32)
        local_entitynodes = getattr(
            elements[:], 'LOCAL_' + self.entity_name.upper() + 'NODES')
        for i, brick in enumerate(elements):
            element_entity_nodes = brick.nodes[local_entitynodes]
            element_entities[i] = [entity_map[tuple(en)]
                                   for en in element_entity_nodes]
        return element_entities

class NumberElementEdges(NumberEntities):
    entity_name = 'edge'
numberElementEdges = NumberElementEdges()

class NumberElementFaces(NumberEntities):
    entity_name = 'face'
numberElementFaces = NumberElementFaces()

class NumberFaceEdges(NumberElementEdges): pass
numberFaceEdges = NumberFaceEdges()

def calcElementConnectivity(brick_nodes):
    n_e = len(brick_nodes)
    n_f = 6                             # Number of faces per brick
    bgd = brick_nodes.brickGridDimension
    ijk_el = ijk_gen(*bgd)
    el_no = node_no_gen(*bgd)
    connect2elem = N.empty((n_e, n_f), dtype=N.int32)
    calc_elcon = lambda ijk: el_no(*ijk) if N.all(ijk > -1) and N.all(ijk < bgd) \
                 else -1
    el_ijks = N.array([ijk_el(el_i) for el_i in range(n_e)], N.int32)
    ijk_offsets = ((-1,0,0), (0,-1,0), (0,0,-1), (1,0,0), (0,1,0), (0,0,1))
    connect2elem = N.array([[calc_elcon(ijk + offset) for offset in ijk_offsets]
                            for ijk in el_ijks], N.int32)
    connect2face = N.where(connect2elem != -1, [3,4,5,0,1,2], -1)
    return connect2elem, connect2face
    

def calcFaceConnectivity(brick_elements, noFaces):
    connect2elem = N.zeros((noFaces, 2), N.int32) -1 # max 2 el connections/face
    for el_i, el in enumerate(brick_elements):
        for face_i, el_j in zip(el.facenos, el.connect2elem):
            if el_i > el_j and el_j >= 0: continue
            connect2elem[face_i] = el_i, el_j

    assert(N.all(connect2elem[:, 0] != -1))
    return connect2elem

def calcEdgeConnectivity(brick_elements, noEdges):
    connect2elem = (N.zeros((noEdges, 4), N.int32) -1).tolist()  # max 4 el connections/edge
    for el_i, el in enumerate(brick_elements):
        for edge_i in el.edgenos:
            e_con2 = connect2elem[edge_i]
            e_con2[e_con2.index(-1)] = el_i
    return N.array(connect2elem, N.int32)
    
class ProxyElements(ProxyList.ProxyList):
    def __init__(self, *names, **kwargs):
        super(ProxyElements, self).__init__(*names, **kwargs)
        try: self.elementNo = self.entity._nodes.brickNo
        except AttributeError: pass

    def list_repr(self, *names, **kwargs):
        lr = super(ProxyElements, self).list_repr(*names, **kwargs)
        lr['gridStepSize']=self.entity.gridStepSize
        return lr


class ProxyFaces(TetMesh.ProxyFaces): pass

class ProxyEdges(TetMesh.ProxyEdges): pass
    
class FaceArea(object):
    def area(self):
        a,b,c,d = self.nodeCoords
        ab = b - a
        ad = d - a
        return N.linalg.norm(scipy.cross(ab,ad))


class Edge(TetMesh.Edge):
    meshtype = 'hex'
    computedAttributes = TetMesh.Edge.computedAttributes + ('gridStepSize',)
    
class Face(TetMesh.BoundaryFace, FaceArea):
    meshtype = 'hex'
    LOCAL_EDGENODES = N.array([[0,1], [1,2], [3,2], [0,3]], N.int32)
    COUNT = Struct(**{'node':4, 'edge':4, 'face':1})

