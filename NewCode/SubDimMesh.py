from __future__ import division

import numpy as N
from itertools import izip
from NewCode.Utilities import Struct
from NewCode import Mesh, ProxyList

ProxySurfaceElement = Mesh.MeshItemClassFactory(
    ('superEdgenos', 'edgenos', 'nodes', 'nodeCoords', 'connect2SuperElems', 'superFacenos',
     'superLocalFacenos',), 'ProxySurfaceElement')

class FaceElement(ProxySurfaceElement, Mesh.FaceArea):
    LOCAL_EDGENODES = N.array([[1,2],
                               [1,3],
                               [2,3],
                               [2,4]], N.int32) - 1
    COUNT=Struct(**{'node':3, 'edge': 3, 'face':1})

    size = property(Mesh.FaceArea.area)
    superNo = ProxySurfaceElement.superFacenos

ProxyEdge = Mesh.MeshItemClassFactory(('nodes', 'nodeCoords', 'superNo'), 'ProxyEdge')

class Edge(ProxyEdge):
    pass


class SubSurface(object):
    MeshModule = Mesh
    FaceElement = FaceElement
    Edge = Edge
    def __init__(self, mesh, faceSelector):
        self.superMesh = mesh
        self.faceSelector = faceSelector
        self.superFacenos = N.fromiter((i for i,f in enumerate(mesh.faces)
                                        if faceSelector(f)), N.int32)
        self.superEdges = mesh.edges
        self.superNodes = mesh.nodes
        # At some stage this will have to be the subset of nodes actually
        # included for if/when nodal basis functions are supported
        self.nodes = mesh.nodes         
        connect2SuperElems = mesh.faces[self.superFacenos].connect2elem
        els = mesh.elements
        def getSuperLocalFaceNo(f_els, f_i):
            return [els[elno].facenos.tolist().index(f_i) if elno > -1 \
                    else -1 for elno in f_els]
        superLocalFacenos = N.array([getSuperLocalFaceNo(f_els,f_i)
                                     for f_els,f_i in zip(
            connect2SuperElems, self.superFacenos)], N.int32)
        superEdgenos = mesh.faces[self.superFacenos].edgenos
        uniqueSuperEdgenos = N.unique(superEdgenos)
        usedSuperEdges = self.superEdges[uniqueSuperEdgenos]
        self.edges = ProxyEdges(self.Edge({
            'nodes': usedSuperEdges.nodes,
            'nodeCoords': self.superMesh.nodes,
            'superNo': uniqueSuperEdgenos,}))
        superToSubEdgenos = dict(izip(uniqueSuperEdgenos, xrange(len(uniqueSuperEdgenos))))
        # Subdim local edgenos are made by looking up the index number of each
        # element's superEdgenos in uniqueSuperEdgenos
        edgenos = N.array([[superToSubEdgenos[i] for i in se] for se in superEdgenos],
                          N.int32)
        self.elements = self.MeshModule.ProxyElements(self.FaceElement({
            'superEdgenos' : superEdgenos,
            'edgenos' : edgenos,
            'nodes' : mesh.faces[self.superFacenos].nodes,
            'nodeCoords': mesh.nodes,
            'connect2SuperElems': connect2SuperElems,
            'superFacenos': self.superFacenos,
            'superLocalFacenos': superLocalFacenos}))
        

class ProxyEdges(Mesh.ProxyEntities): pass
