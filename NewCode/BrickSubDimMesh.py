from __future__ import division

import numpy as N
from itertools import izip
from NewCode.Utilities import Struct
from NewCode import ProxyList
from NewCode.Meshes import BrickMesh
from NewCode import SubDimMesh

ProxySurfaceElement = BrickMesh.MeshItemClassFactory(
    ('superEdgenos', 'edgenos', 'nodes', 'nodeCoords', 'connect2SuperElems', 'superFacenos',
     'superLocalFacenos',), 'ProxySurfaceElement')

class FaceElement(ProxySurfaceElement, BrickMesh.FaceArea):
    LOCAL_EDGENODES = BrickMesh.Face.LOCAL_EDGENODES
    COUNT = BrickMesh.Face.COUNT

    @property
    def size(self): return self.area()

    # This should really be selected from the parent element depending on the
    # face orientation
    gridStepSize = None                     
    superNo = ProxySurfaceElement.superFacenos
    
ProxyEdge = BrickMesh.MeshItemClassFactory(('nodes', 'nodeCoords', 'superNo'), 'ProxyEdge')

class Edge(ProxyEdge):
    pass


class SubSurface(SubDimMesh.SubSurface):
    MeshModule = BrickMesh
    FaceElement = FaceElement
    Edge = Edge
    
