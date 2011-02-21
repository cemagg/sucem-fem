from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.BrickMeshes import FourBricks
from NewCode import BrickSubDimMesh
from NewCode.Meshes import BrickMesh
from NewCode.DifferentialForm import Discretiser

class test_SubMesh(TestCase):
    def setUp(self):
        self.mesh = BrickMesh.Mesh(FourBricks.listmesh)

    def test_init(self):
        fs = lambda face: face.onBoundary
        surfaceSubMesh = BrickSubDimMesh.SubSurface(self.mesh, faceSelector=fs)
        # Number of boundary faces on InscribedTetMesh
        assert_equal(len(surfaceSubMesh.elements), 16) 
        boundaryFaces = [i for i,bf in enumerate(FourBricks.boundaryFaces)
                         if bf]
        assert_equal(surfaceSubMesh.superFacenos, boundaryFaces)
        boundaryFaceEls = self.mesh.faces[boundaryFaces]
        assert_equal(surfaceSubMesh.elements[:].nodeCoords, boundaryFaceEls.nodeCoords)
        assert_equal(surfaceSubMesh.elements[:].nodes, boundaryFaceEls.nodes)
        assert_equal(surfaceSubMesh.elements[:].superEdgenos, boundaryFaceEls.edgenos)
        self.assert_(surfaceSubMesh.nodes is self.mesh.nodes)
        self.assert_(surfaceSubMesh.superNodes is self.mesh.nodes)
        self.assert_(surfaceSubMesh.superEdges is self.mesh.edges)
        
    def test_Edges(self):
        fs = lambda face: face.index == 15 # Top face of el 0/bottom of el 1
        surfaceSubMesh = BrickSubDimMesh.SubSurface(self.mesh, faceSelector=fs)
        assert_equal(surfaceSubMesh.superFacenos, [15])
        assert_equal(surfaceSubMesh.elements[:].superEdgenos,
                     [[10,4,16,1]])
        assert_equal(surfaceSubMesh.edges[:].superNo, [1,4,10,16])
        assert_equal(surfaceSubMesh.edges[:].nodeCoords,
                     self.mesh.edges[[1,4,10,16]].nodeCoords)
        assert_equal(surfaceSubMesh.elements[:].edgenos,
                     [[2,1,3,0]])
        
    def test_super(self):
        fs = lambda face: face.index in [8,9,12,13] # Faces on left/right side
        surfaceSubMesh = BrickSubDimMesh.SubSurface(self.mesh, faceSelector=fs)

        assert_equal(surfaceSubMesh.elements[:].superFacenos, [8,9,12,13])
        assert_equal(surfaceSubMesh.elements[:].connect2SuperElems, [[0,-1],
                                                                     [1,-1],
                                                                     [2,-1],
                                                                     [3,-1]])
        assert_equal(surfaceSubMesh.elements[:].superLocalFacenos,
                     [[1, -1], [1, -1], [4, -1], [4,-1]])
