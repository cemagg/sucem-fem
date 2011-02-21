from __future__ import division

import numpy as N
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
from NewCode.tests.TestMeshes import InscribedTetMesh
from NewCode import Mesh, SubDimMesh
from NewCode.DifferentialForm import Discretiser

class test_SubMesh(TestCase):
    def setUp(self):
        self.mesh = Mesh.Mesh(InscribedTetMesh.listmesh)

    def test_init(self):
        fs = lambda face: face.onBoundary
        surfaceSubMesh = SubDimMesh.SubSurface(self.mesh, faceSelector=fs)
        # Number of boundary faces on InscribedTetMesh
        assert_equal(len(surfaceSubMesh.elements), 12) 
        boundaryFaces = [i for i,bf in enumerate(InscribedTetMesh.BoundaryFaces)
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
        fs = lambda face: 10 in face.connect2elem # Faces of the center-tet
        surfaceSubMesh = SubDimMesh.SubSurface(self.mesh, faceSelector=fs)
        assert_equal(surfaceSubMesh.superFacenos, [24, 25, 26, 27])
        assert_equal(set(surfaceSubMesh.elements[:].superEdgenos.flat),
                     set([5, 10, 14, 18, 21, 23]))
        assert_equal(surfaceSubMesh.edges[:].superNo, [5, 10, 14, 18, 21, 23])
        assert_equal(surfaceSubMesh.edges[:].nodeCoords,
                     self.mesh.edges[[5, 10, 14, 18, 21, 23]].nodeCoords)
        assert_equal(surfaceSubMesh.elements[:].edgenos,
                     [[0,2,5], [3,0,1], [1,2,4], [3,5,4]])
        
    def test_super(self):
        # Select only faces made up of nodes 1,2,3,5, all on one side of larger tet
        fs = lambda face: face.onBoundary and \
             len(set([1,2,3,5]) - set(face.nodes+1)) == 1
        surfaceSubMesh = SubDimMesh.SubSurface(self.mesh, faceSelector=fs)

        assert_equal(surfaceSubMesh.elements[:].superFacenos+1, [5,13,17])
        assert_equal(surfaceSubMesh.elements[:].connect2SuperElems, [[1,-1],
                                                                     [3,-1],
                                                                     [4,-1]])
        assert_equal(surfaceSubMesh.elements[:].superLocalFacenos,
                     [[0, -1], [0, -1], [0, -1]])
