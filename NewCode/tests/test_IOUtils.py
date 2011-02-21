from __future__ import division
from numpy.testing import TestCase, assert_array_equal, assert_almost_equal, assert_equal
import numpy as N

from NewCode.tests.TestMeshes import FlatTet, TwoTets
from NewCode.tests import xfail
from NewCode import IOUtils, Mesh

class test_maxima_mesh(TestCase):
    def test_maxima_mesh(self):
        desired = """
mesh_verts:rationalize(matrix(
[-0.5, 0.5, -0.5]
,[0.5, 0.5, 0.5]
,[0.5, -0.5, -0.5]
,[-0.5, -0.5, 0.5]
,[-0.5, -0.5, -0.5]
))$

tets:[
[1, 2, 3, 4]
,[1, 3, 4, 5]
]$"""
        
        mesh = Mesh.Mesh(TwoTets.listmesh)
        assert_equal(list(IOUtils.maxima_mesh(mesh)), desired.splitlines())
