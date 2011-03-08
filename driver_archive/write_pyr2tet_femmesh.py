from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip
import os, sys
import pickle
import random
import numpy as N
import scipy
from scipy import sparse, linalg
from numpy.testing import assert_almost_equal
#
# Local Imports
#
import NewCode
from NewCode.Utilities import Struct
from NewCode.Meshes import PyramMesh, Conversions, BrickMeshGen

h0 = 1.0001/20.
a,b,c = 29,23,19
mesh = PyramMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [a*h0, b*h0, c*h0]))

outfile = file('/tmp/tst.femmesh', 'w')
Conversions.pyramid2tet_femmesh(mesh, outfile)
outfile.close()
