from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
import os
import pickle
from numpy import *
import numpy as N
import scipy
#
# Local Imports
#
import NewCode
from NewCode.Meshes import PyramMesh, BrickMesh, BrickMeshGen
from NewCode.Meshes import Conversions 
from NewCode.DifferentialForm.HybridMesh import get_tface_nodes
from NewCode.Meshes.MeshIO import FemmeshWriter

#tet_geom_minsize = N.array([0.25,0.25,0.25])
tet_geom_minsize = N.array([0.5,0.5,0.5])
a,b,c, = tet_geom_minsize 
h = 1/16.


listmesh = BrickMeshGen.make_rect_cavity_brick_listmesh(
    a,b,c, [h, h, h], grid_offset=[-a/2,-b/2,-c/2])

pyr_mesh = PyramMesh.Mesh(listmesh)
brick_mesh = BrickMesh.Mesh(listmesh)

tet_bdry_facenodes = N.vstack([N.vstack(get_tface_nodes(hface.nodes))
                           for hface in brick_mesh.faces.onBoundary])

femmesh_file = file('+tst_dipole_pyr.femmesh', 'w') 

tet_mesh = Conversions.pyramid2tetmesh(pyr_mesh)

femmesh_writer = FemmeshWriter(femmesh_file, tet_mesh)
femmesh_writer.addExtraTris(tet_bdry_facenodes,
                            N.ones(len(tet_bdry_facenodes), dtype=N.int32))
femmesh_writer.writeMesh()
femmesh_file.close()


                                           
