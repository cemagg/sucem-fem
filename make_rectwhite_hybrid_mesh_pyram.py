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
from NewCode.Meshes.MeshIO import FemmeshWriter, Gmsh

a,b,c, = geom_size = N.array([29,23,19], N.float64)
h0 = 1/1.
h = geom_size*h0
#tet_geom_size_q = N.array([2,2,2], N.int32)
#tet_geom_size_q = N.array([1/h0/2,1/h0,1/h0], N.int32)
h[0] /= 2 ; tet_geom_size_q = N.array([1,1,1], N.int32)


#tet_offset_q = [1,1,1]
tet_offset_q = [0,0,0]
a_b,b_b,c_b, = brick_size = tet_geom_size_q*h

listmesh = BrickMeshGen.make_rect_cavity_brick_listmesh(
    a_b, b_b, c_b, h, grid_offset=tet_offset_q*h)

pyr_mesh = PyramMesh.Mesh(listmesh)
brick_mesh = BrickMesh.Mesh(listmesh)

print "brick els: ", len(brick_mesh.elements)

tet_bdry_facenodes = N.vstack([N.vstack(get_tface_nodes(hface.nodes))
                           for hface in brick_mesh.faces.onBoundary])

femmesh_file = file('+tst_rectwhite_pyr.femmesh', 'w') 
tet_mesh = Conversions.pyramid2tetmesh(pyr_mesh)

femmesh_writer = FemmeshWriter(femmesh_file, tet_mesh)
femmesh_writer.addExtraTris(tet_bdry_facenodes,
                            N.ones(len(tet_bdry_facenodes), dtype=N.int32))
femmesh_writer.writeMesh()
femmesh_file.close()

gmsh_file = file('+tst_rectwhite_pyr.msh', 'w')
gmsh_file.writelines(Gmsh.meshformat())
gmsh_file.writelines(Gmsh.nodes_to_msh(tet_mesh.nodes))
gmsh_file.writelines(Gmsh.tets_faces_to_msh(tetnodes=tet_mesh.elements[:].nodes,
                                            trinodes=tet_bdry_facenodes))
gmsh_file.close()
