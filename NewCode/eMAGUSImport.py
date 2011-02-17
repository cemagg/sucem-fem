"""
Import geometry structures from eMAGUS data structures.

Module defs
===========

init
  Calls eMAGUS to load the mesh and initialises the data structures

get_listmesh
  Returns a dict with the meshlists in the format acceptable to Mesh.Mesh

get_mesh
  Copies mesh structures from FEMFEKO, and returns a Mesh.Mesh class instance

"""

from numpy import transpose, int32, array
#
# Local Imports
#
from femfeko import femfeko, geomwrap, geometry, parseinput
import os
import Mesh

MAIN_EXECUTED = False
MESH_CWD = ''

def init(cwd=None):
    global MAIN_EXECUTED
    global MESH_CWD
    if MAIN_EXECUTED:
        print "Fortran main routine already executed"
    else:
        old_cwd = os.getcwd()
        if cwd:
            os.chdir(cwd)
            MESH_CWD = os.getcwd()
        else: MESH_CWD = os.getcwd()
        try:
            femfeko.run_femfeko()
        finally:
            os.chdir(old_cwd)
        MAIN_EXECUTED = True

def get_mesh():
    return Mesh.Mesh(get_listmesh())

def get_listmesh():
    geomwrap.init_geom()
    # Make copies so we are not dependent on the original Fortran storage,
    # transpose to get row-major storage and zero-base the indices.  Note it is
    # not neccesary to copy() when subtracting 1, since a copy is created in any
    # case by the subtraction.
    mesh = {}
    mesh['NodeConnect2Element'] = array(geomwrap.node_elements).astype(int32) - 1
    mesh['NodeConnect2ElementPtr'] = array(
        geomwrap.node_element_ptr[1:]).astype(int32).copy()
    mesh['ElementNodes'] = transpose(geomwrap.element_nodes).astype(int32) - 1
    mesh['ElementEdges'] = transpose(geomwrap.element_edges).astype(int32) - 1
    mesh['ElementFaces'] = transpose(geomwrap.element_faces).astype(int32) - 1
    mesh['FaceNodes'] = transpose(geomwrap.face_nodes).astype(int32) - 1
    mesh['EdgeNodes'] = transpose(geomwrap.edge_nodes).astype(int32) - 1
    mesh['Nodes'] = transpose(geomwrap.vertex_coords).copy()
    mesh['ElementConnect2Elem'] = transpose(geomwrap.element_connect2elem).astype(int32) - 1
    mesh['ElementConnect2Face'] = transpose(geomwrap.element_connect2face).astype(int32) - 1
    mesh['FaceConnect2Elem'] = array(geometry.faceconnectelem, int32) - 1
    mesh['EdgeConnect2Elem'] = array(geometry.edgeconnectelem, int32) - 1
    mesh['FemmeshFilename'] = parseinput.extra_meshfilename.tostring().strip()
    mesh['FemmeshDir'] = MESH_CWD
    return mesh
                                   
