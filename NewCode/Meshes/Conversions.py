from __future__ import division

import numpy as N

from NewCode import ProxyList
from NewCode.Utilities import Struct
from NewCode.Meshes import MeshIO

def pyramid_els2tet_els(pyram_element_nodes):
    tet_nodes = []
    for peln in pyram_element_nodes:
        tet_nodes.append(peln[[0,2,3,4]])
        tet_nodes.append(peln[[0,1,3,4]])
    return N.array(tet_nodes, dtype=pyram_element_nodes.dtype)

def pyramid2tet_femmesh(pyram_mesh, outfile):
    fake_Elclass = ProxyList.ItemClassFactory(['nodes'], 'fake_Elclass')
    tet_mesh = Struct(nodes=pyram_mesh.nodes,
                      elements=ProxyList.ProxyList(fake_Elclass(
        {'nodes': pyramid_els2tet_els(pyram_mesh.elements[:].nodes)})))
    
    meshwriter = MeshIO.FemmeshWriter(outfile, tet_mesh)
    meshwriter.writeMesh()
    
def pyramid2tetmesh(pyram_mesh):
    fake_Elclass = ProxyList.ItemClassFactory(['nodes'], 'fake_Elclass')
    tet_mesh = Struct(nodes=pyram_mesh.nodes,
                      elements=ProxyList.ProxyList(fake_Elclass(
        {'nodes': pyramid_els2tet_els(pyram_mesh.elements[:].nodes)})))
    return tet_mesh

def hexfaces2diagnodes(hexfacenodes):
    """Construct diagonal edges on rectangular faces and relattion to non-diagonal edges

    Arguments
    ==========

    hexfacenodes

      Sequence of rectangular faces defined on a hexahedral mesh that form the boundary
      between a hexahedral and tetrahedral mesh region.

    Return Dictionary

      keys: Diagonal edge node numbers

      values: (s_edgenodes, e_edgenodes) where p_edgenodes are the node numbers of
      hex edges that are connected to the starting node of the diagonal and
      n_edgenodes those connected to the ending node of of the diagonal
      """
    # TESTME !!!
    d_edgemap = {}
    for fnodes in hexfacenodes:
        d_edgenodes = fnodes[[0,2]]      # hex face diagonal nodes
        s_edgenodes = N.array([fnodes[[0,1]], fnodes[[0,3]]], N.int32)
        e_edgenodes = N.array([fnodes[[3,2]], fnodes[[1,2]]], N.int32)
        d_edgemap[tuple(d_edgenodes)] = s_edgenodes, e_edgenodes
    return d_edgemap

import dolfin

def listmesh_2_dolfin_mesh(listmesh):
    """Setup dolfin mesh using node and tet data from listmesh"""
    dm = dolfin.Mesh()
    me = dolfin.MeshEditor()
    me.open(dm, 'tetrahedron', 3, 3)
    me.init_vertices(len(listmesh['Nodes']))
    me.init_cells(len(listmesh['ElementNodes']))
    dm.coordinates()[:,:] = listmesh['Nodes']
    dm.cells()[:,:] = listmesh['ElementNodes']
    me.close()
    return dm

def dolfin_mesh_2_listmesh(dolfin_mesh):
    listmesh = {}
    listmesh['Nodes'] = dolfin_mesh.coordinates().copy()
    listmesh['ElementNodes'] = dolfin_mesh.cells().copy()
    return listmesh
