from __future__ import division
"""
Test driver code. Anything of lasting importance should be factored out pronto
"""
from itertools import izip, chain
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
from NewCode import Utilities, ProxyList
from NewCode.Utilities import Struct,  close_to_point
import NewCode.Mesh as Mesh
from NewCode.Meshes import BrickMesh, Conversions, BrickMeshGen
from NewCode.DifferentialForm import Discretiser, allfree

order = 1

# h0 = 1.0001/1.
# ha = 1.0001/2.
a,b,c = 29,23,19
h = a/1.
bfrac = 0.5
# hex_mesh = BrickMesh.Mesh(
#     BrickMeshGen.make_rect_cavity_brick_listmesh(a*bfrac,b,c, [a*ha, b*h0, c*h0]))
hex_mesh = BrickMesh.Mesh(
    BrickMeshGen.make_rect_cavity_brick_listmesh(a*bfrac,b,c,h))

print 'Hex-Mesh elements: ', len(hex_mesh.elements)

g_eps = 1e-10                           # Geometrical tollerance
hybrid_boundary_p = close_to_point(a*bfrac, g_eps)
b_p = close_to_point(b, g_eps)
c_p = close_to_point(c, g_eps)
zero_p = close_to_point(0, g_eps)
on_hbdry = lambda ent: N.all([hybrid_boundary_p(x) for (x,y,z) in ent.nodeCoords])
on_hb_geom_edge_1 = lambda ent: on_hbdry(ent) and \
                    N.all([zero_p(z) for (x,y,z) in ent.nodeCoords])
on_hb_geom_edge_2 = lambda ent: on_hbdry(ent) and \
                    N.all([b_p(y) for (x,y,z) in ent.nodeCoords])
on_hb_geom_edge_3 = lambda ent: on_hbdry(ent) and \
                    N.all([c_p(z) for (x,y,z) in ent.nodeCoords])
on_hb_geom_edge_4 = lambda ent: on_hbdry(ent) and \
                    N.all([zero_p(y) for (x,y,z) in ent.nodeCoords])



tst_geof = file('+tst2.geo', 'w')

def nodes_to_geo(nodes, offset=0,h=1000000.):
    i = offset+1
    for node in nodes:
        yield "Point(%d) = {%.16g,%.16g,%.16g,%.16g};\n" % (
            i, node[0], node[1], node[2], h)
        i += 1

bdry_faces = N.array([i for i,f in enumerate(hex_mesh.faces) if on_hbdry(f)],
                     N.int32)

bdry_edges = N.array([i for i,e in enumerate(hex_mesh.edges) if on_hbdry(e)],
                     N.int32)

# joining_nodes = N.array([N.average(f.nodeCoords, axis=0) + [a*h0/2,0,0]
#                  for f in hex_mesh.faces[bdry_faces]])

joining_nodes = N.array([N.average(f.nodeCoords, axis=0) + [h/2,0,0]
                         for f in hex_mesh.faces[bdry_faces]])
def edges_to_geo(edges, offset=0):
    i = offset+1
    for edge in edges:
        yield "Line(%d) = {%.16g,%.16g};\n" % (i, edge.nodes[0]+1, edge.nodes[1]+1)
        i += 1

def edgenodes_to_geo(edgenodes, offset=0):
    i = offset+1
    for edgends in edgenodes:
        yield "Line(%d) = {%.16g,%.16g};\n" % (i, edgends[0]+1, edgends[1]+1)
        i += 1

def faceedges_to_geo(faceedges, offset=0):
    i = offset+1
    for f in faceedges:
        yield "Line Loop(%d) = {%d,%d,-%d};\n" % (i, f[0]+1,f[1]+1,f[2]+1)
        yield "Plane Surface(%d) = {%d};\n" % (i, i)
        i += 1
                
def hexfaceedges_to_geo(faceeedges, offset=0):
    i = offset+1
    ap1 = lambda x: (N.abs(x)+1)*N.sign(x)
    for f in faceeedges:
        f = ap1(f)
        f[f==0]=1
        substr = "Line Loop(%d) = {" + ','.join('%d' for i in range(len(f))) \
                 + "};\n"
        yield substr % tuple(nn for nn in chain((i,),f))
        yield "Plane Surface(%d) = {%d};\n" % (i, i)
        i += 1

def faces_to_geo(faces, offset=0):
    i = offset+1
    ap1 = lambda x: (N.abs(x)+1)*N.sign(x)
    for f in faces:
        f = ap1(f)
        f[f==0]=1
        substr = "Surface Loop(%d) = {" + ','.join('%d' for i in range(len(f))) \
                 + "};\n"
        yield substr % tuple(nn for nn in chain((i,),f))
        yield "Volume(%d) = {%d};\n" % (i, i)
        
def tetfaces_to_geo(tetfaces, offset=0):
    i = offset+1
    for t in tetfaces:
        yield "Surface Loop(%d) = {%d,%d,%d,%d};\n" % (i, t[0]+1,t[1]+1,t[2]+1,t[3]+1)
        yield "Volume(%d) = {%d};\n" % (i, i)
        i += 1

FakeEdge = Mesh.ProxyEdge

def make_diag_edgenodes(bdry_faces):
    nodes = N.array([f.nodes[[0,2]] for f in bdry_faces], N.int32)
    return nodes

def make_joining_edgenodes(bdry_faces, joiningnode_offset=0):
    nodes_per_face = 4
    nodes = N.zeros((len(bdry_faces)*nodes_per_face, 2), dtype=N.int32)
    ei = 0
    for fi, f in enumerate(bdry_faces):
        for node in f.nodes:
            nodes[ei] = node, fi+joiningnode_offset
            ei+=1
    return nodes

def make_joining_elementnodes(bdry_faces, joiningnode_offset=0):
    nodes_per_tet = 4 ; tets_per_bdryface = 2
    tet_elnodes = N.zeros((len(bdry_faces)*tets_per_bdryface, nodes_per_tet))
    for fi, f in enumerate(bdry_faces):
        fnodes= f.nodes
        j_node = fi+joiningnode_offset
        tet_elnodes[fi*tets_per_bdryface] = N.hstack((fnodes[[0,3,2]], [j_node]))
        tet_elnodes[fi*tets_per_bdryface+1] = N.hstack((fnodes[[0,1,2]], [j_node]))
    return tet_elnodes

def make_joining_facenodes(hex_facenodes):
    tris_per_hex=2 ; nodes_per_face = 3
    facenodes = N.zeros((len(hex_facenodes)*tris_per_hex, nodes_per_face), N.int32)
    fi = 0
    hfaceiter = hex_facenodes.__iter__()
    for nns in hfaceiter:
        facenodes[fi] = nns[[0,3,2]]
        fi += 1
        facenodes[fi] = nns[[0,1,2]]
        fi += 1
    return facenodes

def make_faceedges(facenodes, edgemap):
    edges_per_face = 3
    faceedges = N.zeros((len(facenodes),edges_per_face), N.int32)
    fi = 0
    for fnds in facenodes:
        faceedges[fi] = [edgemap[tuple(fnds[[0,1]])],
                         edgemap[tuple(fnds[[1,2]])],
                         edgemap[tuple(fnds[[0,2]])]]
        fi += 1
    return faceedges

def make_tetfaces(tetnodes, facemap):
    from NewCode.Mesh import Element
    tet_faces = N.zeros((len(tetnodes), Element.COUNT.face), N.int32)
    for ti, tetnd in enumerate(tetnodes):
        for lfi, lfn in enumerate(Element.LOCAL_FACENODES):
            tet_faces[ti, lfi] = facemap[tuple(tetnd[lfn])]

    return tet_faces
            
hyb_nodes = hex_mesh.nodes[:]
hexface_edgenodes = [e.nodes for e in hex_mesh.edges if on_hbdry(e)]
diag_edgenodes = make_diag_edgenodes(hex_mesh.faces[bdry_faces])

hyb_edgenodes = N.vstack((hexface_edgenodes, diag_edgenodes))
hyb_edges = ProxyList.ProxyList(FakeEdge(dict(nodes=hyb_edgenodes,
                                              nodeCoords=hyb_nodes)))

hyb_facenodes = make_joining_facenodes(hex_mesh.faces[bdry_faces].nodes)

hyb_edgemap = dict((tuple(n), i) for i,n in enumerate(hyb_edgenodes))
hyb_facemap = dict((tuple(n), i) for i,n in enumerate(hyb_facenodes))

hyb_faceedges = make_faceedges(hyb_facenodes, hyb_edgemap)

no_hexnodes = len(hex_mesh.nodes[:])
geom_nodes = N.array([[a,0,0], [a,b,0], [a,b,c], [a,0,c]], N.float64)
geom_nodes_nos = gnns = N.arange(len(hyb_nodes), len(hyb_nodes)+len(geom_nodes),
                          dtype=N.int32)
hyb_corner_nodes = N.array([hex_mesh.findClosestNode(r)
                    for r in N.array([[a*bfrac,0,0], [a*bfrac,b,0],
                                      [a*bfrac,b,c], [a*bfrac,0,c]])],
                           N.int32)
geom_edgenodes = N.array([[hyb_corner_nodes[0],gnns[0]],
                          [hyb_corner_nodes[1],gnns[1]],
                          [hyb_corner_nodes[2],gnns[2]],
                          [hyb_corner_nodes[3],gnns[3]],
                          gnns[[0,1]],
                          gnns[[1,2]],
                          gnns[[2,3]],
                          gnns[[0,3]]], N.int32)
geom_edge_offs = len(hyb_edges)
geom_edge_nos = gens = N.arange(len(hyb_edges), len(hyb_edges)+len(geom_edgenodes),
                          dtype=N.int32)

hgeom_lineloop_1 = N.array([i for i, edge in enumerate(hyb_edges)
                            if on_hb_geom_edge_1(edge)], N.int32)
hgeom_lineloop_2 = N.array([i for i, edge in enumerate(hyb_edges)
                            if on_hb_geom_edge_2(edge)], N.int32)
hgeom_lineloop_3 = N.array([i for i, edge in enumerate(hyb_edges)
                            if on_hb_geom_edge_3(edge)], N.int32)
hgeom_lineloop_4 = N.array([i for i, edge in enumerate(hyb_edges)
                            if on_hb_geom_edge_4(edge)], N.int32)
geom_face_offs = len(hyb_faceedges)
hgeom_endface_edges = N.hstack((gens[[4,5,6]], -gens[[7]]))
hgeom_sideface_1 = N.hstack((hgeom_lineloop_1, gens[[1]], -gens[[4,0]]))
hgeom_sideface_2 = N.hstack((hgeom_lineloop_2, gens[[2]], -gens[[1,5]]))
hgeom_sideface_3 = N.hstack((-gens[[3]], gens[[2,6]], hgeom_lineloop_3))
hgeom_sideface_4 = N.hstack((-gens[[0,7]], gens[[3]], hgeom_lineloop_4))
geom_faceedges = [hgeom_endface_edges,
                  hgeom_sideface_1,
                  hgeom_sideface_2,
                  hgeom_sideface_3,
                  hgeom_sideface_4]

geom_face_offset = len(hyb_faceedges)
geom_vol_faces = N.hstack((N.arange(geom_face_offset),
                           N.arange(len(geom_faceedges))+geom_face_offs))


#tst_geof.writelines(nodes_to_geo(hyb_nodes,h=N.sqrt(b**2+c**2)*h0*.75))
tst_geof.writelines(nodes_to_geo(hyb_nodes,h=h))
tst_geof.writelines(edges_to_geo(hyb_edges))
tst_geof.writelines(faceedges_to_geo(hyb_faceedges))
tst_geof.writelines(nodes_to_geo(geom_nodes, offset=len(hyb_nodes), h=h*0.2))
tst_geof.writelines(edgenodes_to_geo(geom_edgenodes, offset=geom_edge_offs))
tst_geof.writelines(hexfaceedges_to_geo(geom_faceedges, offset=geom_face_offs))
tst_geof.writelines(faces_to_geo([geom_vol_faces]))
tst_geof.close()
