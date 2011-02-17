from __future__ import division

from itertools import chain

import numpy as N

def nodes_to_geo(nodes, offset=0,h=1000000.):
    i = offset+1
    for node in nodes:
        yield "Point(%d) = {%.16g,%.16g,%.16g,%.16g};\n" % (
            i, node[0], node[1], node[2], h)
        i += 1

def edges_to_geo(edges, offset=0):
    i = offset+1
    for edge in edges:
        yield "Line(%d) = {%.16g,%.16g};\n" % (i, edge.nodes[0]+1, edge.nodes[1]+1)
        i += 1

def circsegs_to_geo(circ_segs):
    for seg in circ_segs:
        yield "Circle(newc) = {%.16g,%.16g,%.16g};\n" % (seg[0]+1, seg[1]+1, seg[2]+1)
        
def hexfaceedges_to_geo(faceeedges, offset=0):
    i = offset+1
    for f in faceeedges:
        substr = "Line Loop(%d) = {" + ','.join('%d' for i in range(len(f))) \
                 + "};\n"
        yield substr % tuple(nn for nn in chain((i,),f))
        yield "Plane Surface(%d) = {%d};\n" % (i, i)
        i += 1

def faces_to_geo(faces, offset=0):
    i = offset+1
    for f in faces:
        substr = "Surface Loop(%d) = {" + ','.join('%d' for i in range(len(f))) \
                 + "};\n"
        yield substr % tuple(nn for nn in chain((i,),f))
        yield "Volume(%d) = {%d};\n" % (i, i)

def meshformat():
    strings = (
        "$MeshFormat\n",
        "2 0 8\n",
        "$EndMeshFormat\n")
    for s in strings:
        yield s


def node_to_strs(node):
    return [str(u) for u in node]

def nodes_to_msh(nodes, offset=0):
    yield "$Nodes\n"
    yield str(len(nodes)) + '\n'
    for i, n in enumerate(nodes):
        strs = [str(i+1+offset)]
        strs += node_to_strs(n)
        strs += '\n'
        yield " ".join(strs)
    yield "$EndNodes\n"

def tets_faces_to_msh(tetnodes=[], trinodes=[]):
    yield "$Elements\n"
    no_els = 0
    no_els += len(tetnodes)
    no_els += len(trinodes)
    yield str(no_els) + '\n'
    tri_eltype = 2 ; tet_eltype = 4 ; no_tags = 3
    elno = 1
    for tri in trinodes:
        strs = [str(elno), str(tri_eltype)]
        strs += [str(no_tags)] + ['0' for i in range(no_tags)]
        strs += node_to_strs(tri+1)
        strs += '\n'
        yield " ".join(strs)
        elno += 1
    
    for tet in tetnodes:
        strs = [str(elno), str(tet_eltype)]
        strs += [str(no_tags)] + ['0' for i in range(no_tags)]
        strs += node_to_strs(tet+1)
        strs += '\n'
        yield " ".join(strs)
        elno += 1
    
    yield("$EndElements\n")
