from __future__ import division

import numpy as N

from NewCode.Utilities import Struct

def get_femmesh_tris(femmesh_file):
    tris_found = False
    for l in femmesh_file:
        if l == 'BLOCK tris\n':
            tris_found = True
            break
    if not tris_found: raise ValueError('No triangle block found')
    no_tris = int(femmesh_file.next())
    tri_nodes = N.zeros((no_tris, 3), N.int32)
    tri_material = N.zeros(no_tris, N.int32)
    for i,l in enumerate(femmesh_file):
        if l == 'ENDBLOCK\n': break
        lsplit = l.split()
        tri_material[i] = int(lsplit[1])
        tri_nodes[i] = sorted([int(nn) for nn in lsplit[2:]])
    tri_nodes -= 1
    return Struct(nodes=tri_nodes, material=tri_material)
