from __future__ import division

import numpy as N
import re, os

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

def uninit_error(fn):
    def check_fn(self):
        try: return fn(self)
        except AttributeError: raise Exception(
            "Call read_meshfile() method of %s first" % self.__class__)
    return check_fn

class FemmeshReader(object):
    block_start_re = re.compile(r'^BLOCK (.*)')

    def __init__(self, mesh_filename):
        self.mesh_filename = mesh_filename
        self.block_parse_funs = dict(
            nodes=self.parse_nodes,
            tets=self.parse_tets,
            # We ignore trianlgles
            tris=lambda *x: None)
        

    def parse_nodes(self, iter):
        self.no_nodes = int(iter.next().split()[0])
        node_check = set(range(self.no_nodes))
        self.nodes = N.zeros((self.no_nodes, 3))
        for i in range(self.no_nodes):
            sp = iter.next().split()
            node_no = int(sp[0]) -1
            assert node_no >= 0 and node_no < self.no_nodes
            self.nodes[node_no] = [float(x) for x in sp[1:]]
            node_check.discard(node_no)
        assert iter.next().split()[0] == 'ENDBLOCK'
        assert len(node_check) == 0

    def parse_tets(self, iter):
        self.no_tets = int(iter.next().split()[0])
        tet_check = set(range(self.no_tets))
        self.tet_nodes = N.zeros((self.no_tets, 4), dtype=N.int32)
        self.tet_property_nos = N.zeros(self.no_tets, dtype=N.int32)
        for i in range(self.no_tets):
            sp = iter.next().split()
            tet_no = int(sp[0]) -1
            assert tet_no >= 0 and tet_no < self.no_tets
            self.tet_property_nos = int(sp[1])
            self.tet_nodes[tet_no] = sorted([int(x)-1 for x in sp[2:]])
            tet_check.discard(tet_no)
        assert iter.next().split()[0] == 'ENDBLOCK'
        assert len(tet_check) == 0

    def find_block_fun(self, l):
        block_match = self.block_start_re.match(l)
        if not block_match: return None
        block_name = block_match.groups()[0]
        return self.block_parse_funs[block_name]
        
    def read_meshfile(self):
        fo = open(self.mesh_filename)
        for l in fo:
            bfun = self.find_block_fun(l)
            if bfun: bfun(fo)

    @uninit_error
    def get_mesh_filename(self):
        return os.path.basename(self.mesh_filename)

    @uninit_error
    def get_mesh_dirname(self):
        return os.path.dirname(self.mesh_filename)

    @uninit_error
    def get_tet_nodes(self):
        return self.tet_nodes

    @uninit_error
    def get_tet_property_nos(self):
        return self.tet_property_nos

    @uninit_error
    def get_nodes(self):
        return self.nodes


class Femmesh2ListMesh(object):
    def __init__(self, femmesh):
        self.listmesh = dict(
            FemmeshFilename=femmesh.get_mesh_filename(),
            FemmeshDir = femmesh.get_mesh_dirname(),
            Nodes=femmesh.get_nodes(),
            ElementNodes=femmesh.get_tet_nodes())            
    
    def get_listmesh(self):
        return self.listmesh
        
def get_femmesh_as_listmesh(meshfile):
    """Helper function that reads a femmesh file and returns a listmesh"""
    femmesh_reader = FemmeshReader(meshfile)
    femmesh_reader.read_meshfile()
    femmesh2listmesh = Femmesh2ListMesh(femmesh_reader)
    return femmesh2listmesh.get_listmesh()
