# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as N
import re, os

from FenicsCode.Utilities import MeshConverters

class UninitialisedException(Exception):
    pass

def uninit_error(fn):
    def check_fn(self):
        try: return fn(self)
        except AttributeError: raise Exception(
            "Call read_meshfile() method of %s first" % self.__class__)
    return check_fn


class FemmeshReader(object):
    block_start_re = re.compile(r'^BLOCK (.*)')
    read = False
    
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
        self.read = True

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

    
def femmesh_2_dolfin_mesh(femmesh_file):
    """Convert a femmesh file to a dolfin mesh object"""
    femmesh_reader = FemmeshReader(femmesh_file)
    return MeshConverters.femmesh_reader_2_dolfin_mesh(femmesh_reader)
