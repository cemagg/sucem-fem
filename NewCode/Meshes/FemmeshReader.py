from __future__ import division

import numpy as N
import re

#def nodes

class FemmeshReader(object):
    block_start_re = re.compile(r'^BLOCK (.*)')

    def __init__(self, mesh_filename):
        self.mesh_filename = mesh_filename
        self.block_parse_funs = dict(
            nodes=self.parse_nodes,
            tets=self.parse_tets)
        

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
        pass

    def find_block_fun(self, l):
        block_match = self.block_start_re.match(l)
        if not block_match: return None
        block_name = block_match.groups()[0]
        
        
    def read_meshfile(self):
        fo = open(self.mesh_filename)
        for l in fo:
            bfun = self.find_block(l)
            if bfun: bfun(fo)
            

    def get_listmesh(self):
        try: return self.listmesh
        except AttributeError: raise Exception("Call read_meshfile() method first")
