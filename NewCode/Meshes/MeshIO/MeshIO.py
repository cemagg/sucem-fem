from __future__ import division

from itertools import izip

class FemmeshWriter(object):
    def __init__(self, outmeshfile, inmesh):
        self.outmeshfile = outmeshfile
        self.inmesh = inmesh

    def writeMesh(self):
        self.writeNodes()
        if hasattr(self, 'extra_trinodes'):
            self.writeExtraTris()
        self.writeTets()

    def addExtraTris(self, trinodes, mat_property):
        self.extra_trinodes = trinodes
        self.extra_tri_mat_property = mat_property

    def writeExtraTris(self):
        outmeshfile = self.outmeshfile
        trinodes = self.extra_trinodes
        tri_mats = self.extra_tri_mat_property
        print >> outmeshfile, "BLOCK tris"
        print >> outmeshfile, len(trinodes)

        tri_index = 1
        for tri, tri_mat in izip(trinodes, tri_mats):
            nodes = [i + 1 for i in tri] # femmesh format is 1-based, not 0
            print >> outmeshfile, tri_index, tri_mat,\
                  nodes[0], nodes[1], nodes[2]
            tri_index += 1

        print >> outmeshfile, "ENDBLOCK"
    

    def writeTets(self):
        """
        Writes tetraheders to the femmesh output file
        """

        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK tets"
        print >> outmeshfile, len(inmesh.elements)

        tet_index = 1
        material = 0
        for tet in inmesh.elements:
            nodes = [i + 1 for i in tet.nodes] # femmesh format is 1-based, not 0
            print >> outmeshfile, tet_index, material,\
                  nodes[0], nodes[1], nodes[2], nodes[3]
            tet_index += 1

        print >> outmeshfile, "ENDBLOCK"
    

    def writeNodes(self):
        """
        Writes nodes to the femmesh output file
        """
        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK nodes"
        print >> outmeshfile, len(inmesh.nodes)

        node_index = 1
        for node in inmesh.nodes:
            print >> outmeshfile, node_index, node[0], node[1], node[2]
            node_index += 1

        print >> outmeshfile, "ENDBLOCK"

