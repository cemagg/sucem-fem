from __future__ import division
import dolfin

class CalculateConnectivity(object):
    def set_listmesh(self,input_listmesh):
        self.input_listmesh = input_listmesh
        self.dolfin_mesh = None

    def setup_mesh(self):
        """Setup dolfin mesh using node and tet data from self.input_listmesh"""
        dm = self.dolfin_mesh = dolfin.Mesh()
        me = dolfin.MeshEditor()
        me.open(dm, 'tetrahedron', 3, 3)
        me.init_vertices(len(self.input_listmesh['Nodes']))
        me.init_cells(len(self.input_listmesh['ElementNodes']))
        dm.coordinates()[:,:] = self.input_listmesh['Nodes']
        dm.cells()[:,:] = self.input_listmesh['ElementNodes']
        me.close()

    ### FIXME use @uninit_error -like decorator to check initialisation
    def get_nodes(self):
        return self.dolfin_mesh.coordinates()

    def get_tet_nodes(self):
        return self.dolfin_mesh.cells()
