from __future__ import division
import dolfin
import numpy

class EnsureInitialised(object):
    """Ensures that a dolfin.mesh object connectivity info is initialised"""
    dim_entity_map = {
        0:dolfin.Vertex,
        1:dolfin.Edge,
        2:dolfin.Face,
        3:dolfin.Cell
        }
    
    def __init__(self, dolfin_mesh):
        self.dolfin_mesh = dolfin_mesh

    def __call__(self, dim0, dim1):
        """Ensure that connectivity from dim0 to dim1 is initialised

        A test is done to determine wether the connectivity info has
        been initialised, and initalised the connectivity if needed
        """
        ent0 = self.dim_entity_map[dim0](self.dolfin_mesh, 0)
        if len(ent0.entities(dim1)) == 0: self.dolfin_mesh.init(dim0, dim1)
    
class CalculateConnectivity(object):
    def set_input_listmesh(self,input_listmesh):
        self.input_listmesh = input_listmesh
        self.output_listmesh = dict(input_listmesh)
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
        self.ensure_initialised = EnsureInitialised(dm)

    ### FIXME use @uninit_error -like decorator to check initialisation
    def get_nodes(self):
        return self.dolfin_mesh.coordinates()

    def get_tet_nodes(self):
        return self.dolfin_mesh.cells()

    def calc_node_element_connectivity(self):
        """Caluculate NodeConnect2Element and NodeConnect2ElementPtr"""
        self.ensure_initialised(1,3)
        vert_cells = list(vert.entities(3) for vert in dolfin.vertices(self.dolfin_mesh))
        self.node_connect_2_element = numpy.hstack((vert_cells))
        self.node_connect_2_element_ptr = numpy.cumsum(
            [0] + [len(x) for x in vert_cells])

        self.output_listmesh['NodeConnect2Element'] = self.node_connect_2_element
        self.output_listmesh['NodeConnect2ElementPtr'] = self.node_connect_2_element_ptr

    def get_node_connect_2_element(self):
        return self.node_connect_2_element

    def get_node_connect_2_element_ptr(self):
        return self.node_connect_2_element_ptr

    def calc_element_edge_connectivity(self):
        """Calculate ElementEdges"""
        self.ensure_initialised(3,1)
        # As defined in eMAGUS
        target_local_edgenodes = [tuple(nodepair) for nodepair in numpy.array(
            [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], numpy.int32) - 1
                                  ]
        # As defined in UFC user manual
        UFC_local_edgenodes_map = dict(
            (nodepair,ind) for ind,nodepair in enumerate(
                ((2,3), (1,3), (1,2), (0,3), (0,2), (0,1))) )
        # Permutation of UFC -> eMAGUS local edge ordering
        perm = [UFC_local_edgenodes_map[nodepair]
                for nodepair in target_local_edgenodes]
        self.element_connect_2_edge = numpy.vstack(
            [cell.entities(1)[perm] for cell in dolfin.cells(self.dolfin_mesh)])
                                                    
    def get_element_connect_2_edge(self):
        return self.element_connect_2_edge

    def calc_edge_node_connectivity(self):
        """Calculate EdgeNodes"""
        self.ensure_initialised(1,0)
        self.edge_connect_2_node = numpy.vstack(
            [edge.entities(0) for edge in dolfin.edges(self.dolfin_mesh)])

    def get_edge_connect_2_node(self):
        return self.edge_connect_2_node

    def calc_face_node_connectivity(self):
        """Calculate FaceNodes"""
        self.ensure_initialised(2,0)
        self.face_connect_2_node = numpy.vstack(
            [face.entities(0) for face in dolfin.faces(self.dolfin_mesh)])

    def get_face_connect_2_node(self):
        return self.face_connect_2_node

    def calc_element_face_connectivity(self):
        """Calculate ElementFaces"""
        self.ensure_initialised(3,2)
        # As defined in eMAGUS
        target_local_facenodes = [tuple(nodepair) for nodepair in numpy.array(
            [[1,2,3], [1,2,4], [1,3,4], [2,3,4]], numpy.int32) - 1
                                  ]
        # As defined in UFC user manual
        UFC_local_facenodes_map = dict(
            (nodepair,ind) for ind,nodepair in enumerate(
                ((1,2,3), (0,2,3), (0,1,3), (0,1,2) )) )
        # Permutation of UFC -> eMAGUS local face ordering
        perm = [UFC_local_facenodes_map[nodepair]
                for nodepair in target_local_facenodes]
        self.element_connect_2_face = numpy.vstack(
            [cell.entities(2)[perm] for cell in dolfin.cells(self.dolfin_mesh)])
                                                    
    def get_element_connect_2_face(self):
        return self.element_connect_2_face

    def get_output_listmesh(self):
        return self.output_listmesh
