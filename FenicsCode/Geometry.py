# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import dolfin

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
        been initialised, and initalises the connectivity if needed
        """
        ent0 = self.dim_entity_map[dim0](self.dolfin_mesh, 0)
        if len(ent0.entities(dim1)) == 0: self.dolfin_mesh.init(dim0, dim1)


class BoundaryEdges(object):
    def __init__(self, mesh, boundary_facefun=None, boundary_value=1):
        """Find and mark the edges that occur on a set of boundary faces
    
        @param mesh: A dolfin mesh object
        @keyword boundary_facefun: a face mesh function that describes the
            faces to be considered 
            (default: None).
        @keyword boundary_value: the mesh function value marking the boundary faces
            (default: 1)
        """
        self.mesh = mesh
        if boundary_facefun is None:
            boundary_facefun = dolfin.FaceFunction('uint', mesh)
            boundary_facefun.set_all(0)
            domain_boundary = dolfin.DomainBoundary()
            domain_boundary.mark(boundary_facefun, boundary_value)
        self.boundary_facefun = boundary_facefun
        self.boundary_value = boundary_value
        self.ensure_initialised = EnsureInitialised(self.mesh)
        
    def mark(self, edge_meshfun, value):
        edge_array = edge_meshfun.array()
        self.ensure_initialised(2,1)
        self.ensure_initialised(1,0)
        for face in dolfin.SubsetIterator(self.boundary_facefun,
                                          self.boundary_value):
            edge_array[face.entities(1)] = value

class BoundaryEdgeCells(object):
    """Find and mark cells connected the domain boundary through edges"""
    def __init__(self, mesh):
        self.edge_boundary_value = 1
        boundary_edgefun = dolfin.EdgeFunction('uint', mesh)
        boundary_edgefun.set_all(0)
        boundary_edges = BoundaryEdges(mesh)
        boundary_edges.mark(boundary_edgefun, self.edge_boundary_value)
        self.boundary_edgefun = boundary_edgefun

    def mark(self, cell_meshfun, value):
        cells_connected2edges = CellsConnected2Edges(
            self.boundary_edgefun, self.edge_boundary_value)
        cells_connected2edges.mark(cell_meshfun, value)
        

class CellsConnected2Edges(object):
    def __init__(self, boundary_edgefun, boundary_value):
        self.boundary_edgefun = boundary_edgefun
        self.boundary_value = boundary_value
        self.mesh = boundary_edgefun.mesh()
        self.ensure_initialised = EnsureInitialised(self.mesh)

    def mark(self, cell_meshfun, value):
        cell_array = cell_meshfun.array()
        self.ensure_initialised(1,3)    # edge -> cell connectivity
        for edge in dolfin.SubsetIterator(self.boundary_edgefun,
                                          self.boundary_value):
            cell_array[edge.entities(3)] = value
         

