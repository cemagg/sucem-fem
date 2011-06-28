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
    """Find and mark the edges that occur on a set of boundary faces
    Parameters
    ==========

    mesh -- dolfin mesh object

    boundary_facefun -- optional face mesh function that describes the
        faces to be considered

    boundary_value -- mesh function value marking the boundary faces
    """
    def __init__(self, mesh, boundary_facefun=None, boundary_value=1):
        self.mesh = mesh
        if boundary_facefun is None:
            boundary_facefun = dolfin.FaceFunction('uint', mesh)
            boundary_facefun.set_all(0)
        self.boundary_facefun = boundary_facefun
        self.boundary_value = boundary_value
        domain_boundary = dolfin.DomainBoundary()
        domain_boundary.mark(boundary_facefun, boundary_value)
        self.ensure_initialised = EnsureInitialised(self.mesh)
        
    def mark(self, edge_meshfun, value):
        edge_array = edge_meshfun.array()
        self.ensure_initialised(2,1)
        self.ensure_initialised(1,0)
        for face in dolfin.SubsetIterator(self.boundary_facefun,
                                          self.boundary_value):
            edge_array[face.entities(1)] = value

