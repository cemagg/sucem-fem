from __future__ import division

import numpy as N
import dolfin

def get_centred_cube(domain_size, max_edge_len, centred_element_coordinate=None):
    domain_subdivisions = N.array(N.ceil(N.sqrt(2)*domain_size/max_edge_len), N.uint)
    mesh = dolfin.UnitCube(*domain_subdivisions)
    # Transform mesh to correct dimensions
    mesh.coordinates()[:] *= domain_size
    mesh.coordinates()[:] -= domain_size/2
    if centred_element_coordinate is not None:
        ## Translate mesh slightly so that source coordinate lies at
        ## centroid of an element
        centred_element_point = dolfin.Point(*centred_element_coordinate)
        source_elnos = mesh.all_intersected_entities(centred_element_point)
        closest_elno = source_elnos[(N.argmin(
            [centred_element_point.distance(dolfin.Cell(mesh, i).midpoint())
             for i in source_elnos]))]
        centre_pt = dolfin.Cell(mesh, closest_elno).midpoint()
        centre_coord = N.array([centre_pt.x(), centre_pt.y(), centre_pt.z()])
        # The mesh intersect operator caches certain data
        # structures. We need to clear them if the mesh coordinates
        # are changed after calling any mesh intersection methods.
        mesh.intersection_operator().clear()
        mesh.coordinates()[:] -= centre_coord
    return mesh


    
    

