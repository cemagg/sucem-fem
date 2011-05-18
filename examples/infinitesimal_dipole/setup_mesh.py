from __future__ import division

import numpy as N
import dolfin

def get_centred_cube(domain_size, max_edge_len, centred_element_coordinate=None):
    domain_subdivisions = N.array(N.ceil(N.sqrt(2)*domain_size/max_edge_len), N.uint)
    mesh = dolfin.UnitCube(*domain_subdivisions)
    # Transform mesh to correct dimensions
    mesh.coordinates()[:] *= domain_size
    mesh.coordinates()[:] -= domain_size
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
        # There seems to be an issue with the intersect operator if the
        # mesh coordinates are changed after calling it for the first
        # time. Since we called it to find the centroid, we should init a
        # new mesh
        mesh_coords = mesh.coordinates().copy()
        mesh = dolfin.UnitCube(*domain_subdivisions)
        mesh.coordinates()[:] = mesh_coords
        mesh.coordinates()[:] -= centre_coord
        ##
    #import pdb ; pdb.set_trace()
    return mesh


    
    

