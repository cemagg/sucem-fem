## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import numpy as N
import dolfin

def get_centred_cube(domain_size, max_edge_len, centred_element_coordinate=None):
    """ Generate a cube mesh centred arount [0,0,0]

    optionally translates mesh slightly such that
    centred_element_coordinate is at the centre of an element.
    """
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


    
    

