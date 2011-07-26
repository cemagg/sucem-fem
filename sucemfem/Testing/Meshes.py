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
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
from __future__ import division
import numpy as N
import dolfin

class InscribedTet(object):
    def __init__(self):
        self.coordinates = N.array(
            [[-1./2 ,     1./2 ,    -1./2],
             [ 1./2 ,     1./2 ,     1./2],
             [ 1./2 ,    -1./2 ,    -1./2],
             [-1./2 ,    -1./2 ,     1./2],  
             [ 1./6 ,     1./6 ,    -1./6],
             [-1./6 ,     1./6 ,     1./6],
             [-1./6 ,    -1./6 ,    -1./6],
             [ 1./6 ,    -1./6 ,     1./6]],
            N.float64)
        self.element_nodes = N.array(
            [[1, 4, 6, 7],
             [1, 3, 5, 7],          
             [3, 4, 7, 8],          
             [1, 2, 5, 6],          
             [2, 3, 5, 8],          
             [2, 4, 6, 8],          
             [4, 6, 7, 8],          
             [1, 5, 6, 7],          
             [3, 5, 7, 8],          
             [2, 5, 6, 8],          
             [5, 6, 7, 8]], N.uint32) - 1

    def get_dolfin_mesh(self):
        mesh = dolfin.Mesh()
        me = dolfin.MeshEditor()
        me.open(mesh,'tetrahedron', 3, 3)
        me.init_vertices(len(self.coordinates))
        me.init_cells(len(self.element_nodes))
        for i, coord in enumerate(self.coordinates):
            me.add_vertex(N.uint(i), *coord)
        for i, el_nodes in enumerate(self.element_nodes):
            me.add_cell(i, *el_nodes)
        me.close()
        return mesh
