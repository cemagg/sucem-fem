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

import unittest
import dolfin
import numpy as np

from FenicsCode.Testing.Meshes import InscribedTet
# Module under test
from FenicsCode import Geometry 

class test_BoundaryEdges(unittest.TestCase):
    def setUp(self):
        self.mesh = dolfin.UnitCube(1,1,1)
        self.no_boundary_edges = 18 # for 1x1x1 UnitCube meshed with 2
                                    # tris on each face
        self.DUT = Geometry.BoundaryEdges(self.mesh)
        
    def test_mark(self):
        edge_meshfunc = dolfin.EdgeFunction('uint', self.mesh)
        edge_meshfunc.set_all(0)
        test_val = 11
        self.DUT.mark(edge_meshfunc, test_val)
        no_edges = len(edge_meshfunc.array())
        self.assertTrue(
            sum(edge_meshfunc.array() == test_val) == self.no_boundary_edges)
        for edge in dolfin.SubsetIterator(edge_meshfunc, test_val):
            self.assertTrue(self._check_on_boundary(edge))

    def _check_on_boundary(self, entity):
        """Check if an entity is on boundary using numerical node positions

        For the unit cube geometry an entity is on a boundary face if
        all the nodes have one coordinate value that is equal to 0 or
        one. Eg. if all the x values are equal to 1 that means the
        entity lies on the x=0 boundary face.
        """
        node_coords = entity.mesh().coordinates()[entity.entities(0)]
        for coordvals in node_coords.T:
            if np.allclose(coordvals,0) or np.allclose(coordvals, 1):
                return True
        return False

class test_BoundaryEdgeCells(unittest.TestCase):
    def setUp(self):
        self.mesh = InscribedTet().get_dolfin_mesh()
        self.DUT = Geometry.BoundaryEdgeCells(self.mesh)

    def test_mark(self):
        mark_value = 3
        desired_edgecells = np.ones(self.mesh.num_cells())*mark_value
        # All except last tet are edge-connected to the boundary for
        # the inscribed tet mesh
        desired_edgecells[-1] = 0
        cell_fn = dolfin.CellFunction('uint', self.mesh)
        self.DUT.mark(cell_fn, mark_value)
        self.assertTrue(np.all(cell_fn.array() == desired_edgecells))
        

    
