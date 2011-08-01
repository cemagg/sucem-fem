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
# Authors
# Neilen Marais <nmarais@gmail.com>
from __future__ import division

import itertools
import dolfin
import numpy as N
from sucemfem.Sources.current_source import CurrentSource

class PointCurrentSource(CurrentSource):
    def set_position(self, position):
        """Set point source position. Expects an array with x,y,z coordinates
        """
        self.position = N.array(position, dtype=N.float64)

    def set_value(self, value):
        """Set point value. Expects an array with x,y,z components of current
        """
        iscomplex = N.any(N.iscomplex(value))
        dtype = N.complex128 if iscomplex else N.float64
        self.value = N.array(value, dtype=dtype)

    def get_contribution(self):
        """Get source contribution dofnos and value

        See documentation of calc_pointsource_contrib for more detail
        on the return values
        """
        return calc_pointsource_contrib(self.function_space, self.position, self.value)

def calc_pointsource_contrib(V, source_coords, source_value):
    """Calculate the RHS contribution of a current point source (i.e. electric dipole)
    
    @param V: dolfin FunctionSpace object
    @param source_coords: length 3 array with x,y,z coordinates of point source
    @param source_value: length 3 array with x,y,z componets of source current 
    
    @rtype: (C{numpy.array}, C{numpy.array})
    @return: (dofnos, rhs_contribs) -- An array containing the indices of the degrees of freedom associated with
        the source, and the numerical values of the contributions of the current source.
        
        C{RHS[dofnos] += rhs_contribs} will add the current source to the system's RHS.
    """
    source_coords = N.asarray(source_coords, dtype=N.float64)
    source_value = N.asarray(source_value)
    dm = V.dofmap()
    dofnos = N.zeros(dm.max_cell_dimension(), dtype=N.uintc)
    source_pt = dolfin.Point(*source_coords)
    try:
        cell_index = V.mesh().any_intersected_entity(source_pt)
    except StandardError:
        # CGAL as used by dolfin to implement intersection searches
        # seems to break with 1-element meshes
        if dolfin.Cell(V.mesh(), 0).intersects(source_pt):
            cell_index = 0
        else: raise
    c = dolfin.Cell(V.mesh(), cell_index)
    # Check that the source point is in this element    
    assert(c.intersects_exactly(source_pt)) 

    dm.tabulate_dofs(dofnos,  c)
    finite_element = V.dolfin_element()
    no_basis_fns = finite_element.space_dimension()
    # Vector valued elements have rank of 1
    assert(finite_element.value_rank() == 1)
    # Vector valued elements have only one rank (i.e. 0) along which
    # dimensions are defined. This is the dimension that the basis
    # function value vector is. Since we have 3D Nedelec elements here
    # this should be 3
    bf_value_dimension = finite_element.value_dimension(0)
    el_basis_vals = N.zeros((no_basis_fns, bf_value_dimension), dtype=N.float64)
    finite_element.evaluate_basis_all(el_basis_vals, source_coords, c)
    rhs_contribs = N.sum(el_basis_vals*source_value, axis=1)
    return dofnos, rhs_contribs

