from __future__ import division

import dolfin
import numpy as N

class CurrentSources(object):
    def __init__(self):
        self.sources = []

    def set_function_space(self, function_space):
        self.function_space = function_space

    def add_source(self, source):
        source.set_function_space(self.function_space)
        self.sources.append(source)
        
    def get_source_contributions(self):
        pass

class PointCurrentSource(object):
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

    def set_function_space(self, function_space):
        """Set function space that the source is to be applied to"""
        self.function_space = function_space

    def get_contribution(self):
        """Get source contribution dofnos and value

        See documentation of calc_pointsource_contrib for more detail
        on the return values
        """
        return calc_pointsource_contrib(self.function_space, self.position, self.value)

def calc_pointsource_contrib(V, source_coords, source_value):
    """Calculate the RHS contribution of a current point source (i.e. electric dipole)
    Input Values
    -------------
    @param V: dolfin FunctionSpace object
    @param source_coords: length 3 array with x,y,z coordinates of point source
    @param source_value: length 3 array with x,y,z componets of source current 
    Return Values
    -------------
    (dofnos, rhs_contribs) with

    dofnos -- Array of degree of freedom indices of the source contribution

    rhs_contribs -- Numerical values of RHS contribution, such that
        RHS[dofnos] += rhs_contribs will add the current source to the system.

    """
    source_coords = N.asarray(source_coords, dtype=N.float64)
    source_value = N.asarray(source_value)
    dm = V.dofmap()
    dofnos = N.zeros(dm.max_cell_dimension(), dtype=N.uint32)
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

