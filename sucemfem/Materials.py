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

import numpy 
from sucemfem import Consts
import dolfin

mu0 = Consts.mu0
eps0 = Consts.eps0

class MaterialProperties(object):
    allowed_properties_defaults = {'eps_r':1.0, 'mu_r':1.0}

    def __init__(self):
        for prop, val in self.allowed_properties_defaults.items():
            setattr(self, prop, val)

    def init_values(self, **kwargs):
        for prop, val in kwargs.items():
            if prop not in self.allowed_properties_defaults:
                raise ValueError('Unknown material property')
            setattr(self, prop, val)
            
    def set_eps_r(self, eps_r):
        self.eps_r = eps_r

    def set_mu_r(self, mu_r):
        self.mu_r = mu_r

    def get_eps_r(self):
        return self.eps_r

    def get_mu_r(self):
        return self.mu_r

    def get_mu_r_inv(self):
        return 1/self.get_mu_r()

    def get_eps(self):
        return eps0*self.eps_r

    def get_mu(self):
        return mu0*self.mu_r

class MaterialPropertiesFactory(object):
    def __init__(self, material_regions):
        """Set up a material region number to MaterialProperties lookup

        @param material_regions: A dictionary of with keys equal to material
            region number, and sub-dictionaries of parameter name and numeric
            values. E.g.
    
            {100:{'eps_r':1, 'mu_r':1}, 101:{'eps_r':2, 'mu_r':1}}
    
            Allowed material parameter names are defined by the class
            MaterialProperties. Default values for unspecified materials
            are also handled by that class
    
            If material_regions is None, the whole domain defaults to
            region 0, which is typicaly freespace
        
        """

        self.material_properties = {}
        if material_regions is None:
            # default to default material properties in region 0
            matprop = MaterialProperties()
            self.material_properties[0] = matprop
        else:
            for region_no, property_values in material_regions.items():
                matprop = MaterialProperties()
                matprop.init_values(**property_values)
                self.material_properties[region_no] = matprop

    def get_material_properties(self):
        return self.material_properties

class MaterialFunctionFactory(object):

    def __init__(self, region_material_properties, region_meshfunction, mesh):
        """Set up material property functions
    
        @param region_material_properties: A dict with key region_no, value
            MaterialProperties object for that region number
        @param region_meshfunction: A dolfin MeshFunction mapping elements to
            region numbers. If it is None, all elements are set to region 0
        @param mesh: A dolfin Mesh object relating to region_meshfunction. The
            material functions will be defined on this mesh
        """
        self.region_material_properties = region_material_properties
        # if the meshfunction is not defined then initialise to a zero mesh function
        if region_meshfunction is None:
            region_meshfunction = dolfin.MeshFunction (
                'uint', mesh, mesh.topology().dim() )
            region_meshfunction.set_all (0)
        self.region_meshfunction = region_meshfunction
        self.mesh = mesh

    def get_material_functions(self, *property_names):
        """Return material functions for the requested properties

        @param property_names: A sequence of property names for which to
            calculate material functions
        
            Note, the property names have to be defined by the
            region_material_properties objects passed to the class
            constructor

        @return: a dict with elements {property_name:property_value_function}
        """
        
        mat_fns = {}
        # Lowest order discontinuous Galerkin space gives you
        # element-wise constant functions, which is how material
        # properties are defined.
        mat_funcspace = dolfin.FunctionSpace(self.mesh, 'DG', 0)
        region_meshfn = self.region_meshfunction
        mats = self.region_material_properties
        for pname in property_names:
            valfunc = 'get_'+pname
            # Simple dictionary of region number to value for property pname
            region_valuemap = dict((k, getattr(mat, valfunc)())
                                   for k, mat in mats.items())
            matfn = dolfin.Function(mat_funcspace)
            try:
                matfn.vector()[:] = numpy.array([region_valuemap[int(i)]
                                                 for i in region_meshfn.array()])
            except KeyError, s:
                raise ValueError('Material number %d not found' % s.args[0])
            mat_fns[pname] = matfn
        return mat_fns
        
