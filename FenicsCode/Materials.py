from __future__ import division

import numpy 
from NewCode import Consts
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
    """Set up a material region number to MaterialProperties lookup

    Input Parameters
    ----------------

    material_regions -- Dictionary of with keys equal to material
        region number, and sub-dictionaries of parameter name and numeric
        values. E.g.

        {100:{'eps_r':1, 'mu_r':1}, 101:{'eps_r':2, 'mu_r':1}}

        Allowed material parameter names are defined by the class
        MaterialProperties. Default values for unspecified materials
        are also handled by that class
    
    """
    def __init__(self, material_regions):
        self.material_properties = {}
        for region_no, property_values in material_regions.items():
            matprop = MaterialProperties()
            matprop.init_values(**property_values)
            self.material_properties[region_no] = matprop

    def get_material_properties(self):
        return self.material_properties

class MaterialFunctionFactory(object):
    """Set up material property functions

    Input Parameters
    ----------------

    region_material_properties -- dict with key region_no, value
        MaterialProperties object for that region number

    region_meshfunction -- a dolfin MeshFunction mapping elements to
        region numbers

    mesh -- dolfin Mesh object relating to region_meshfunction. The
        material functions will be defined on this mesh

    """

    def __init__(self, region_material_properties, region_meshfunction, mesh):
        self.region_material_properties = region_material_properties
        self.region_meshfunction = region_meshfunction
        self.mesh = mesh

    def get_material_functions(self, *property_names):
        """Return material functions for the requested properties

        Input Parameters
        ----------------

        property_names -- Sequence of property names for which to
            calculate material functions
            
        Note, the property names have to be defined by the
        region_material_properties objects passed to the class
        constructor

        Return Value

        dict with {property_name:property_value_function}

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
            matfn.vector()[:] = numpy.array([region_valuemap[int(i)]
                                             for i in region_meshfn.array()])
            mat_fns[pname] = matfn
        return mat_fns
        
