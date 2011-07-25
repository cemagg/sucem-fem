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
from __future__ import division

from sucemfem.Sources.fillament_current_source import FillamentCurrentSource
from sucemfem.Utilities.Geometry import unit_vector, vector_length

class FillamentSource(object):
    """High level current fillament source representation"""
    def __init__(self, function_space):
        self.function_space = function_space

    def set_source_parameters(self, source_parameters):
        """Set current fillament source paerameters using dict

        source_parameters is a dict with keys:

        I -- Source current in amperes
        endpoints -- 2x3 array with start and end coordinates of the source
        """
        self.source_parameters = source_parameters
        endpoints = source_parameters['endpoints']
        delta = endpoints[1] - endpoints[0]
        self.direction = unit_vector(delta)
        self.length = vector_length(delta)

    def get_current_source(self):
        """get FillamentCurrentSource instance

        Each call results in a new instance being instantiated.
        """
        cs = FillamentCurrentSource()
        cs.set_function_space(self.function_space)
        cs.set_source_endpoints(self.source_parameters['endpoints'])
        cs.set_value(self.source_parameters['I'])
        return cs
    
