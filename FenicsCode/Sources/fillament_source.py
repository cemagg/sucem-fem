from __future__ import division

from FenicsCode.Sources.fillament_current_source import FillamentCurrentSource
from FenicsCode.Utilities.Geometry import unit_vector, vector_length

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
    
