# Authors
# Neilen Marais <nmarais@gmail.com>
from __future__ import division
import collections
import itertools
import numpy as N


class CurrentSources(object):
    def __init__(self):
        self.sources = []

    def set_function_space(self, function_space):
        self.function_space = function_space

    def add_source(self, source):
        self.sources.append(source)

    def init_sources(self):
        for src in self.sources:
            src.set_function_space(self.function_space)
        
    def get_source_contributions(self):
        """Get and return the RHS contribution of the current source
        """
        contribs = collections.defaultdict(lambda : 0.)
        for src in self.sources:
            dofnos, values = src.get_contribution()
            for dn, v in itertools.izip(dofnos, values):
                contribs[dn] = contribs[dn] + v

        dofnos = N.array(contribs.keys(), dtype=N.uint)
        rhs_contribs = N.array(contribs.values())
        return dofnos, rhs_contribs

class CurrentSource(object):
    """Abstract base class for current sources"""

    def set_function_space(self, function_space):
        """Set function space that the source is to be applied to"""
        self.function_space = function_space

    def get_contribution(self):
        """Get and return the RHS contribution of the current source
        
        @raise NotImplementedError: This method should be implemented in a sub-class.
        @rtype: (C{numpy.array}, C{numpy.array})
        @return: (dofnos, rhs_contribs) -- An array containing the indices of the degrees of freedom associated with
            the source, and the numerical values of the contributions of the current source.
            
            C{RHS[dofnos] += rhs_contribs} will add the current source to the system's RHS.
        """
        raise NotImplementedError('User subclass should implement get_contribution()')
