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
        """Calculate the RHS contribution of all current contained sources 

        Return Values
        -------------
        (dofnos, rhs_contribs) with

        dofnos -- Array of degree of freedom indices of the source contribution

        rhs_contribs -- Numerical values of RHS contribution

        Calculated such that RHS[dofnos] += rhs_contribs will add the
        current sources' contribution to the system.

        """
        contribs = collections.defaultdict(lambda : 0.)
        for src in self.sources:
            dofnos, values = src.get_contribution()
            for dn, v in itertools.izip(dofnos, values):
                contribs[dn] = contribs[dn] + v

        dofnos = N.array(contribs.keys(), dtype=N.uint)
        rhs_contribs = N.array(contribs.values())
        return dofnos, rhs_contribs
