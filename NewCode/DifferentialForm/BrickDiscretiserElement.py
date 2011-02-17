from __future__ import division

import numpy as N

from NewCode.DifferentialForm import DiscretiserElement
from NewCode.Meshes import BrickMesh

NdimPformElement = DiscretiserElement.NdimPformElement

class PformElement(DiscretiserElement.BasePformElement,
                   BrickMesh.Element): 
    def _get_local_numbering(self):
        return dict(edge=(self.FACEDUAL_EDGE_NUMBERING,),
                    face=(self.FACEDUAL_FACE_NUMBERING,))

    def refVals(self):
        """
        Evaluate all basis functions on the reference element at all integration points

        Return Value
        ============

        Assuming n basisfunctions, bf1, ... bfn, a list is returned with one
        element per bf, containing the given bf evaluated at each integration
        point pt1, .., ptm eg:

        [[bf1(pt1), ... bf1(ptm)], ... [bfn(pt1), ... bfn(ptm)]]

        If the integration rule requires different points for the x/y/z
        directed basisfunctions, the appropriate pointset is chosed by quering
        basisfunction.direction
        
        """
        
        try: return self._refVals
        except AttributeError: pass
        
        intg_points = self.rule.evalPoints()
        # Call standard refVals routine for standard integration rules
        if not hasattr(intg_points, 'xx'):
            return DiscretiserElement.BasePformElement.refVals(self)
        coord = {(1.,0.,0.):0,
                 (0.,1.,0.):1,
                 (0.,0.,1.):2}
        vals_xx = self.refValsAtPoints(intg_points.xx)
        vals_yy = self.refValsAtPoints(intg_points.yy)
        vals_zz = self.refValsAtPoints(intg_points.zz)
        vals = (vals_xx, vals_yy, vals_zz)
        
        self._refVals = N.array([vals[coord[tuple(bf.direction)]][i]
                                 for i, bf in enumerate(self.basisSet)],
                                N.float64)

        return self._refVals

    def alt_refVals(self):
        try: return self._alt_refVals
        except AttributeError: pass
        
        intg_points = self.alt_rule.evalPoints()
        coord = {(1.,0.,0.):0,
                 (0.,1.,0.):1,
                 (0.,0.,1.):2}
        vals_xx = self.refValsAtPoints(intg_points.xx)
        vals_yy = self.refValsAtPoints(intg_points.yy)
        vals_zz = self.refValsAtPoints(intg_points.zz)
        self._alt_refVals = (vals_xx, vals_yy, vals_zz)
        
        return self._alt_refVals

    def alt_physVals(self):
        return N.array(zip(*[self.physValsFromRefVals(rvs)
                             for rvs in self.alt_refVals()]))

    def set_alt_IntegrationRule(self, intg_rule):
        self.alt_rule = intg_rule
        try: del self._alt_refVals               # Kill any cached values
        except AttributeError: pass


class OneformElement(PformElement):
    entityNames = PformElement.entityNames[1:]
    coordVecs = PformElement.covBaseVecs
    
class OneformElement_D(OneformElement):
    coordVecs = PformElement.conBaseVecs

class TwoformElement(PformElement):
    entityNames = PformElement.entityNames[2:]
    coordVecs = PformElement.conBaseVecs

class TwoformElement_D(TwoformElement):
    def physValsFromRefVals (self, refVals):
        return N.array(refVals, N.float64)/N.linalg.det(self.J())
