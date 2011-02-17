from __future__ import division

import numpy as N

from NewCode.DifferentialForm import DiscretiserElement
from NewCode.Meshes import PyramMesh

NdimPformElement = DiscretiserElement.NdimPformElement

class PformElement(DiscretiserElement.BasePformElement,
                   PyramMesh.Element): 
    entityNames = ('node', 'edge', 'baseface', 'apexface', 'vol')
    def _get_local_numbering(self):
        edge_info = (self.FACEDUAL_EDGE_NUMBERING, self.FACEDUAL_EDGE_SENSE)
        return dict(edge=edge_info,
                    apexface=(self.LOCAL_APEXFACEEDGES,)+edge_info,
                    baseface=edge_info, vol=edge_info)


class OneformElement(PformElement):
    entityNames = PformElement.entityNames[1:]
    coordVecs = PformElement.covBaseVecs
    
class OneformElement_D(OneformElement):
    coordVecs = PformElement.conBaseVecs

    
