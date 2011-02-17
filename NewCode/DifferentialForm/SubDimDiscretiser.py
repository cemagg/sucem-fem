from __future__ import division

import numpy as N
from NewCode import ProxyList, Integration
from NewCode.Utilities import Struct, memoized2
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiserElement

class BasePformSubDimDiscretiser(object):
    defaultIntegrationOrder = 8
    elementEntity = 'face'
    Integrator = None
    def __init__(self, p, mesh, geomEntities, permuter, superDisc, *names, **kwargs):
        self.superMesh = mesh.superMesh
        self.superDisc = superDisc
        super(BasePformSubDimDiscretiser, self).__init__(
            p, mesh, geomEntities, permuter, *names, **kwargs)
        self.elements.entity.setBasisFunctions(superDisc.elements.entity)
        self.initPermutation(permuter, self.elements.entity)

    def _setElements(self, freefun=None):
        try: ElClass=self.discElClasses[self.p]
        except KeyError:
            raise NotImplementedError, "%d-form not fully implemented" % self.p
        self.elements = ProxyList.ProxyList(ElClass(
            self.mesh.elements.list_repr(), superMesh=self.superMesh,
            mesh=self.mesh, freefun=freefun))
        # A hack, assming that subdimg faces are always free FIXME
        self.elements.noFree = len(self.elements)

class GenericPformSubDimDiscretiser(BasePformSubDimDiscretiser):
    @memoized2
    def D(self):
        superDisc_D = self.superDisc.D()
        return self.DClass(self.p, self.mesh, self.geomEntities,
                           self.permuter, superDisc_D)

class GenericPformSubDimDiscretiser_D(BasePformSubDimDiscretiser):
    def D(self):
        raise Exception("Repeated Application of exterior derivative gives zero")


class PformSubDimDiscretiser_D(GenericPformSubDimDiscretiser_D,
                               Discretiser.BasePformDiscretiser): 
    discElClasses = {1:SubDimDiscretiserElement.OneformSubDimElement_D}
    Integrator = Integration.TriIntegrator
    from DiscretiserDOFs import PformDiscretiserDOFs_D as DiscDOFsClass
    
class PformSubDimDiscretiser(GenericPformSubDimDiscretiser,
                             Discretiser.BasePformDiscretiser): 
    DClass = PformSubDimDiscretiser_D
    discElClasses = {1:SubDimDiscretiserElement.OneformSubDimElement}
    Integrator = Integration.TriIntegrator
    from DiscretiserDOFs import PformDiscretiserDOFs as DiscDOFsClass

class GenericPformOutwardSubDimDiscretiser(GenericPformSubDimDiscretiser):
    def set_interior_selector(self, insidep):
        super_els = self.elements[:].connect2SuperElems
        super_localfacenos = self.elements[:].superLocalFacenos
        el_inside = N.array([insidep(N.average(el.nodeCoords, axis=0))
                             for el in self.superMesh.elements[
            super_els[:,0]]], N.bool)
        
        self.elements.entity._superElement = N.where(
            el_inside, super_els[:,0], super_els[:,1])
        self.elements.entity._superLocalFaceno = N.where(
            el_inside, super_localfacenos[:,0], super_localfacenos[:,1])

class PformOutwardSubDimDiscretiser(GenericPformOutwardSubDimDiscretiser,
                                    Discretiser.BasePformDiscretiser): 
    DClass = None
    discElClasses = {1:SubDimDiscretiserElement.OutwardOneformSubDimElement}
    Integrator = Integration.TriIntegrator
    from DiscretiserDOFs import PformDiscretiserDOFs as DiscDOFsClass
