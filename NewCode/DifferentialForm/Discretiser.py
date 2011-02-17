from __future__ import division

import new
import numpy as N
from NewCode import ProxyList
from NewCode.DifferentialForm.DiscretiserMatrices import DiscretiserMatrices
from NewCode.Utilities import Struct, isBound, memoized2
from NewCode.DifferentialForm import DiscretiserEntities
from NewCode.DifferentialForm.Permuter import Permuter
from DiscretiserEntities import geom_entity_names
from NewCode.DifferentialForm import DiscretiserElement
from BoundaryConditions import allfree

class BasePformDiscretiser(object):
    discElClasses = {1:DiscretiserElement.OneformElement,
                     2:DiscretiserElement.TwoformElement}
    from NewCode.Integration import TetIntegrator as Integrator
    defaultIntegrationOrder = 8
    elementEntity = 'vol'
    identicalElements = False           # is each element is locally identical?
    diagonalMetric = False              # Is the unit metric (mass) matrix diag?
    
    def __init__(self, p, mesh, geomEntities, *names, **kwargs):
        self.p = p
        self.mesh = mesh
        try:
            el_ent = geomEntities[self.elementEntity]
            freefun = el_ent.entity.freefun
        except KeyError: freefun=None
        self._setElements(freefun=freefun)
        self.setIntegrationRule(self.defaultIntegrationOrder)
        if self.elementEntity in geomEntities:
            geomEntities[self.elementEntity] = self.elements
            self.elements.entity.entity_type = self.elementEntity
            self.elements.mesh = mesh

        # Check that the same mesh geomEntities' meshes are consistent
        for ent in geomEntities.values():
            assert(mesh is ent.mesh)
        self.geomEntities = geomEntities
        self.matrix = DiscretiserMatrices(self)
        
        super(BasePformDiscretiser, self).__init__(p, mesh, geomEntities,
                                                   *names, **kwargs)

    def _setElements(self, freefun=None):
        try: ElClass=self.discElClasses[self.p]
        except KeyError:
            raise NotImplementedError, "%d-form not fully implemented" % self.p
        proxy_el = ElClass(attrs=self.mesh.elements.list_repr(partial=True),
                           mesh=self.mesh, freefun=freefun)
        self.elements = DiscretiserEntities.DiscretiserEntityListNoBoundary(
            proxy_el)

    def setIntegrationRule(self, integrationOrder):
        self.elements.entity.setIntegrationRule(self.Integrator(integrationOrder))

    def set_integrator(self, integrator):
        self.Integrator = integrator

    def setFaceIntegrationRule(self, rule):
        self.elements.entity.setFaceIntegrationRule(rule)

    def initPermutation(self, permuter, parent_element, *names, **kwargs):
        self.permuter = isBound(permuter.__init__) and permuter \
                        or permuter(parent_element, self.geomEntities)
        parent_element.permutation = new.instancemethod(
            self.permuter.permuteElement, parent_element, self.permuter)
        self.newDOFsArray = self.permuter.newDOFs
        self.totalDOFs = self.permuter.totalDOFs

    def newDOFs(self, dtype=N.float64):
        return self.DiscDOFsClass(self, dtype)        
        
class PformDiscretiserSetsElementBasis(BasePformDiscretiser):
    def initPermutation(self, permuter, parent_element, basisFuns, *names, **kwargs):
        parent_element.setBasisFunctions(basisFuns)
        super(PformDiscretiserSetsElementBasis, self).initPermutation(
            permuter, parent_element, *names, **kwargs)

class PformDiscretiser_D(PformDiscretiserSetsElementBasis):
    discElClasses={1: DiscretiserElement.OneformElement_D,
                   2: DiscretiserElement.TwoformElement_D}
    from DiscretiserDOFs import PformDiscretiserDOFs_D as DiscDOFsClass
                     
    def __init__(self, p, mesh, geomEntities, permuter, basisSet, *names, **kwargs):
        super(PformDiscretiser_D, self).__init__(p, mesh, geomEntities,
                                                 *names, **kwargs)
        self.initPermutation(permuter, self.elements.entity, basisSet.fns_D)


class PformDiscretiser(PformDiscretiserSetsElementBasis):
    Discretiser_D = PformDiscretiser_D
    from DiscretiserDOFs import PformDiscretiserDOFs as DiscDOFsClass

    def __init__(self, p, mesh, geomEntities, Permuter, basisSet, *names, **kwargs):
        """

        """
        self.basisSet = basisSet
        super(PformDiscretiser, self).__init__(p, mesh, geomEntities,
                                                 *names, **kwargs)
        self.initPermutation(Permuter, self.elements.entity, basisSet.fns)

    @memoized2
    def D(self):
        """
        Perform exterior derivative on basis functions and return equivalent discriteser
        """
        # Check that the basisfunctions support the exterior derivative
        assert(hasattr(self.basisSet, 'fns_D'))
        disc_D =  self.Discretiser_D(self.p, self.mesh, self.geomEntities,
                                     self.permuter, self.basisSet)
        disc_D.elements.entity.dead_elements = self.elements.entity.dead_elements
        return disc_D

def setup_PformDiscretiser(
    mesh, form, order=1, mixed=True, freeFun=allfree, vol_freeFun=None,
    el_in_disc=None, btype=None):
    """
    Set up a discretiser, using defaults for basis function and Permuter selection

    Frees the user from having to setup the basisSet and DiscretiserEntities
    """
    from NewCode.DifferentialForm.BasisFunction import basisForm
    basisSet = basisForm[form].basis_set(order=order, mixed=mixed)
    geomEntities = DiscretiserEntities.make_geomEntities(
        mesh, basisSet, freeFun, vol_freeFun, el_in_disc)
    return  PformDiscretiser(
        form, mesh, geomEntities, Permuter, basisSet)


