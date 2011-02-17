from __future__ import division

from NewCode.DifferentialForm import Discretiser, PyramDiscretiserEntities
from NewCode.DifferentialForm import PyramDiscretiserElement
from NewCode.DifferentialForm.Permuter import Permuter
from Discretiser import allfree

class BasePyramDiscretiser(object):
    from NewCode.Integration import PyramIntegrator as Integrator
    defaultIntegrationOrder = (3,4)         # Actually no of intg pts

class PformDiscretiser_D(BasePyramDiscretiser,
                         Discretiser.PformDiscretiser_D):
    discElClasses = {1:PyramDiscretiserElement.OneformElement_D,
                     }
    

class PformDiscretiser(BasePyramDiscretiser,
                       Discretiser.PformDiscretiser):
    Discretiser_D = PformDiscretiser_D
    discElClasses = {1:PyramDiscretiserElement.OneformElement,
                     }

def setup_PformDiscretiser(
    mesh, form, order=1, mixed=True,
    freeFun=allfree, vol_freeFun=None, btype='graglia99'):
    """
    Set up a discretiser, using defaults for basis function and Permuter selection

    Frees the user from having to setup the basisSet and DiscretiserEntities
    """
    from NewCode.DifferentialForm.BasisFunction import pyramBasisForm
    basisSet = pyramBasisForm[form].basis_set(order=order, mixed=mixed, btype=btype)
    geomEntities = PyramDiscretiserEntities.make_geomEntities(
        mesh, basisSet, freeFun, vol_freeFun)
    return  PformDiscretiser(
        form, mesh, geomEntities, Permuter, basisSet)
