from __future__ import division

from NewCode import Integration
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser, DiscretiserDOFs
from NewCode.DifferentialForm import BrickSubDimDiscretiserElement


class PformSubDimDiscretiser_D(SubDimDiscretiser.GenericPformSubDimDiscretiser_D,
                               Discretiser.BasePformDiscretiser): 
    discElClasses = {1:BrickSubDimDiscretiserElement.OneformSubDimElement_D}
    Integrator = Integration.QuadIntegrator
    from DiscretiserDOFs import PformDiscretiserDOFs_D as DiscDOFsClass
    
class PformSubDimDiscretiser(SubDimDiscretiser.GenericPformSubDimDiscretiser,
                             Discretiser.BasePformDiscretiser): 
    discElClasses = {1:BrickSubDimDiscretiserElement.OneformSubDimElement}
    Integrator = Integration.QuadIntegrator
    DClass = PformSubDimDiscretiser_D
    from DiscretiserDOFs import PformDiscretiserDOFs as DiscDOFsClass    

class PformOutwardSubDimDiscretiser(SubDimDiscretiser.GenericPformOutwardSubDimDiscretiser,
                             Discretiser.BasePformDiscretiser): 
    discElClasses = {1:BrickSubDimDiscretiserElement.OutwardOneformSubDimElement}
    Integrator = Integration.QuadIntegrator
    DClass = PformSubDimDiscretiser_D
    from DiscretiserDOFs import PformDiscretiserDOFs as DiscDOFsClass    
