from NewCode import BrickSubDimMesh
from NewCode.DifferentialForm import SubDimDiscretiserElement
from NewCode.DifferentialForm import BrickDiscretiserElement

class PformSubDimElement(SubDimDiscretiserElement.BasePformSubDimElement,
                         BrickDiscretiserElement.NdimPformElement,
                         BrickSubDimMesh.FaceElement):
    SubDimMeshModule = BrickSubDimMesh
    facenos = property(BrickDiscretiserElement.NdimPformElement.entityNo)

class OneformSubDimElement(SubDimDiscretiserElement.BaseOneformSubDimElement,
                           PformSubDimElement): pass

class OutwardOneformSubDimElement(
    SubDimDiscretiserElement.BaseOneformOutwardSubDimElement,
    PformSubDimElement): pass

class OneformSubDimElement_D(
    SubDimDiscretiserElement.OneformSubDimElement_DPhysvals,
    OneformSubDimElement): pass

