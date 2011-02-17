from NewCode.Utilities import partial, Struct
import numpy as N
from numpy.linalg import det, inv
from coord_vecs import conv_coord_vecs
import Oneform

zero_vec1 = N.zeros(1, N.float64)
# Make this set of zeros immutable so only explicit rebinding can change it
zero_vec1.flags.writeable = False
def zero_D(*names):
    def zero_D_(coords):
        return zero_vec1
    return zero_D_

def facefun0(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    def facefun0_(coord):
        try:
            coords = coord[local_facenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')    
        # 2*(lambda_A*(grad(lambda_B) x grad(lambda_C)) +
        #    lambda_B*(grad(lambda_C) x grad(lambda_A)) +
        #    lambda_C*(grad(lambda_A) x grad(lambda_B)))
        return 2*N.sum(coords[i]*conv_coord_vecs[local_facenodes[(i+1)%3],
                                                 local_facenodes[(i+2)%3]]
                       for i in range(3))
    return facefun0_

def facefun0_D(local_faceno, local_facenodes):
    def facefun0_D_(coord):
        return N.array([6*(-1)**(local_faceno % 2)], N.float64)
    return facefun0_D_

def basis_set(order, mixed=True):
# Perhaps we should look to adapting the one from Oneform.py
    assert(order > 0)
    if (order > 1 and mixed) or order > 3:
        raise NotImplementedError
    face_toget = order + (1 if not (order==1 and mixed) else 0)
    facefuns, facefuns_D = [reduce(lambda x,y: x+y, all_funs[0:face_toget])
                            for all_funs in (all_facefuns, all_facefuns_D)]
    order = order - (mixed and 0.5 or 0)
    fns=dict(face=facefuns); fns_D=dict(face=facefuns_D)
    if order >=2:
        vol_toget_R = order
        volfuns, volfuns_D = ([], [])
        for i in range(2, vol_toget_R+1):
            volfuns.extend(all_volfunsR[i] )       
            volfuns_D.extend(all_volfunsR_D[i])
        fns['vol'] = volfuns ; fns_D['vol'] = volfuns_D
    # ..._solh since all the higher order functions are solenoidal.
    return Struct(fns=fns, fns_D=fns_D,
                  info=Struct(form=2, type='mmbotha06_solh', order=order))

facefuns0 = tuple([partial(facefun0, i) for i in range(4)])
facefuns0_D = tuple([partial(facefun0_D, i) for i in range(4)])

facefunsR1 = Oneform.facefunsR2_D
facefunsR1_D = (zero_D,)*8

facefunsR2 = Oneform.facefunsR3_D
facefunsR2_D = (zero_D,)*len(facefunsR2)

facefunsR3 = Oneform.facefunsR4_D
facefunsR3_D = (zero_D,)*len(facefunsR3)

volfunsR2 = Oneform.volfunsR3_D
volfunsR2_D = (zero_D,)*len(volfunsR2)

volfunsR3 = Oneform.volfunsR4_D
volfunsR3_D = (zero_D,)*len(volfunsR3)


all_facefuns = [facefuns0, facefunsR1, facefunsR2, facefunsR3]
all_facefuns_D = [facefuns0_D, facefunsR1_D, facefunsR2_D, facefunsR3_D]

all_volfunsR = {
    2:volfunsR2,
    3:volfunsR3,
    }
all_volfunsR_D = {
    2:volfunsR2_D,
    3:volfunsR3_D,
                  }
