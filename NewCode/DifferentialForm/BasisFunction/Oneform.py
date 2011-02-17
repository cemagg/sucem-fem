import numpy as N

from NewCode.Utilities import partial, Struct
from coord_vecs import cov_coord_vecs, conv_coord_vecs

zero_vec6 = N.zeros(6, N.float64)
# Make this set of zeros immutable so only explicit rebinding can change it
zero_vec6.flags.writeable = False
def zero_D(*names):
    def zero_D_(coords):
        return zero_vec6
    return zero_D_

def edgefun0(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (grad0, grad1) = cov_coord_vecs[local_edgenodes]
    def edgefun0_(coord):
        try:
            coords = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return coords[0]*grad1 - coords[1]*grad0
    return edgefun0_

def edgefun0_D(local_edgeno, local_edgenodes):
    l_edge0, l_edge1 = local_edgenodes = local_edgenodes[local_edgeno] 
    def edgefun0_D_(coord):
        try:
            coords = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return 2*conv_coord_vecs[l_edge0, l_edge1]
    return edgefun0_D_

def edgefun1(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (grad0, grad1) = cov_coord_vecs[local_edgenodes]
    def edgefun1_(coord):
        try:
            coords = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return coords[0]*grad1 + coords[1]*grad0
    return edgefun1_

def edgefun2(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (grad1, grad2) = cov_coord_vecs[local_edgenodes]
    def edgefun2_(coord):
        try:
            (lam1, lam2) = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return  -lam2*(lam2-2*lam1)*grad1 - lam1*(2*lam2-lam1)*grad2
    return edgefun2_

def edgefun3(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (gv1, gv2) = cov_coord_vecs[local_edgenodes]
    def edgefun3_(coord):
        try:
            (l1, l2) = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return gv2*l1*(l2-l1)*(3*l2-l1) + gv1*l2*(l2-3*l1)*(l2-l1)
    return edgefun3_

def edgefun4(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (gv1, gv2) = cov_coord_vecs[local_edgenodes]
    def edgefun4_(coord):
        try:
            (l1, l2) = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return -gv2*l1*(l2-l1)**2*(4*l2-l1) - gv1*l2*(l2-4*l1)*(l2-l1)**2
    return edgefun4_

def edgefun5(local_edgeno, local_edgenodes):
    local_edgenodes = local_edgenodes[local_edgeno]
    (gv1, gv2) = cov_coord_vecs[local_edgenodes]
    def edgefun5_(coord):
        try:
            (l1, l2) = coord[local_edgenodes]
        except TypeError:
            raise TypeError('coord must be a numpy array')
        return gv2*l1*(l2-l1)**3*(5*l2-l1)+gv1*l2*(l2-5*l1)*(l2-l1)**3
    return edgefun5_

def facefunR20_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    grads = cov_coord_vecs[local_facenodes]
    def facefunR20_1_(coord):
        try: coords = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return coords[1]*coords[2]*grads[0] + coords[0]*coords[2]*grads[1] \
               - 2*coords[0]*coords[1]*grads[2]
    return facefunR20_1_

def facefunR20_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    grads = cov_coord_vecs[local_facenodes]
    def facefunR20_2_(coord):
        try: coords = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -2*coords[1]*coords[2]*grads[0] + coords[0]*coords[2]*grads[1] \
               + coords[0]*coords[1]*grads[2]
    return facefunR20_2_

def facefunR30_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR30_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return l1*l2*(l1-l2)*gv3
    return facefunR30_1_

def facefunR30_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR30_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return l2*l3*(l2-l3)*gv1
    return facefunR30_2_

def facefunR30_3(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR30_3_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return l1*l3*(l3-l1)*gv2
    return facefunR30_3_

def facefunR40_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR40_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv3*l1*l2*(l2-l1)**2
    return facefunR40_1_

def facefunR40_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR40_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv1*l2*l3*(l3-l2)**2
    return facefunR40_2_

def facefunR40_3(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR40_3_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv2*l1*l3*(l3-l1)**2
    return facefunR40_3_

def facefunR41_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR41_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -gv1*l1*l2*l3*(l3-l2)+gv2*l1*l2*l3*(l3-l1)-gv3*l1*l2*(l2-l1)*l3
    return facefunR41_1_

def facefunR50_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR50_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -gv3*l1*l2*(l2-l1)**3
    return facefunR50_1_

def facefunR50_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR50_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -gv1*l2*l3*(l3-l2)**3
    return facefunR50_2_

def facefunR50_3(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR50_3_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv2*l1*l3*(l3-l1)**3
    return facefunR50_3_

def facefunR51_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR51_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv1*l1*l2**2*l3**2+gv2*l1**2*l2*l3**2-2*gv3*l1**2*l2**2*l3
    return facefunR51_1_

def facefunR51_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunR51_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -2*gv1*l1*l2**2*l3**2+gv2*l1**2*l2*l3**2+gv3*l1**2*l2**2*l3
    return facefunR51_2_


def facefunG20_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG20_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return l2*l3*gv1 + l1*l3*gv2 + l1*l2*gv3
    return facefunG20_1_

def facefunG30_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG30_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return  -gv2*l1*(2*l2-l1)*l3-gv1*l2*(l2-2*l1)*l3-gv3*l1*l2*(l2-l1)
    return facefunG30_1_

def facefunG30_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG30_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return   -gv3*l1*l2*(2*l3-l2)-gv1*l2*l3*(l3-l2)-gv2*l1*l3*(l3-2*l2)
    return facefunG30_2_

def facefunG40_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG40_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return 2*gv1*l1*l2**2*l3+2*gv2*l1**2*l2*l3+gv3*l1**2*l2**2
    return facefunG40_1_

def facefunG40_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG40_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv1*l2**2*l3**2+2*gv2*l1*l2*l3**2+2*gv3*l1*l2**2*l3
    return facefunG40_2_

def facefunG40_3(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG40_3_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return 2*gv1*l1*l2*l3**2+gv2*l1**2*l3**2+2*gv3*l1**2*l2*l3
    return facefunG40_3_

def facefunG50_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG50_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -gv2*l1**2*l2*(3*l2-2*l1)*l3-gv1*l1*l2**2*(2*l2-3*l1)*l3 \
               -gv3*l1**2*l2**2*(l2-l1)
    return facefunG50_1_

def facefunG50_2(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG50_2_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -gv3*l1*l2**2*l3*(3*l3-2*l2)-gv2*l1*l2*l3**2*(2*l3-3*l2) \
               -gv1*l2**2*l3**2*(l3-l2)
    return facefunG50_2_

def facefunG50_3(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG50_3_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return gv3*l1**2*l2*l3*(3*l3-2*l1)+gv1*l1*l2*l3**2*(2*l3-3*l1) \
               +gv2*l1**2*l3**2*(l3-l1)
    return facefunG50_3_

def facefunG51_1(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    (gv1, gv2, gv3) = cov_coord_vecs[local_facenodes]
    def facefunG51_1_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return 2*gv1*l1*l2**2*l3**2+2*gv2*l1**2*l2*l3**2+2*gv3*l1**2*l2**2*l3
    return facefunG51_1_


def facefunR20_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    def facefunR20_1_D_(coord):
        try: (l0,l1,l2) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -3*l0*cv23 - 3*l1*cv13
    return facefunR20_1_D_

def facefunR20_2_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv02 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    cv01 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    def facefunR20_2_D_(coord):
        try: (l0,l1,l2) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return 3*l2*cv01 + 3*l1*cv02
    return facefunR20_2_D_

def facefunR30_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR30_1_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return  l1*(l1 - 2*l2)*cv23 + l2*(2*l1 - l2)*cv13
    return facefunR30_1_D_

def facefunR30_2_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR30_2_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return l2*(2*l3 - l2)*cv13 + l3*(l3 - 2*l2)*cv12
    return facefunR30_2_D_

def facefunR30_3_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    def facefunR30_3_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return  l3*(l3-2*l1)*cv12 + l1*(l1-2*l3)*cv23
    return facefunR30_3_D_

def facefunR40_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR40_1_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return cv23*l1*(l2-l1)*(3*l2-l1)+cv13*l2*(l2-3*l1)*(l2-l1)
    return facefunR40_1_D_

def facefunR40_2_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR40_2_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -cv13*l2*(l3-l2)*(3*l3-l2)-cv12*l3*(l3-3*l2)*(l3-l2)
    return facefunR40_2_D_

def facefunR40_3_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    def facefunR40_3_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return cv12*l3*(l3-3*l1)*(l3-l1)-cv23*l1*(l3-l1)*(3*l3-l1)
    return facefunR40_3_D_

def facefunR41_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    def facefunR41_1_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -cv23*l1*(4*l2*l3-l1*l3-l1*l2) + \
               cv12*l3*(l2*l3+l1*l3-4*l1*l2) - cv13*l2*(l2*l3-4*l1*l3+l1*l2)
    return facefunR41_1_D_

def facefunR50_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR50_1_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -cv23*l1*(l2-l1)**2*(4*l2-l1)-cv13*l2*(l2-4*l1)*(l2-l1)**2
    return facefunR50_1_D_

def facefunR50_2_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR50_2_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return cv13*l2*(l3-l2)**2*(4*l3-l2)+cv12*l3*(l3-4*l2)*(l3-l2)**2
    return facefunR50_2_D_

def facefunR50_3_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR50_3_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return cv12*l3*(l3-4*l1)*(l3-l1)**2-cv23*l1*(l3-l1)**2*(4*l3-l1)
    return facefunR50_3_D_

def facefunR51_1_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR51_1_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return -6*cv13*l1*l2**2*l3-6*cv23*l1**2*l2*l3
    return facefunR51_1_D_

def facefunR51_2_D(local_faceno, local_facenodes):
    local_facenodes = local_facenodes[local_faceno]
    cv23 = conv_coord_vecs[tuple(local_facenodes[[1,2]])]
    cv12 = conv_coord_vecs[tuple(local_facenodes[[0,1]])]
    cv13 = conv_coord_vecs[tuple(local_facenodes[[0,2]])]
    def facefunR51_2_D_(coord):
        try: (l1,l2,l3) = coord[local_facenodes]
        except TypeError: raise TypeError('coord must be a numpy array')
        return 6*cv12*l1*l2*l3**2+6*cv13*l1*l2**2*l3
    return facefunR51_2_D_



(gv1, gv2, gv3, gv4) = cov_coord_vecs
def volfunR3_1(coord):
    (l1,l2,l3,l4) = coord
    return  gv4*l1*l2*l3
def volfunR3_2(coord):
    (l1,l2,l3,l4) = coord
    return  gv1*l2*l3*l4
def volfunR3_3(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3*l4

def volfunR4_1(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1**2*l2*l3      

def volfunR4_2(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**2*l3*l4      

def volfunR4_3(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3**2*l4      

def volfunR4_4(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1*l2**2*l3      

def volfunR4_5(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1*l2*l3**2      

def volfunR4_6(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3**2*l4      

def volfunR4_7(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3*l4**2      

def volfunR4_8(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3*l4**2

def volfunR5_1(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1**3*l2*l3

def volfunR5_2(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**3*l3*l4

def volfunR5_3(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3**3*l4

def volfunR5_4(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1**2*l2**2*l3

def volfunR5_5(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1**2*l2*l3**2

def volfunR5_6(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**2*l3**2*l4

def volfunR5_7(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**2*l3*l4**2

def volfunR5_8(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3**2*l4**2

def volfunR5_9(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1*l2**3*l3

def volfunR5_10(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1*l2**2*l3**2

def volfunR5_11(coord):
    (l1,l2,l3,l4) = coord
    return gv4*l1*l2*l3**3

def volfunR5_12(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3**3*l4

def volfunR5_13(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3**2*l4**2

def volfunR5_14(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3*l4**3

def volfunR5_15(coord):
    (l1,l2,l3,l4) = coord
    return gv2*l1*l3*l4**3



(cv12, cv13, cv14, cv23, cv24, cv34) = [conv_coord_vecs[i,j]
                                        for i in range(4) for j in range(4)
                                        if i < j]
def volfunR3_1_D(coord):
    (l1,l2,l3,l4) = coord
    return  cv14*l2*l3+cv24*l1*l3+cv34*l1*l2

def volfunR3_2_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3*l4-cv13*l2*l4-cv14*l2*l3

def volfunR3_3_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3*l4-cv23*l1*l4-cv24*l1*l3

def volfunR4_1_D(coord):
    (l1,l2,l3,l4) = coord
    return 2*cv14*l1*l2*l3+cv24*l1**2*l3+cv34*l1**2*l2

def volfunR4_2_D(coord):
    (l1,l2,l3,l4) = coord
    return -2*cv12*l2*l3*l4-cv13*l2**2*l4-cv14*l2**2*l3

def volfunR4_3_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3**2*l4-2*cv23*l1*l3*l4-cv24*l1*l3**2

def volfunR4_4_D(coord):
    (l1,l2,l3,l4) = coord
    return cv14*l2**2*l3+2*cv24*l1*l2*l3+cv34*l1*l2**2

def volfunR4_5_D(coord):
    (l1,l2,l3,l4) = coord
    return cv14*l2*l3**2+cv24*l1*l3**2+2*cv34*l1*l2*l3

def volfunR4_6_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3**2*l4-2*cv13*l2*l3*l4-cv14*l2*l3**2

def volfunR4_7_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3*l4**2-cv13*l2*l4**2-2*cv14*l2*l3*l4

def volfunR4_8_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3*l4**2-cv23*l1*l4**2-2*cv24*l1*l3*l4

def volfunR5_1_D(coord):
    (l1,l2,l3,l4) = coord
    return 3*cv14*l1**2*l2*l3+cv24*l1**3*l3+cv34*l1**3*l2

def volfunR5_2_D(coord):
    (l1,l2,l3,l4) = coord
    return -3*cv12*l2**2*l3*l4-cv13*l2**3*l4-cv14*l2**3*l3

def volfunR5_3_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3**3*l4-3*cv23*l1*l3**2*l4-cv24*l1*l3**3

def volfunR5_4_D(coord):
    (l1,l2,l3,l4) = coord
    return 2*cv14*l1*l2**2*l3+2*cv24*l1**2*l2*l3+cv34*l1**2*l2**2

def volfunR5_5_D(coord):
    (l1,l2,l3,l4) = coord
    return 2*cv14*l1*l2*l3**2+cv24*l1**2*l3**2+2*cv34*l1**2*l2*l3

def volfunR5_6_D(coord):
    (l1,l2,l3,l4) = coord
    return -2*cv12*l2*l3**2*l4-2*cv13*l2**2*l3*l4-cv14*l2**2*l3**2

def volfunR5_7_D(coord):
    (l1,l2,l3,l4) = coord
    return -2*cv12*l2*l3*l4**2-cv13*l2**2*l4**2-2*cv14*l2**2*l3*l4

def volfunR5_8_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3**2*l4**2-2*cv23*l1*l3*l4**2-2*cv24*l1*l3**2*l4

def volfunR5_9_D(coord):
    (l1,l2,l3,l4) = coord
    return cv14*l2**3*l3+3*cv24*l1*l2**2*l3+cv34*l1*l2**3

def volfunR5_10_D(coord):
    (l1,l2,l3,l4) = coord
    return cv14*l2**2*l3**2+2*cv24*l1*l2*l3**2+2*cv34*l1*l2**2*l3

def volfunR5_11_D(coord):
    (l1,l2,l3,l4) = coord
    return cv14*l2*l3**3+cv24*l1*l3**3+3*cv34*l1*l2*l3**2

def volfunR5_12_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3**3*l4-3*cv13*l2*l3**2*l4-cv14*l2*l3**3

def volfunR5_13_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3**2*l4**2-2*cv13*l2*l3*l4**2-2*cv14*l2*l3**2*l4

def volfunR5_14_D(coord):
    (l1,l2,l3,l4) = coord
    return -cv12*l3*l4**3-cv13*l2*l4**3-3*cv14*l2*l3*l4**2

def volfunR5_15_D(coord):
    (l1,l2,l3,l4) = coord
    return cv12*l3*l4**3-cv23*l1*l4**3-3*cv24*l1*l3*l4**2

def volfunG3(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3*l4 + gv2*l1*l3*l4 + gv3*l1*l2*l4 + gv4*l1*l2*l3

def volfunG4_1(coord):
    (l1,l2,l3,l4) = coord
    return 2*gv1*l1*l2*l3*l4+gv2*l1**2*l3*l4+gv3*l1**2*l2*l4+gv4*l1**2*l2*l3

def volfunG4_2(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**2*l3*l4+2*gv2*l1*l2*l3*l4+gv3*l1*l2**2*l4+gv4*l1*l2**2*l3

def volfunG4_3(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3**2*l4+gv2*l1*l3**2*l4+2*gv3*l1*l2*l3*l4+gv4*l1*l2*l3**2

def volfunG5_1(coord):
    (l1,l2,l3,l4) = coord
    return 3*gv1*l1**2*l2*l3*l4+gv2*l1**3*l3*l4+gv3*l1**3*l2*l4+gv4*l1**3*l2*l3

def volfunG5_2(coord):
    (l1,l2,l3,l4) = coord
    return 2*gv1*l1*l2**2*l3*l4+2*gv2*l1**2*l2*l3*l4+gv3*l1**2*l2**2*l4+gv4*l1**2*l2**2*l3

def volfunG5_3(coord):
    (l1,l2,l3,l4) = coord
    return 2*gv1*l1*l2*l3**2*l4+gv2*l1**2*l3**2*l4+2*gv3*l1**2*l2*l3*l4+gv4*l1**2*l2*l3**2

def volfunG5_4(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**3*l3*l4+3*gv2*l1*l2**2*l3*l4+gv3*l1*l2**3*l4+gv4*l1*l2**3*l3

def volfunG5_5(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2**2*l3**2*l4+2*gv2*l1*l2*l3**2*l4+2*gv3*l1*l2**2*l3*l4+gv4*l1*l2**2*l3**2

def volfunG5_6(coord):
    (l1,l2,l3,l4) = coord
    return gv1*l2*l3**3*l4+gv2*l1*l3**3*l4+3*gv3*l1*l2*l3**2*l4+gv4*l1*l2*l3**3
      

def basis_set(order, mixed=True, btype=None):
    assert(order > 0)
    if order > 5:
        raise NotImplementedError
    edge_toget = order
    if not mixed: edge_toget += 1
    # Concatenate all the required edge functions
    edgefuns, edgefuns_D = [reduce(lambda x,y: x+y, all_funs[0:edge_toget])
                            for all_funs in (all_edgefuns, all_edgefuns_D)]
    basisSet=Struct(fns=dict(edge=edgefuns), fns_D=dict(edge=edgefuns_D))
    if order >= 2:                      # Face functions only for order >= 2
        face_toget_G = face_toget_R = order
        if mixed: face_toget_G -= 1     # Gradient functions excluded for mixed
        facefuns, facefuns_D = ([], [])
        for i in range(2,face_toget_R+1):
            facefuns.extend(all_facefunsR[i])
            facefuns_D.extend(all_facefunsR_D[i])
            if i <= face_toget_G:
                facefuns.extend(all_facefunsG[i])
                facefuns_D.extend(all_facefunsG_D[i])
        basisSet.fns['face'] = facefuns
        basisSet.fns_D['face'] = facefuns_D
    if order >= 3:                      # Volume funcs only for order >= 3
        vol_toget_G = vol_toget_R = order
        if mixed: vol_toget_G -= 1
        volfuns, volfuns_D = ([], [])
        for i in range(3, vol_toget_R+1):
            volfuns.extend(all_volfunsR[i])
            volfuns_D.extend(all_volfunsR_D[i])
            if i <= vol_toget_G:
                volfuns.extend(all_volfunsG[i])
                volfuns_D.extend(all_volfunsG_D[i])
        basisSet.fns['vol'] = volfuns
        basisSet.fns_D['vol'] = volfuns_D
    if mixed: order -= 0.5
    basisSet.info = Struct(form=1, type='webb99', order=order)
    return basisSet
    
edgefuns0 = tuple(partial(edgefun0, i) for i in range(6))
edgefuns0_D = tuple(partial(edgefun0_D, i) for i in range(6))
edgefuns1 = tuple(partial(edgefun1, i) for i in range(6))
edgefuns1_D = (zero_D,)*6
edgefuns2 = tuple(partial(edgefun2, i) for i in range(6))
edgefuns2_D = (zero_D,)*6
edgefuns3 = tuple(partial(edgefun3, i) for i in range(6))
edgefuns3_D = (zero_D,)*6
edgefuns4 = tuple(partial(edgefun4, i) for i in range(6))
edgefuns4_D = (zero_D,)*6
edgefuns5 = tuple(partial(edgefun5, i) for i in range(6))
edgefuns5_D = (zero_D,)*6

all_edgefuns = [edgefuns0, edgefuns1, edgefuns2, edgefuns3, edgefuns4,
                edgefuns5]
all_edgefuns_D = [edgefuns0_D, edgefuns1_D, edgefuns2_D, edgefuns3_D,
                  edgefuns4_D, edgefuns5_D]

ffs_chain = lambda *ffs: tuple(partial(fun, i) for fun in ffs for i in range(4))

facefunsR2 = ffs_chain(facefunR20_1, facefunR20_2)
facefunsR2_D = ffs_chain(facefunR20_1_D, facefunR20_2_D)

facefunsR3 = ffs_chain(facefunR30_1, facefunR30_2, facefunR30_3)
facefunsR3_D = ffs_chain(facefunR30_1_D, facefunR30_2_D, facefunR30_3_D)
facefunsR4 = ffs_chain(facefunR40_1, facefunR40_2, facefunR40_3,
                       facefunR41_1)
facefunsR4_D = ffs_chain(facefunR40_1_D, facefunR40_2_D, facefunR40_3_D,
                         facefunR41_1_D)
facefunsR5 = ffs_chain(facefunR50_1, facefunR50_2, facefunR50_3,
                       facefunR51_1, facefunR51_2)
facefunsR5_D = ffs_chain(facefunR50_1_D, facefunR50_2_D, facefunR50_3_D,
                       facefunR51_1_D, facefunR51_2_D)


facefunsG2 = ffs_chain(facefunG20_1)
facefunsG2_D = (zero_D,)*4
facefunsG3 = ffs_chain(facefunG30_1, facefunG30_2)
facefunsG3_D = (zero_D,)*8

facefunsG4 = ffs_chain(facefunG40_1, facefunG40_2,facefunG40_3)
facefunsG4_D = (zero_D,)*len(facefunsG4)
facefunsG5 = ffs_chain(facefunG50_1, facefunG50_2,facefunG50_3,
                       facefunG51_1)
facefunsG5_D = (zero_D,)*len(facefunsG5)

volfunsR3 = (volfunR3_1,volfunR3_2,volfunR3_3)
volfunsR3_D = (volfunR3_1_D,volfunR3_2_D,volfunR3_3_D)

volfunsR4 = (volfunR4_1, volfunR4_2, volfunR4_3, volfunR4_4, volfunR4_5,
             volfunR4_6, volfunR4_7, volfunR4_8, )
volfunsR4_D = (volfunR4_1_D, volfunR4_2_D, volfunR4_3_D, volfunR4_4_D, volfunR4_5_D,
               volfunR4_6_D, volfunR4_7_D, volfunR4_8_D, )
volfunsR5 = (volfunR5_1, volfunR5_2, volfunR5_3, volfunR5_4, volfunR5_5,
             volfunR5_6, volfunR5_7, volfunR5_8, volfunR5_9, volfunR5_10,
             volfunR5_11, volfunR5_12, volfunR5_13, volfunR5_14,
             volfunR5_15)
volfunsR5_D = (volfunR5_1_D, volfunR5_2_D, volfunR5_3_D, volfunR5_4_D, volfunR5_5_D,
             volfunR5_6_D, volfunR5_7_D, volfunR5_8_D, volfunR5_9_D,
             volfunR5_10_D, volfunR5_11_D, volfunR5_12_D, volfunR5_13_D,
             volfunR5_14_D, volfunR5_15_D)

volfunsG3 = (volfunG3,)
volfunsG3_D = (zero_D(),)
volfunsG4 = (volfunG4_1,volfunG4_2,volfunG4_3,)
volfunsG4_D = (zero_D(),)*3
volfunsG5 = (volfunG5_1, volfunG5_2, volfunG5_3,
             volfunG5_4, volfunG5_5, volfunG5_6)
volfunsG5_D = (zero_D(),)*6

all_facefunsR = {2:facefunsR2,
                 3:facefunsR3,
                 4:facefunsR4,
                 5:facefunsR5,
                 }
all_facefunsR_D = {2:facefunsR2_D,
                   3:facefunsR3_D,
                   4:facefunsR4_D,
                   5:facefunsR5_D,
                   }
all_facefunsG = {2:facefunsG2,
                 3:facefunsG3,
                 4:facefunsG4,
                 5:facefunsG5,
                 }
all_facefunsG_D = {2:facefunsG2_D,
                   3:facefunsG3_D,
                   4:facefunsG4_D,
                   5:facefunsG5_D,
                   }

all_volfunsR = {3:volfunsR3,
                4:volfunsR4,
                5:volfunsR5,
                }
all_volfunsR_D = {3:volfunsR3_D,
                  4:volfunsR4_D,
                  5:volfunsR5_D,
                  }

all_volfunsG = {3:volfunsG3,
                4:volfunsG4,
                5:volfunsG5,
                }
all_volfunsG_D = {3:volfunsG3_D,
                  4:volfunsG4_D,
                  5:volfunsG5_D,
                  }
