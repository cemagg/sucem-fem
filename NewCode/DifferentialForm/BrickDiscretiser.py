from __future__ import division

from itertools import izip
import numpy as N
from NewCode.DifferentialForm import Discretiser, BrickDiscretiserEntities
from NewCode.DifferentialForm import BrickDiscretiserElement
from NewCode.DifferentialForm.Permuter import Permuter
from NewCode.Integration import HexaLobattoIntegrator
from NewCode.Exceptions import AllZero
from Discretiser import allfree

space_dims = 3
def calc_interpolation_pointscoords(disc):
    int_pts = N.zeros((disc.totalDOFs, space_dims),
                       dtype=disc.mesh.nodes.dtype)
    int_coords = N.zeros(disc.totalDOFs, dtype=N.int8)
    coordmap = {(1.,0.,0.):0,
                (0.,1.,0.):1,
                (0.,0.,1.):2}

    el0 = disc.elements[0]
    l_int_pts = [bf.interp_pt for bf in el0.basisSet]
    bf_cdir_perel = [N.array([coordmap[tuple(bf.direction)]
                              for bf in el0.basisSet], N.int8)
                     ]*len(disc.elements)
    g_int_pts_perel = N.array([disc.elements[:].local2global(l)
                               for l in l_int_pts]).swapaxes(0,1)
    for el, pts, cdir in izip(disc.elements, g_int_pts_perel, bf_cdir_perel):
        try: l,g = el.permutation()
        except AllZero: continue
        int_pts[g] = pts[l]
        int_coords[g] = cdir[l] # x/y/z directed BF
    return (int_pts, int_coords)

class BaseBrickDiscretiser(object):
    from NewCode.Integration import HexaIntegrator as Integrator
    defaultIntegrationOrder = 5         # Actually no of intg pts?
    identicalElements = True

    def __init__(self, *names, **kwargs):
        if 'dead_elements' in kwargs:
            dead_elements = kwargs['dead_elements']
            del kwargs['dead_elements']
        else: dead_elements=set()
        super(BaseBrickDiscretiser, self).__init__(*names, **kwargs)
        if dead_elements is None: dead_elements = set()
        self.elements.entity.dead_elements = dead_elements
        
    def interpolationPointsCoords(self):
        try: return self._interpolationPointsCoords
        except AttributeError: pass
        self._interpolationPointsCoords = calc_interpolation_pointscoords(self)
        return self._interpolationPointsCoords
        
    def diagonalise(self):
        binfo = self.basisSet.info
        assert binfo.type == 'cohen98'
        order = int(N.ceil(binfo.order))
        self.set_integrator(HexaLobattoIntegrator)
        self.setIntegrationRule(order+1)
        self.D().set_integrator(HexaLobattoIntegrator)
        self.D().setIntegrationRule(order+1)
        self.diagonalMetric = True
        #self.matrix.sparseType = sparse.dia_matrix
        
    def set_alt_integrator(self, alt_integrator):
        self.alt_integrator = alt_integrator

    def set_alt_IntegrationRule(self, integrationOrder):
        self.elements.entity.set_alt_IntegrationRule(
            self.alt_integrator(integrationOrder))

class PformDiscretiser_D(BaseBrickDiscretiser,
                         Discretiser.PformDiscretiser_D):
    discElClasses = {1:BrickDiscretiserElement.OneformElement_D,
                     2:BrickDiscretiserElement.TwoformElement_D}
    
class PformDiscretiser(BaseBrickDiscretiser,
                       Discretiser.PformDiscretiser):
    Discretiser_D = PformDiscretiser_D
    discElClasses = {1:BrickDiscretiserElement.OneformElement,
                     2:BrickDiscretiserElement.TwoformElement}


def setup_PformDiscretiser(
    mesh, form, order=1, mixed=True,
    freeFun=allfree, vol_freeFun=None, el_in_disc=None, btype='rieben04',
    dead_elements=None):
    """
    Set up a discretiser, using defaults for basis function and Permuter selection

    Frees the user from having to setup the basisSet and DiscretiserEntities
    """
    from NewCode.DifferentialForm.BasisFunction import brickBasisForm
    basisSet = brickBasisForm[form].basis_set(order=order, mixed=mixed, btype=btype)
    geomEntities = BrickDiscretiserEntities.make_geomEntities(
        mesh, basisSet, freeFun, vol_freeFun, el_in_disc)
    return  PformDiscretiser(
        form, mesh, geomEntities, Permuter, basisSet, dead_elements=dead_elements)
