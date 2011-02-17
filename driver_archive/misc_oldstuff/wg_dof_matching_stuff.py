from __future__ import division
import numpy as N
from scipy import integrate

from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, DifferentialForm
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities

cb = DifferentialForm.constrained_on_boundary
free = DifferentialForm.allfree

class RoughWGMatcher(object):
    DiscretiserModule = Discretiser
    SubDimDiscretiserModule = SubDimDiscretiser
    SubDimDiscretiserEntitiesModule = SubDimDiscretiserEntities
    SubDimMeshModule = SubDimMesh
    eps = 1e-10
    a = 1.0
    b = 0.25
    l = 10.
    edge1_x = staticmethod(close_to_point(0, eps))
    edge2_x = staticmethod(close_to_point(a, eps))
    edge3_y = staticmethod(close_to_point(b, eps))
    edge4_y = staticmethod(close_to_point(0, eps))
    z_port = staticmethod(close_to_point(0, eps))

    fs = lambda self, face: N.all(self.z_port(face.nodeCoords[:,2]))
    onPort = fs
    subFree = lambda self, edge: not N.all([
        self.edge1_x(x) or self.edge2_x(x) or self.edge3_y(y) or self.edge4_y(y)
        for (x,y,z) in edge.nodeCoords])
    
    superFree = lambda self, edge: cb(edge) 

    def __init__(self, mesh):
        self.mesh = mesh
        subsurf = self.SubDimMeshModule.SubSurface(mesh, self.fs)
        disc = self.DiscretiserModule.setup_PformDiscretiser(mesh, 1, freeFun=self.superFree)
        subGeomEntities = {
            'edge': self.SubDimDiscretiserEntitiesModule.SubDimDiscretiserEntityList(
            self.SubDimDiscretiserEntitiesModule.Edge(mesh=subsurf, freefun=self.subFree,
                                           attrs=subsurf.edges.list_repr())),}
        subDisc = self.SubDimDiscretiserModule.PformSubDimDiscretiser(
            1, subsurf, subGeomEntities, self.DiscretiserModule.Permuter, disc)

        subdisc_edges = subDisc.geomEntities['edge']
        work_edgenos = N.array([e.superNo for e in subdisc_edges
                                if e.freeNo >= 0])
        work_edgeno_set = set(work_edgenos)
        work_edges = disc.geomEntities['edge'][work_edgenos]
        
        dofs = N.zeros(len(work_edgenos), N.float64)

        def E_wg(r):
            return [0,N.sin(N.pi*r[0]/self.a),0]

        for i,e in enumerate(work_edges):
            (l0, l1) = e.nodeCoords
            evec = l1-l0
            f = lambda t: N.dot(E_wg(l0 + (l1-l0)*t), evec)
            dofs[i] = integrate.quad(f, 0, 1)[0]

        self.dofs = dofs
        self.work_edgenos = work_edgenos
        self.work_edgeno_set = work_edgeno_set
        

class RoughCoaxMatcher(object):
    eps = 1e-10
    r_outter = 1.
    r_inner = 2/3.
    E0 = 1.
    on_outter_r = staticmethod(close_to_point(r_outter, eps))
    on_inner_r = staticmethod(close_to_point(r_inner, eps))

    def on_outter(self, pos):
        r = N.linalg.norm(pos[0:2])
        eps = self.eps
        assert (r <= self.r_outter+eps and r >= self.r_inner-eps)
        return self.on_outter_r(r)

    def on_inner(self, pos):
        r = N.linalg.norm(pos[0:2])
        eps = self.eps
        assert (r <= self.r_outter+eps and r >= self.r_inner-eps)
        return self.on_inner_r(r)

    def E_val(self,pos):
        r = N.linalg.norm(pos[0:2])
        r_hat = (pos)/r
        r_hat[2] = 0
        return self.E0/r*r_hat

    z_port = 0.
    on_z_port = staticmethod(close_to_point(z_port, eps))
    fs = lambda self, face: N.all(self.on_z_port(face.nodeCoords[:,2]))
    onPort = fs
    def subFree(self, edge):
        on_o, on_i = [], []
        for pos in edge.nodeCoords:
            on_o.append(self.on_outter(pos))
            on_i.append(self.on_inner(pos))
        return not (N.all(on_o) or N.all(on_i))
        
    superFree = lambda self, edge: cb(edge) 

    def __init__(self, mesh):
        self.mesh = mesh
        subsurf = SubDimMesh.SubSurface(mesh, self.fs)
        disc = Discretiser.setup_PformDiscretiser(mesh, 1, freeFun=self.superFree)
        subGeomEntities = {'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            SubDimDiscretiserEntities.Edge(mesh=subsurf, freefun=self.subFree,
                                           attrs=subsurf.edges.list_repr())),}
        subDisc = SubDimDiscretiser.PformSubDimDiscretiser(
            1, subsurf, subGeomEntities, Discretiser.Permuter, disc)

        subdisc_edges = subDisc.geomEntities['edge']
        work_edgenos = N.array([e.superNo for e in subdisc_edges
                                if e.freeNo >= 0])
        work_edgeno_set = set(work_edgenos)
        work_edges = disc.geomEntities['edge'][work_edgenos]
        
        dofs = N.zeros(len(work_edgenos), N.float64)

        for i,e in enumerate(work_edges):
            (l0, l1) = e.nodeCoords
            evec = l1-l0
            f = lambda t: N.dot(self.E_val(l0 + (l1-l0)*t), evec)
            dofs[i] = integrate.quad(f, 0, 1)[0]

        self.dofs = dofs
        self.work_edgenos = work_edgenos
        self.work_edgeno_set = work_edgeno_set
        self.subDisc = subDisc
        self.subSurf = subsurf
