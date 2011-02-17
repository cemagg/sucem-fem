from __future__ import division
import numpy as N
import scipy 
from scipy.sparse.linalg.eigen.arpack import speigs
from scipy import integrate

from NewCode.Utilities import Struct, partial, close_to_point
from NewCode import SubDimMesh, BrickSubDimMesh, DifferentialForm
from NewCode.DifferentialForm import Discretiser, SubDimDiscretiser,SubDimDiscretiserEntities
from NewCode.DifferentialForm import BrickSubDimDiscretiser, BrickSubDimDiscretiserEntities
from NewCode.DifferentialForm import BrickDiscretiserEntities, allfree, allconstrained


class WaveguideEigenMatcher(object):
    Permuter = Discretiser.Permuter
    PformSubDimDiscretiser = SubDimDiscretiser.PformSubDimDiscretiser
    SubDimDiscretiserEntityList = SubDimDiscretiserEntities.SubDimDiscretiserEntityList
    SubDimEdge = SubDimDiscretiserEntities.Edge
    SubDimFace = SubDimDiscretiserEntities.Face
    SubSurface= SubDimMesh.SubSurface
    eps = 1e-10
    a = 1.0
    b = 0.25
    edge1_x = staticmethod(close_to_point(0, eps))
    edge2_x = staticmethod(close_to_point(a, eps))
    edge3_y = staticmethod(close_to_point(b, eps))
    edge4_y = staticmethod(close_to_point(0, eps))
    z_port = staticmethod(close_to_point(0, eps))
    sigma = 1.

    fs = lambda self, face: N.all(self.z_port(face.nodeCoords[:,2]))
    onPort = fs
    subFree = lambda self, edge: not N.all([
        self.edge1_x(x) or self.edge2_x(x) or self.edge3_y(y) or self.edge4_y(y)
        for (x,y,z) in edge.nodeCoords])
    
    def __init__(self, disc):
        subsurf = self.SubSurface(disc.mesh, self.fs)
        subGeomEntities = {
            'edge': self.SubDimDiscretiserEntityList(
            self.SubDimEdge(mesh=subsurf, freefun=self.subFree,
                            attrs=subsurf.edges.list_repr()))}
        if 'face' in disc.basisSet.fns:
            subGeomEntities['face'] = self.SubDimDiscretiserEntityList(
                self.SubDimFace(mesh=subsurf, freefun=self.subFree,
                                attrs=subsurf.elements.list_repr()))
        
        self.subDisc = subDisc = self.PformSubDimDiscretiser(
            1, subsurf, subGeomEntities, self.Permuter, disc)

    def solve(self, no_modes=1):
        M = self.subDisc.matrix.mass()
        S = self.subDisc.matrix.stiffness()
        sigma_solve = scipy.sparse.linalg.dsolve.factorized(S - self.sigma*M)
        w,v = speigs.ARPACK_gen_eigs(M.matvec, sigma_solve, M.shape[0],
                                     self.sigma, no_modes)
        sort_ind = N.argsort(w)
        return w[sort_ind], v.T[sort_ind]
        
    def get_dofs(self):
        return self.solve()[1][0]


def gen_E_TE01(a, b):
    def E_TE01(r):
        return [0,N.sin(N.pi*r[0]/a),0]
    return E_TE01


class SurfaceFieldMatcher(object):
    Permuter = Discretiser.Permuter
    PformSubDimDiscretiser = SubDimDiscretiser.PformSubDimDiscretiser
    SubDimDiscretiserEntityList = SubDimDiscretiserEntities.SubDimDiscretiserEntityList
    SubDimEdge = SubDimDiscretiserEntities.Edge
    SubDimFace = SubDimDiscretiserEntities.Face
    SubSurface= SubDimMesh.SubSurface

    def initSubdim(self, disc, fs, sub_freefun, useLU=True, dtype=N.float64):
        self.fs = fs
        self.subSurf = subsurf = self.SubSurface(disc.mesh, self.fs)
        self.subFree = sub_freefun
        self.subGeomEntities = {
            'edge': self.SubDimDiscretiserEntityList(
            self.SubDimEdge(mesh=subsurf, freefun=self.subFree,
                            attrs=subsurf.edges.list_repr()))}
        if 'face' in disc.basisSet.fns:
            self.subGeomEntities['face'] = self.SubDimDiscretiserEntityList(
                self.SubDimFace(mesh=subsurf, freefun=self.subFree,
                                attrs=subsurf.elements.list_repr()))
        self.subDisc = subDisc = self.PformSubDimDiscretiser(
            1, subsurf, self.subGeomEntities, self.Permuter, disc)
        self.DOFs = self.subDisc.newDOFs(dtype=dtype)
        self.DOFs.useLU = useLU

    def matchKnown(self, matchfun):
        self.DOFs.matchFunction(matchfun)
        return self.DOFs.dofArray.copy()


class BrickSurfaceFieldMatcher(SurfaceFieldMatcher):
    PformSubDimDiscretiser = BrickSubDimDiscretiser.PformSubDimDiscretiser
    SubDimDiscretiserEntityList = BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList
    SubDimEdge = BrickSubDimDiscretiserEntities.Edge
    SubDimFace = BrickSubDimDiscretiserEntities.Face
    SubSurface = BrickSubDimMesh.SubSurface

class TangentialPlaneWaveMatcher(object):
    def __init__(self, z_ref, c, drv_fun):
        """z+ traveling x+ polarised plainwave, with zero phase/time ref at z=z_ref"""
        self.z_ref, self.c, self.drv_fun = z_ref, c, drv_fun
        self.E_dir = N.array([1.,0,0], N.float64)

    def set_matcher(self, matcher_fun):
        """Set matcher function, such that direch_dofarr = matcher_fun(matchfun)

        matchfun is supplied by self and represents the plane wave E-field in space
        """
        self.matcher = matcher_fun

    def E0(self, dt, n, r):
        z = r[2]
        return self.E_dir*self.drv_fun(dt*n-(z-self.z_ref)/self.c, 1)

    def get_direchdofs(self, dt, n):
        return self.matcher(partial(self.E0, dt, n))
        
class TangentialFuncProjSurfaceIntegral(SurfaceFieldMatcher):
    # Will only work properly for 1-forms in current form, since the tangential
    # subDimDiscretiser is used to do the surface integration, and also only
    # tet meshes. Also TESTME!!
    epsrel = 1e-11
    espabs = 1e-11
    def __init__(self, disc, fs, sub_freefun):
        self.initSubdim(disc, fs, sub_freefun)
        self.disc = disc
        p_sub = self.subDisc.permuter
        p_super = disc.permuter
        self.superDOFMap = N.zeros(len(self.DOFs.dofArray), N.int32)
        for ent_name, ents in self.subDisc.geomEntities.iteritems():
            ent_glob_super_dofnos = p_super.globalEntityPermutationTable(ent_name)
            ent_glob_sub_dofnos = p_sub.globalEntityPermutationTable(ent_name)
            for i, ent in enumerate(ents):
                subsuper = N.array([
                    (sub_dofno, super_dofno) for sub_dofno, super_dofno
                    in zip(ent_glob_sub_dofnos[i],
                           ent_glob_super_dofnos[ent.superNo])
                    if sub_dofno >= 0], N.int32)
                sub_dofnos, super_dofnos = [subsuper[:,0], subsuper[:,1]] \
                                           if len(subsuper) > 0 else ([],[])
                self.superDOFMap[sub_dofnos] = super_dofnos

    def calcIntegrationVector(self, func):
        
        a,b = 0,1.                      # 1/2 Unit area triangle 
        gfun = lambda x: 0.             # integration domain
        hfun = lambda x: b-x
        for el in self.subDisc.elements:
            size = el.size*2            # Factor 2 for 1/2 unit integration domain
            local_intg = N.array([integrate.dblquad(lambda y,x: size*N.dot(
                el.bfPhysValAtPoint(bf,N.array([x,y,1-x-y], N.float64)),
                func(el.local2global([x,y,1-x-y]))),
                                                    a,b,gfun,hfun,
                                                    epsrel=self.epsrel,
                                                    epsabs=self.epsabs)[0]
                                  for bf in el.basisSet], N.float64)
            p_l, p_g = el.permutation()
            self.DOFs.dofArray[p_g] += local_intg[p_l]

        superDofArray = self.disc.newDOFs().dofArray
        superDofArray[self.superDOFMap] = self.DOFs.dofArray
        return superDofArray

class BrickTangentialFuncProjSurfaceIntegral(BrickSurfaceFieldMatcher,
                                             TangentialFuncProjSurfaceIntegral):
    def calcIntegrationVector(self, func):
        
        a,b = 0,1.                      # Unit area square
        gfun = lambda x: 0.             # integration domain
        hfun = lambda x: 1.
        for el in self.subDisc.elements:
            size = el.size            
            local_intg = N.array([integrate.dblquad(lambda y,x: size*N.dot(
                el.bfPhysValAtPoint(bf,N.array([x,y,1-x,1-y], N.float64)),
                func(el.local2global([x,y,1-x,1-y]))),
                                                    a,b,gfun,hfun,
                                                    epsrel=self.epsrel,
                                                    epsabs=self.epsabs)[0]
                                  for bf in el.basisSet], N.float64)
            p_l, p_g = el.permutation()
            self.DOFs.dofArray[p_g] += local_intg[p_l]

        superDofArray = self.disc.newDOFs().dofArray
        superDofArray[self.superDOFMap] = self.DOFs.dofArray
        return superDofArray


def find_WG_inc_B_facesvols(mesh, E_edges, E_faces, scatteredside_p):
    """Find B face and volume numbers to log for TF/SF incident WG feed

    E_edges and E_faces are the edge and face numbers of mesh that make up the
    E-field incident condition surface.

    scatteredside_p(pt) must return true if pt is on the scattered-field side
    of the log surface.

    Returns Struct(faces= face numbers of mesh for B logging,
                   vols= volume numbers of mesh for B logging.)
    """

    possible_els = N.unique(mesh.edges[E_edges].connect2elem)
    possible_els = possible_els[possible_els >= 0]
    meshels = mesh.elements
    meshfaces = mesh.faces
    B_els = N.array([i for i in possible_els if scatteredside_p(
        N.average(meshels[i].nodeCoords, axis=0))], N.int32)
    #
    E_edgeset = set(E_edges)
    E_faceset = set(E_faces)
    possible_B_faces = N.unique(mesh.elements[B_els].facenos)
    possible_B_faces = possible_B_faces[possible_B_faces >=0]
    B_faces = N.array([i for i in possible_B_faces
                       if i not in E_faceset \
                       and N.any([ei in E_edgeset for ei in meshfaces[i].edgenos])],
                      N.int32)
    return Struct(faces=B_faces, vols=B_els)
                        
                       
class TFSFIncidentField(object):
    def __init__(self, discs, bdry_facesel, in_scattered_p):
        self.discs = discs
        self.mesh = discs.E.mesh
        assert (self.mesh is self.discs.B.mesh)
        self.bdry_facesel = bdry_facesel
        self.in_scattered_p = in_scattered_p

    def init_incsurf(self, incsurf_free):
        self.E_surf_integrator = tif = BrickTangentialFuncProjSurfaceIntegral(
            self.discs.E, self.bdry_facesel, incsurf_free)

    def get_inc_info(self):
        tif = self.E_surf_integrator
        E_inc_dofnos = tif.superDOFMap
        E_inc_facenos = tif.subDisc.elements[:].superFacenos
        E_inc_edgenos = tif.subDisc.mesh.edges[:].superNo
        B_inc_ents = find_WG_inc_B_facesvols(
            self.mesh, E_inc_edgenos, E_inc_facenos, self.in_scattered_p)
        B_disc = self.discs.B
        B_entdofs = []
        for ent in sorted(B_disc.geomEntities.keys()):
            tmp_entdofs = N.unique(
                B_disc.permuter.partialGlobalEntityPermutationTable(
                ent, B_inc_ents[ent+'s']))
            tmp_entdofs = tmp_entdofs[tmp_entdofs >= 0]
            B_entdofs.append(tmp_entdofs)
        B_inc_dofnos = N.hstack(B_entdofs)
        return Struct(E=Struct(dofnos=E_inc_dofnos, ents=Struct(
            faces=E_inc_facenos, edges=E_inc_edgenos)),
                      B=Struct(dofnos=B_inc_dofnos, ents=B_inc_ents))
        
# class TFSFDriveDiscs(object):
#     def __init__(self, inc_field, discs):
#         self.glob_discs = discs
#         def get_disc(quant, disc):
#             assert(disc.diagonalMetric == True)
#             assert(disc.basisSet.info.type == 'cohen98')
#             return self._setup_drv_disc(disc, inc_field[quant].ents)

#             self.inc_discs = Struct((quant, get_disc(quant, disc))
#                                     for quant, disc in discs.items())
            

#     @staticmethod
#     def _setup_drv_disc(disc, entnos):
#         basisSet = disc.basisSet
#         glob_geomEntities = disc.geomEntities
#         order = int(N.ceil(basisSet.info.order))
#         assert (basisSet.info.order < order) # I.e. mixed
#         drv_entsets = Struct((ent, set(entnos[ent+'s']))
#                              for ent in ('edge', 'face', 'vol')
#                              if ent+'s' in entnos)
#         freefun = lambda ent: ent.index in drv_entsets[ent.entity_type] and \
#                   glob_geomEntities[ent.entity_type][ent.index].isFree
#         drv_disc = BrickDiscretiser.setup_PformDiscretiser(
#             disc.mesh, disc.p, order, mixed=True, freeFun=freefun,
#             btype='cohen98', dead_elements=disc.elements.entity.dead_elements)

#         return drv_disc
    
        
