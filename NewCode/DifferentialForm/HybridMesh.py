from __future__ import division

from itertools import izip, chain

from NewCode.DifferentialForm import Discretiser, BrickDiscretiser, allfree, allconstrained
from NewCode.Integration import TriIntegrator, TetProdIntegrator
from NewCode import SubDimMesh, BrickSubDimMesh
from NewCode.DifferentialForm import SubDimDiscretiserEntities, SubDimDiscretiser
from NewCode.DifferentialForm import BrickSubDimDiscretiser, BrickSubDimDiscretiserEntities
from NewCode.SystemMatrix import local_self_projection_matrix, insert_global, eps
from NewCode.Utilities import Struct

import numpy as N
from scipy import sparse, linalg

def gen_hc_dofnos(edges, d_edges):
    no_edges = len(edges) 
    hc_edge_dofnos = N.zeros((no_edges,1), dtype=N.int32)
    dofno_n = 0
    for i, edge in enumerate(edges):
        if i in d_edges: hc_edge_dofnos[i,0] = -1
        else:
            hc_edge_dofnos[i,0] = dofno_n
            dofno_n += 1
    return hc_edge_dofnos

def tet2hc_tet(edges, d_edges, d_edgemap, edge_dofnos):
    """Transform fully linear tet DOFs into CT/LN hex-compatible tets

    Return Value
    ============

    (T, hc_edge_dofnos) where

    T -- Transform matrix, T*orig_dofs -> hex_compat_dofs

    hc_edge_dofnos

      n x m array of edge dofnos for the hex-compat dofs, where n is the number
      of edges in the mesh and m the number of DOFs per edge. Diagonal edges
      are all constrained here.

    """

    edgemap = edges.nodemap
    hc_edge_dofnos = gen_hc_dofnos(edges, d_edges)
    no_newdofs = N.max(hc_edge_dofnos)+1
    no_olddofs = N.max(edge_dofnos)+1

    T = sparse.lil_matrix(shape=(no_newdofs, no_olddofs), dtype=N.float64)

    for i, edge in enumerate(edges):
        if i in d_edges:
            (pen_1, pen_2), (nen_1, nen_2) = d_edgemap[tuple(edge.nodes)]
            penos = [edgemap[tuple(pen_1)], edgemap[tuple(pen_2)]]
            nenos = [edgemap[tuple(nen_1)], edgemap[tuple(nen_2)]]
            olddofno_d0, olddofno_d1 = edge_dofnos[i]
            pd1 = hc_edge_dofnos[penos[0],0]
            pd2 = hc_edge_dofnos[penos[1],0]
            nd1 = hc_edge_dofnos[nenos[0],0]
            nd2 = hc_edge_dofnos[nenos[1],0]
            T[[pd1, pd2, nd1, nd2], olddofno_d0] = 1/2
            T[[pd1, pd2], olddofno_d1] = 1/2
            T[[nd1, nd2], olddofno_d1] = -1/2
            continue
        T[hc_edge_dofnos[i,0],edge_dofnos[i,0]] = 1
    return (T, hc_edge_dofnos)

def make_hybmesh_discs(hyb_mesh, order, freefun, dead_elements=None,
                       constrained_mtet_vols=True):
    on_hbdry = hyb_mesh.on_hbdry
    if dead_elements == None: dead_elements = set()
    otet_freefun = lambda ent: not on_hbdry(ent) and freefun(ent)
    mtet_order = order*2-1              # mtet for Matching Tet

    mtet_volfree = allconstrained if constrained_mtet_vols else allfree
    mtet_disc = Discretiser.setup_PformDiscretiser( 
        hyb_mesh.tet, form=1, order=mtet_order, mixed=False,
        freeFun=freefun,vol_freeFun=mtet_volfree)
    try:
        mtet_disc.setIntegrationRule(mtet_order*2)
        mtet_disc.D().setIntegrationRule(mtet_order*2)
    except KeyError:
        mtet_disc.set_integrator(TetProdIntegrator)
        mtet_disc.D().set_integrator(TetProdIntegrator)
        mtet_disc.setIntegrationRule(mtet_order*2)
        mtet_disc.D().setIntegrationRule(mtet_order*2)
        
    mtet_disc.setFaceIntegrationRule(TriIntegrator(mtet_order*2))
    otet_disc = Discretiser.setup_PformDiscretiser( # otet for Ordinary Tet
        hyb_mesh.tet, form=1, order=order, mixed=True, freeFun=otet_freefun)
    brick_volfree = lambda el: el.index not in dead_elements
    def brick_free(ent):
        con2 = ent.connect2elem
        in_dis = N.array([con in dead_elements for con in con2], N.bool8)
        con_mesh_bdry = con2 <= -1
        return freefun(ent) and not N.all(in_dis | con_mesh_bdry)
    
    brick_disc = BrickDiscretiser.setup_PformDiscretiser(
        hyb_mesh.brick, form=1, order=order, mixed=True, freeFun=brick_free,
        vol_freeFun = brick_volfree, btype='cohen98', dead_elements=dead_elements)
    brick_disc.diagonalise()
    return Struct(mtet=mtet_disc, otet=otet_disc, brick=brick_disc)

def make_hyb_submeshes(hyb_mesh):
    mtet_submesh = SubDimMesh.SubSurface(
        hyb_mesh.tet, faceSelector=hyb_mesh.on_hbdry)
    brick_submesh = BrickSubDimMesh.SubSurface(
        hyb_mesh.brick, faceSelector=hyb_mesh.on_hbdry)
    return Struct(mtet=mtet_submesh, brick=brick_submesh)

def make_hybsubmesh_discs(hyb_submesh, hyb_disc):
    mtet_subGeomEntities = {
        'edge': SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
        SubDimDiscretiserEntities.Edge(mesh=hyb_submesh.mtet, freefun=allfree,
                                       attrs=hyb_submesh.mtet.edges.list_repr()))}
    if 'face' in hyb_disc.mtet.basisSet.fns:
        mtet_subGeomEntities['face'] = SubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            SubDimDiscretiserEntities.Face(mesh=hyb_submesh.mtet, freefun=allfree,
                                           attrs=hyb_submesh.mtet.elements.list_repr()))
    
    mtet_subdisc = SubDimDiscretiser.PformSubDimDiscretiser(
        1, hyb_submesh.mtet, mtet_subGeomEntities, Discretiser.Permuter, hyb_disc.mtet)
    
    brick_subGeomEntities = {
        'edge': BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
        BrickSubDimDiscretiserEntities.Edge(mesh=hyb_submesh.brick, freefun=allfree,
                                            attrs=hyb_submesh.brick.edges.list_repr()))}
    if 'face' in hyb_disc.brick.basisSet.fns:
        brick_subGeomEntities['face'] = BrickSubDimDiscretiserEntities.SubDimDiscretiserEntityList(
            BrickSubDimDiscretiserEntities.Face(mesh=hyb_submesh.brick, freefun=allfree,
                                                attrs=hyb_submesh.brick.elements.list_repr()))
    
    brick_subdisc = BrickSubDimDiscretiser.PformSubDimDiscretiser(
        1, hyb_submesh.brick, brick_subGeomEntities, Discretiser.Permuter, hyb_disc.brick)
    
    return Struct(mtet=mtet_subdisc, brick=brick_subdisc)


def make_hface_permutation(tri_els):
    glob_perms = [el.permutation() for el in tri_els]
    unique_glob_dofs = N.unique(N.hstack([ep[1] for ep in glob_perms]))
    unique_glob_dofs = unique_glob_dofs[unique_glob_dofs >= 0]
    glob2hface = dict((g,i) for i,g in enumerate(unique_glob_dofs))
    tri_eldofs = [N.array([glob2hface[g] for g in elperm[1]], N.int32)
                     for elperm in glob_perms]
    glob_perm = [N.arange(len(unique_glob_dofs)), unique_glob_dofs]
    return [(gp[0], tfp) for gp, tfp
            in zip(glob_perms, tri_eldofs)], glob_perm
    

def rem_constr_perm(perm):
    pl, pg = perm
    pl_n = pl[pg >= 0]
    pg_n = pg[pg >= 0]
    return [pl_n, pg_n]

def get_d_nodes(hface_nodes):
    return hface_nodes[[0,2]]

def get_tface_nodes(hface_nodes):
    return (hface_nodes[[0,3,2]], hface_nodes[[0,1,2]])

class HybridMeshOneformDiscretiser(object):
    def __init__(self, hyb_mesh):
        self.hyb_mesh = hyb_mesh

    def init_discs(self,order,freefun, dead_elements=None, constrained_mtet_vols=True):
        self.discs = make_hybmesh_discs(self.hyb_mesh, order, freefun,
                                        dead_elements=dead_elements,
                                        constrained_mtet_vols=constrained_mtet_vols)

    def init_subdiscs(self):
        self.submeshes = make_hyb_submeshes(self.hyb_mesh)
        self.subdiscs = make_hybsubmesh_discs(self.submeshes, self.discs)

    def init_entity_dofnos(self):
        self.entity_dofnos = Struct()
        mt_perm = self.discs.mtet.permuter ; b_perm = self.discs.brick.permuter
        self.entity_dofnos.mtet = Struct((k,mt_perm.globalEntityPermutationTable(k))
                                         for k in self.discs.mtet.geomEntities)
        self.entity_dofnos.brick = Struct((k,b_perm.globalEntityPermutationTable(k))
                                          for k in self.discs.brick.geomEntities)

    def init_nodofs(self):
        self.nodofs = Struct()
        self.nodofs.mtet = self.discs.mtet.totalDOFs
        self.nodofs.brick = self.discs.brick.totalDOFs
        self.nodofs.hyb = self.subdiscs.brick.totalDOFs
        self.nodofs.otet = self.discs.otet.totalDOFs 
        self.nodofs.tot = self.nodofs.otet + self.nodofs.brick 

class HybridSubdim2GlobDofmaps(object):
    subdim_ent_attr = dict(edge='edges', face='elements')


    def _make_dofmap(self, subdisc, glob_ent_dofnos):
        subdim_ent_attr = self.subdim_ent_attr
        subdisc_glob_dofnos = Struct((k, glob_ent_dofnos[k][
            getattr(subdisc.mesh, subdim_ent_attr[k])[:].superNo])
                                     for k in subdisc.geomEntities)
        subdisc_dofnos = Struct((k,subdisc.permuter.globalEntityPermutationTable(k))
                                for k in subdisc.geomEntities)
        return (subdisc_glob_dofnos, subdisc_dofnos)
    
    def init_dofmaps(self, hyb_subdiscs, glob_entity_dofnos):
        self.dofmap = Struct()
        for k, subdisc in hyb_subdiscs.items():
            subdisc_glob_dofnos, subdisc_dofnos = self._make_dofmap(
                subdisc, glob_entity_dofnos[k])
            self.dofmap[k] = dict((sd,gd) for sd,gd in izip(
                chain(*(nos.flat for nos in subdisc_dofnos.values())),
                chain(*(nos.flat for nos in subdisc_glob_dofnos.values()))))



def make_hface_locmats(t_el, b_el):
    l_M = local_self_projection_matrix(t_el)
    eval_r = t_el.physEvalPoints()
    eval_b_el_l = [b_el.global2local(r) for r in eval_r]
    b_el_physvals = b_el.physValsAtPoints(eval_b_el_l)
    intg = t_el.rule.integrateFun
    l_P = N.array([[intg(N.sum(fn_i*fn_j, axis=1))
                    for fn_j in b_el_physvals]
                   for fn_i in t_el.physVals()], N.float64)
    l_P *= t_el.size
    l_P[N.abs(l_P) < eps] = 0
    return l_M, l_P

def make_hface_locmat(tri_els, b_el):
    saved_indices = tri_els.index        # Hack around ProxyList limitation
    perm_tri_els, hface2glob_tet_perm = make_hface_permutation(
        tri_els)                         # Hack around ProxyList limitation
    tri_els.index = saved_indices
    no_2tf_dofs = len(hface2glob_tet_perm[0])
    lh_M = N.zeros(shape=(no_2tf_dofs, no_2tf_dofs),
                   dtype=N.float64)
    lh_P = N.zeros(shape=(no_2tf_dofs, b_el.noDOFs.element),
                   dtype=N.float64)
    for t_el_i, t_el in enumerate(tri_els):
        l_M, l_P = make_hface_locmats(t_el, b_el)
        tl_perm = t_el.permutation()[0]
        insert_global(lh_M, l_M, perm_tri_els[t_el_i])
        insert_global(lh_P, l_P, perm_tri_els[t_el_i],
                      [N.arange(l_P.shape[1]),
                       N.arange(l_P.shape[1])])

    lh_M_inv = linalg.inv(lh_M)
    lh_hc_T = (N.dot(lh_M_inv, lh_P)).T
    lh_hc_T[N.abs(lh_hc_T) < eps] = 0

    return lh_hc_T, hface2glob_tet_perm

normal_dirs = {(1.,0.,0.):0, (0.,1.,0.):1, (0.,0.,1.):2}
def calc_face_normaldir(face_coords):
    bfc = face_coords
    cp = N.abs(N.cross(bfc[1]-bfc[0], bfc[2]-bfc[0]))
    return normal_dirs[tuple(N.round(cp/N.linalg.norm(cp)))]

class HexTriFaceMatcher(object):
    def __init__(self, submeshes):
        self.trism, self.hexsm = submeshes.mtet, submeshes.brick
        # Search radius should be half of the brick diagonal. Making it just
        # over half for incase :)
        self.search_radius = N.linalg.norm(
            self.hexsm.superMesh.gridStepSize)/1.999
        self.init_kdtrees()
        self.init_hextri_map()
        
    def init_kdtrees(self):
        from Bio.KDTree import KDTree
        hface_centroids = [[],[],[]]
        self.hface_indices = hface_indices = [[],[],[]]
        for hface in self.hexsm.elements:
            hface_coords = hface.nodeCoords
            midpt = N.average(hface_coords, axis=0).astype(N.float32)
            ndir = calc_face_normaldir(hface_coords)
            hface_centroids[ndir].append(midpt)
            hface_indices[ndir].append(hface.index)
        self.kdtrees = kdtrees = []
        for centroids in hface_centroids :
            kdtrees.append(KDTree(dim=3))
            if len(centroids) == 0: continue
            kdtrees[-1].set_coords(N.array(centroids, N.float32))

    def init_hextri_map(self):
        r = self.search_radius ; kdtrees = self.kdtrees
        self.hex_tris = [[] for hface in self.hexsm.elements]
        for tri in self.trism.elements:
            tface_coords = tri.nodeCoords
            mid_pt = N.average(tface_coords, axis=0).astype(N.float32)
            ndir = calc_face_normaldir(tface_coords)
            kdt = kdtrees[ndir]
            kdt.search(mid_pt, r)
            inds = kdt.get_indices() ; radii = kdt.get_radii()
            closest_ind = inds[N.argmin(radii)]
            hface_ind = self.hface_indices[ndir][closest_ind]
            self.hex_tris[hface_ind].append(tri.index)
