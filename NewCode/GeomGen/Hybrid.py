from __future__ import division
from itertools import izip
import numpy as N

from NewCode import ProxyList
from NewCode.Meshes import BrickMesh, BrickMeshGen
import NewCode.Mesh as Mesh
from NewCode.Meshes.MeshIO import Gmsh
from NewCode.Utilities import in_box

FakeEdge = Mesh.ProxyEdge

def make_brick_disabled_set_rect(hexmesh, tetrange_n, tetrange_p):
    disabled_bricks = set()
    in_tetregionp = in_box(tetrange_n, tetrange_p)
    for el in hexmesh.elements:
        centroid = N.average(el.nodeCoords, axis=0)
        if  in_tetregionp(centroid):
            disabled_bricks.add(el.index)
    return disabled_bricks

class DeadElements(object):
    has_extra_geom_nodes = False
    has_extra_geom_edges = False
    has_extra_geom_faces = False
    has_extra_geom_circsegs = False
    def __init__(self, hex_mesh):
        self.hex_mesh = hex_mesh
        self.background_mesh = hex_mesh
        self.h = hex_mesh.gridStepSize
        self.offset = hex_mesh.offset
    
    def on_hbdry_hex_geom(self, ent):
        con2 = ent.connect2elem
        in_dis = N.array([con in self.dead_background_elements
                          for con in con2], N.bool8)
        return not N.all(in_dis) and N.any(in_dis) 

    def on_hbdry_hex(self, ent):
        con2 = ent.connect2elem
        in_dis = N.array([con in self.dead_background_elements
                          for con in con2], N.bool8)
        con_mesh_bdry = con2 <= -1
        return  not N.all(in_dis | con_mesh_bdry) and N.any(in_dis)

    def init_geom(self, include_meshbdry_els=True):
        hex_mesh = self.background_mesh 
        on_hbdry = self.on_hbdry_hex_geom if include_meshbdry_els \
                   else self.on_hbdry_hex
        self.hyb_nodes = hex_mesh.nodes[:]
        self.hyb_edgenodes = [e.nodes for e in hex_mesh.edges if on_hbdry(e)]
        self.hyb_edges = ProxyList.ProxyList(
            FakeEdge(dict(nodes=self.hyb_edgenodes,
                          nodeCoords=self.hyb_nodes)))
        hyb_edgemap = dict((tuple(n), i) for i,n in enumerate(self.hyb_edgenodes))
        self.bdry_faces = bdry_faces = N.array(
            [i for i,f in enumerate(hex_mesh.faces) if on_hbdry(f)], N.int32)
        self.hyb_faceedges = N.array([[hyb_edgemap[tuple(hex_mesh.edges[edgeno].nodes)]
                          for edgeno in hf.edgenos]
                         for hf in hex_mesh.faces[bdry_faces]], N.int32)
        self.hyb_faceedges[:] += 1
        self.hyb_faceedges[:,2:] *= -1
        
        self.hyb_volfaces = N.arange(len(self.hyb_faceedges))+1

    def add_extra_geom_nodes(self, extra_nodes, h=None):
        self.has_extra_geom_nodes = True
        self.extra_node_offs=len(self.hyb_nodes)
        self.extra_nodes = extra_nodes
        self.extra_h = h
        
    def add_extra_geom_edges(self, extra_edges):
        self.has_extra_geom_edges = True
        self.extra_edge_offset = len(self.hyb_edgenodes)
        self.extra_edges = ProxyList.ProxyList(
            FakeEdge(dict(nodes=extra_edges+self.extra_node_offs,
                          nodeCoords=self.extra_nodes)))

    def add_extra_geom_circsegs(self, circ_segs):
        self.has_extra_geom_circsegs = True
        self.extra_circ_segs = circ_segs + self.extra_node_offs

    def add_extra_geom_faces(self, extra_faces ):
        self.has_extra_geom_faces = True
        self.extra_faces = extra_faces

    def write_geom(self, fileo, h_geo=None, write_vols=True):
        """Write geometry in gmsh format to file-like object fileo
        """
        h_mesh = self.h
        if h_geo is None: h_geo = N.average(h_mesh)
        fileo.writelines(Gmsh.nodes_to_geo(self.hyb_nodes,h=h_geo))
        fileo.writelines(Gmsh.edges_to_geo(self.hyb_edges))
        fileo.writelines(Gmsh.hexfaceedges_to_geo(self.hyb_faceedges))
        if write_vols: fileo.writelines(Gmsh.faces_to_geo([self.hyb_volfaces]))
        fileo.write("Physical Surface(1) = {%s};" %
                       ','.join(map(str, self.hyb_volfaces)) + '\n' )
        if write_vols: fileo.write("Physical Volume(1) = {1};\n")
        self.write_extra_geom(fileo, h_geo)

    def write_extra_geom(self, fileo, h_geo):
        if self.has_extra_geom_nodes:
            extra_h = self.extra_h if self.extra_h is not None else h_geo
            fileo.writelines(Gmsh.nodes_to_geo(
                self.extra_nodes, offset=self.extra_node_offs, h=extra_h))
        if self.has_extra_geom_edges:
            fileo.writelines(Gmsh.edges_to_geo(self.extra_edges,
                                               offset=self.extra_edge_offset))
            
        if self.has_extra_geom_circsegs:
            fileo.writelines(Gmsh.circsegs_to_geo(self.extra_circ_segs))

    def init_dead_background_element_set(self):
        deadel_fn = self.dead_element_selector
        self.dead_background_elements = set(
            el.index for el in self.hex_mesh.elements
            if deadel_fn(el))
            
class NodeInDeadElements(DeadElements):
    def set_dead_node_selector(self, dead_node_fn):
        self.dead_node_selector = dead_node_fn

    def dead_element_selector(self, el):
        dns = self.dead_node_selector
        return N.any([dns(n) for n in el.nodeCoords])

class ElementCloseDeadElements(NodeInDeadElements):
    def set_dead_node_selector(self, lowres_dead_node_fn, highres_dead_node_fn,
                               no_sub_divs=[4,4,4]):
        """
        lowres_dead_node_fn --- true for nodes up to 1/2 of the element
                                diagonal from the geometry.
        highres_dead_node_fn --- true for nodes within the required minimum
                                 distance
        """
        self.lowres_dead_node_selector = lowres_dead_node_fn
        self.highres_dead_node_selector = highres_dead_node_fn
        self.sub_divs = [N.linspace(0, 1., nosub) for nosub in no_sub_divs]

    def dead_element_selector(self, el):
        l_dns = self.lowres_dead_node_selector
        h_dns = self.highres_dead_node_selector
        el_nds = el.nodeCoords
        mid_pt = N.average(el_nds, axis=0)
        if not l_dns(mid_pt): return False
        if h_dns(mid_pt): return True
        x_divs, y_divs, z_divs = self.sub_divs
        coord_min = el_nds[0]
        coord_max = el_nds[-1]
        for xd in x_divs:
            for yd in y_divs:
                for zd in z_divs:
                    d = N.array([xd,yd,zd], N.float64)
                    if h_dns(coord_min*(1-d) + coord_max*(d)): return True

class SimpleSurroundedRect(DeadElements):

    def set_tet_range(self, tetrange_n, tetrange_p):
        tr_n_unoff = tetrange_n - self.offset
        tr_p_unoff = tetrange_p - self.offset
        h = self.h
        tr_n_unoff_q = N.floor(tr_n_unoff/h)*h
        tr_p_unoff_q = N.ceil(tr_p_unoff/h)*h
        self.tetrange_n = tr_n_unoff_q + self.offset
        self.tetrange_p = tr_p_unoff_q + self.offset
    
    def init_dead_background_element_set(self):
        self.dead_background_elements = make_brick_disabled_set_rect(
            self.hex_mesh, self.tetrange_n, self.tetrange_p)
        


class OffsetRect(DeadElements):
    def __init__(self, tet_geom_size_q, volsize_q, h,
                 tet_offset_q=[0,0,0]):
        self.tet_geom_size_q = tet_geom_size_q
        self.volsize_q = volsize_q
        self.h = h
        self.tet_offset_q=N.array(tet_offset_q, N.int32)
        
    def init_dead_background_element_set(self):
        self.tet_geom_range_n = self.tet_offset_q*self.h
        self.tet_geom_range_p = self.tet_geom_range_n + self.h*self.tet_geom_size_q
        
        self.dead_background_elements = make_brick_disabled_set_rect(
            self.background_mesh, self.tet_geom_range_n, self.tet_geom_range_p)
    
    def init_background_mesh(self, pmlcells=0):
        a,b,c = self.volsize_q*self.h
        self.background_mesh = BrickMesh.Mesh(
            BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, self.h))
        
        print 'background_mesh elements: ', len(self.background_mesh.elements)
        

class SurroundedRect(DeadElements):
    """Construct hybrid geometry of tet-meshed rectangle surrounded by hex elements

    Note formats are suitable for use with Gmsh, i.e. 1 based rather than
    0-based indices
    """

    def __init__(self, tet_geom_minsize, freespace_minsize, h):
        self.tet_geom_minsize = tet_geom_minsize
        self.freespace_minsize = freespace_minsize
        self.h = h

    def qsize(self, size):
        """Quantises size to hex-mesh h such that qsize*h >= size
        """
        return N.int32(N.ceil(size/self.h))


    def init_background_mesh(self, pmlcells=0):
        h = self.h ; qsize = self.qsize
        self.tet_geom_sizeq = qsize(self.tet_geom_minsize)
        self.buffsize_q = qsize(self.freespace_minsize - self.tet_geom_sizeq*h)
        # Must be even, i.e. buffer on both sides of freespace
        self.buffsize_q[:] = [s if s%2 == 0 else s+1 for s in self.buffsize_q]
        assert(N.all(self.buffsize_q > 0))
        self.freespace_sizeq = self.tet_geom_sizeq+self.buffsize_q
        self.tet_geom_range_p = self.tet_geom_sizeq*h - self.tet_geom_sizeq*h/2
        self.tet_geom_range_n = - self.tet_geom_range_p
        self.background_mesh_size = a,b,c = (
            self.freespace_sizeq*h + pmlcells*h*2)
        self.offs = -N.array([a,b,c])/2
        if N.all(pmlcells > 0):
            self.pml_pos = self.freespace_sizeq*h/2

        self.background_mesh = BrickMesh.Mesh(
            BrickMeshGen.make_rect_cavity_brick_listmesh(a,b,c, [h, h, h],
                                                 grid_offset=self.offs))
        
        print 'background_mesh elements: ', len(self.background_mesh.elements)

    def init_dead_background_element_set(self):
        self.dead_background_elements = make_brick_disabled_set_rect(
            self.background_mesh, self.tet_geom_range_n, self.tet_geom_range_p)
        


def femmesh_tris2ent_nodesets(femmesh_tris):
    tri_nodes = femmesh_tris.nodes ; tri_lables = femmesh_tris.material
    nodesets = dict( (int(label), set()) for label in N.unique(tri_lables))
    def tri2edgenodes(trinodes):
        return set(e_nds for e_nds in [(t_nds[0], t_nds[1]),
                                       (t_nds[0], t_nds[2]),
                                       (t_nds[1], t_nds[2])])

    for t_nds, t_lbl in izip(tri_nodes, tri_lables):
        ns = nodesets[int(t_lbl)]
        ns.add(tuple(t_nds))
        ns.update(tri2edgenodes(t_nds))

    return nodesets

class WGSizeCalc(object):
    def __init__(self, h, hybrid_geomlen, geom_min_buflen):
        self.h = h
        self.buflen_q = int(N.max([2, N.ceil(geom_min_buflen/h)]))
        self.hybrid_len_q = int(N.ceil(hybrid_geomlen/h))
        self.TF_len_q = 2*self.buflen_q + self.hybrid_len_q
        self.TF_len = self.TF_len_q*h

class MultiPinMinGeomlens(object):
    def __init__(self, p_cs, r_ps, r_es):
        self.p_cs = p_cs; selfr_ps = r_ps; self.r_es = r_es
        self.p_z_maxind = p_z_maxind = N.argmax(p_cs[:,2])
        self.p_z_minind =p_z_minind = N.argmin(p_cs[:,2])
        self.p_z_max, self.p_z_min = p_cs[[p_z_maxind, p_z_minind], 2]
        self.g_z_max = self.p_z_max + r_ps[p_z_maxind] + r_es[p_z_maxind] 
        self.g_z_min = self.p_z_min - r_ps[p_z_minind] - r_es[p_z_minind] 
        self.min_geomlen = self.g_z_max - self.g_z_min
