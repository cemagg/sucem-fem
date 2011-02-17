from __future__ import division
from itertools import izip
import numpy as N

def rect2polar(x2, x1):
    r = N.sqrt(x1**2 + x2**2)
    theta = N.arctan2(x2, x1)
    return r, theta

class CylindricalSection(object):
    # Angle covered by metallic part in negative direction from phi_0
    seg_phi = N.pi
    spaceDim = 3
    def __init__(self, p_c, r_p, phi_0, x1=2, x2=0):
        self.p_c = N.array(p_c, N.float64) ; self.r_p = r_p ; self.phi_0 = phi_0
        self.x1 = x1 ; self.x2 = x2
        x3_coords = range(self.spaceDim)
        x3_coords.remove(x1) ; x3_coords.remove(x2)
        self.x3 = x3_coords[0]
        
    def to_polar(self, xyz):
        xyz = xyz - self.p_c
        r, phi = rect2polar(xyz[self.x2], xyz[self.x1])
        return r, phi

    def to_xyz(self, r, phi, x3):
        xyz = self.p_c.copy()
        xyz[self.x1] += N.cos(phi)*r
        xyz[self.x2] += N.sin(phi)*r
        xyz[self.x3] += x3
        return xyz
    
class PointInCylindricalSection(object):
    def __init__(self, pin, r_e=0):
        self.pin = pin ; self.r_e = r_e 
        self.r_o_max = pin.r_p + r_e
        self.phi_e = N.arcsin(r_e/self.r_o_max)
        
    def r_o(self, phi):
        phi_0, phi_e, r_o = self.pin.phi_0, self.phi_e, self.r_o_max
        seg_phi = self.pin.seg_phi
        phi = N.mod(phi - phi_0, 2*N.pi)
        return r_o if (phi > (2*N.pi-seg_phi) - phi_e) or (phi < phi_e) \
               else self.r_e/N.sin(phi)

    def in_pin(self, xyz):
        r, phi = self.pin.to_polar(xyz)
        r_o = self.r_o(phi)
        return r < r_o

class PointInCylindricalSections(object):
    def __init__(self, pins, r_es=None):
        if r_es is None: r_es = N.zeros(len(pins))
        self._in_pins = [PointInCylindricalSection(pin, r_e) for
                        pin, r_e in izip(pins, r_es)]

    def in_pins(self, xyz):
        for pin in self._in_pins:
            if pin.in_pin(xyz): return True
        return False

class GenPinGeo(object):
    """
    Generate pin geometry nodes/edges/faces suitable for export to a mesher
    (Gmsh mainly)
    """
    g_eps = 1e-10
    def __init__(self, pin, x3_t, x3_b=0.):
        self.pin = pin
        self.x3_t = x3_t ; self.x3_b = x3_b
    
    def calc_nodes(self):
        pin = self.pin
        self.no_segs = no_segs = int(N.ceil(pin.seg_phi/(N.pi/2)))
        self.d_phi = d_phi = pin.seg_phi/no_segs
        phis = pin.phi_0 - d_phi*N.arange(no_segs+1, dtype=N.float64)
        r = pin.r_p
        self.top_nodes = tn =N.array([pin.to_xyz(r, phi, self.x3_t) for phi in phis], N.float64)
        tn[N.abs(tn) < self.g_eps] = 0
        self.bot_nodes = bn = tn.copy()
        bn[:,pin.x3] = self.x3_b
        self.t_c_node = pin.to_xyz(0, 0, self.x3_t)
        self.b_c_node = pin.to_xyz(0, 0, self.x3_b)
        c_nodes = [self.t_c_node, self.b_c_node]
        self.nodes = N.vstack((tn, bn, c_nodes))
        self.t_node_nos = N.arange(len(tn))
        self.b_node_nos = N.arange(len(bn)) + len(tn)
        self.b_c_nodeno = len(self.nodes)-1
        self.t_c_nodeno = self.b_c_nodeno-1

    def calc_edges(self):
        tnn = self.t_node_nos ; bnn = self.b_node_nos
        tn_s = tnn[0] ; tn_e = tnn[-1]
        bn_s = bnn[0] ; bn_e = bnn[-1]
        horz_edges = N.array([[tn_s, tn_e], [bn_s, bn_e]], N.int32)
        vert_edges = N.array([[bn, tn] for tn, bn in zip(tnn, bnn)], N.int32)
        self.vert_edgenos = N.arange(len(vert_edges))
        self.horz_edgenos = N.arange(len(horz_edges)) + len(self.vert_edgenos)
        self.edges = N.vstack((vert_edges, horz_edges))
        
    def calc_segs(self):
        tnn = self.t_node_nos ; bnn = self.b_node_nos
        t_segs = N.array([[tnn[i], self.t_c_nodeno, tnn[i+1]]
                          for i in range(len(tnn)-1)], N.int32)
        b_segs = N.array([[bnn[i], self.b_c_nodeno, bnn[i+1]]
                          for i in range(len(bnn)-1)], N.int32)
        self.segs = N.vstack((t_segs, b_segs))

class GenPinsGeo(object):
    def __init__(self, pins, x3_t, x3_b=0.):
        self.pin_geom_gens = [GenPinGeo(pin, x3_t, x3_b)
                              for pin in pins]
    def calc_nodes(self):
        nds = []
        for pin_g in self.pin_geom_gens:
            pin_g.calc_nodes()
            nds.append(pin_g.nodes)
        self.pin_node_offs = N.cumsum([0] + [len(nd) for nd in nds[:-1]])
        self.nodes = N.vstack(nds)

    def calc_edges(self):
        edgs = []
        for i, pin_g in enumerate(self.pin_geom_gens):
            pin_g.calc_edges()
            edgs.append(pin_g.edges)
            edgs[-1] += self.pin_node_offs[i]
        self.edges = N.vstack(edgs)

    def calc_segs(self):
        sgs = []
        for i, pin_g in enumerate(self.pin_geom_gens):
            pin_g.calc_segs()
            sgs.append(pin_g.segs)
            sgs[-1] += self.pin_node_offs[i]
        self.segs = N.vstack(sgs)
