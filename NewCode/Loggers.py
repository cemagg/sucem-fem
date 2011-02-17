from __future__ import division

from itertools import izip
import pickle
import numpy as N

from NewCode.Utilities import partial, close_to_point, Struct, \
     in_box_vec, on_box_surf

class NTFFCollecter(object):
    g_eps = 1e-10
    def __init__(self, discs, Gamma_C_corners, g_eps=None):
        if g_eps is not None: self.g_eps = g_eps
        self.discs = discs
        self.mesh = discs.E.mesh
        assert(self.mesh is discs.B.mesh)
        self.Gamma_C_corners = Gamma_C_corners
        self.on_Gamma_Cp = on_box_surf(Gamma_C_corners[0], Gamma_C_corners[1], self.g_eps)
        self.in_Gamma_Cp = in_box_vec(*Gamma_C_corners)
        self.C_faces = C_faces = set(
            [f.index for f in self.mesh.faces if self.on_Gamma_Cp(f.nodeCoords)])
        self.Gamma_C_fs = lambda face: face.index in C_faces
        self.init_logentities()
        self.init_dofnos()

    def init_logentities(self):
        f_con_els = self.mesh.faces[list(self.C_faces)].connect2elem.flatten()
        self.f_con_els = f_con_els = N.unique(f_con_els[f_con_els != -1])
        # Elements connected to the collection surface  that have mid-pts inside it.
        log_elements = self.f_con_els[self.in_Gamma_Cp(N.average(
            self.mesh.elements[self.f_con_els].nodeCoords, axis=1))]
        log_edges = N.unique(self.mesh.elements[log_elements].edgenos.flatten())
        log_faces = N.unique(self.mesh.elements[log_elements].facenos.flatten())
        self.log_ents = Struct(vol=log_elements, face=log_faces, edge=log_edges)
        
    def init_dofnos(self):
        globtables = Struct()
        B_perm = self.discs.B.permuter
        B_globtables = Struct(
            (ent, B_perm.globalEntityPermutationTable(ent)[self.log_ents[ent]])
            for ent in self.discs.B.geomEntities)
        E_perm = self.discs.E.permuter
        E_globtables = Struct(
            (ent, E_perm.globalEntityPermutationTable(ent)[self.log_ents[ent]])
            for ent in self.discs.E.geomEntities)
        self.globtables = Struct(E=E_globtables, B=B_globtables)
        for quant, doftbl in self.globtables.items():
            for ent, entdoftbl in doftbl.items():
                assert(not N.any(entdoftbl.flatten() < 0))
        self.logdofs = Struct((quant, N.hstack(
            [self.globtables[quant][ent].flatten() for ent
             in sorted(self.globtables[quant].keys())]))
                              for quant in ('E', 'B'))
            
