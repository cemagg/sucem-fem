from __future__ import division

from itertools import izip, chain
import numpy as N
from scipy import sparse

from NewCode.Utilities import Struct, CacheLast
from NewCode.DifferentialForm import allconstrained

__all__ = ('gen_group_freefuns', 'gen_elgroups')

def gen_group_freefuns(mesh, elgroups, global_freefun):
    """Generate freedom-functions for the various implicit/explicit DOF groups

    Input Params
    ============
    mesh -- Mesh of the geometry
    elgroups -- Implict/Explicit element groups as calculated by gen_elgroups
    global_freefun -- Freedom function for whole problem, will be ANDed

    Output Value
    ============
    Freedom function tuples (general, volume) in Struct with keys
    a -- DOFs interior to implicit region
    b -- Fully implicit DOFs in implicit/explicit hybrid element
    c -- Imp/Explicit hybrid DOFs on bdry between implicit and explicit region
    d -- Fully explicit DOFs in explicit/implicit hybrid element
    e -- DOFs interior to explicit region
    """
    ii = elgroups.ii ; ie = elgroups.ie ; ei = elgroups.ei ; ee = elgroups.ee
    iiie = ii|ie ;  eiee = ei|ee
    constr= set([-1])
    def free_b(ent):
        connect2 = set(ent.connect2elem) - constr
        return (global_freefun(ent) and not 
                connect2 - iiie and connect2 - ii)

    def free_c(ent):
        connect2 = set(ent.connect2elem) - constr
        return (global_freefun(ent) and connect2 & ie and connect2 & ei)

    def free_d(ent):
        connect2 = set(ent.connect2elem) - constr
        return (global_freefun(ent) and not 
                connect2 - eiee - constr and connect2 - ee)

    return Struct(
        a=(lambda ent: (global_freefun(ent) and not 
                        set(ent.connect2elem) - ii - constr),
           lambda ent: ent.index in ii), 
        b=(free_b, lambda ent: ent.index in ie),
        c=(free_c, allconstrained),
        d=(free_d, lambda ent: ent.index in ei),
        e=(lambda ent: (global_freefun(ent) and not
                        set(ent.connect2elem) - ee - constr),
           lambda ent: ent.index in ee))

def gen_elgroups(mesh, imp_elset):
    """
    Generate implict, implicit/explicit boundary and explicit element groups

    Input Params
    ============

    mesh -- Mesh of the geometry
    imp_elset -- set of the element numbers that are implicit. Others assumed
                 explicit

    Output Value
    ============
    Struct with keys
    ii -- Implicit elements interior to implicit region
    ie -- Implicit elements connected to explicit element(s) at
          implicit/explicit boundary
    ei -- Explicit elements connected to implicit element(s) in ie
    ee -- Explicit element interior to explicit region
    """
    all_elset = set(range(len(mesh.elements)))
    exp_elset = all_elset - imp_elset
    ii=set(); ie=set(); ei=set(); ee=set()
    for elno, el in enumerate(mesh.elements):
        if elno in imp_elset:
            if not N.any([ce in exp_elset for ce in el.connect2elem]):
                ii.add(elno)
            else: ie.add(elno)
        else:
            if not N.any([ce in imp_elset for ce in el.connect2elem]):
                ee.add(elno)
            else: ei.add(elno)
    return Struct(ii=ii, ie=ie, ei=ei, ee=ee)

    
def gen_hybmesh_group_freefuns(meshes, on_hbdry, global_freefun,
                               elgroups=None, dead_elements=None):
    """Generate freedom-functions for the various implicit/explicit DOF groups

    Input Params
    ============
    meshes -- Struct with attrs imp -> implicit mesh, exp -> explicit mesh
    global_freefun -- Freedom function for whole problem, will be ANDed
    on_hbdry -- Function; returns True for geom entities on hybrid boundary
    elgroups -- Optional saved value of ie_set and ei_set

    Output Value
    ============
    (FreeTuplesStruct, elgroups)

    FreeTuplesStruct contains tuples (general, volume) in Struct with keys
    a -- DOFs interior to implicit region
    b -- Fully implicit DOFs in implicit/explicit hybrid element
    c -- Imp/Explicit hybrid DOFs on bdry between implicit and explicit region
    d -- Fully explicit DOFs in explicit/implicit hybrid element
    e -- DOFs interior to explicit region
    """
    if not elgroups:
        ie_set = set(chain(*(e.connect2elem for e in meshes.imp.edges
                             if on_hbdry(e))))
        ei_set = set(chain(*(e.connect2elem for e in meshes.exp.edges
                             if on_hbdry(e))))
    else:
        ie_set, ei_set = elgroups.ie, elgroups.ei
    ie_set.discard(-1) ; ei_set.discard(-1)
    if dead_elements is None:
        e_ent = lambda ent: (global_freefun(ent) and not
                             set(ent.connect2elem) & ei_set)
        e_vol = lambda vol: vol.index not in ei_set
    else:
        def e_ent(ent):
            conset = set(ent.connect2elem)
            return (global_freefun(ent) and not
                    conset & ei_set and not
                    conset & dead_elements)
        def e_vol(vol):
            return (vol.index not in ei_set and
                    vol.index not in dead_elements)

    return (Struct(a=(lambda ent: (global_freefun(ent) and not
                                   set(ent.connect2elem) & ie_set),
                      lambda vol: vol.index not in ie_set),
                   b=(lambda ent: (global_freefun(ent) and not on_hbdry(ent) and
                                   bool(set(ent.connect2elem) & ie_set)),
                      lambda vol: vol.index in ie_set),
                   c=(lambda ent: (global_freefun(ent) and on_hbdry(ent)),
                                        lambda vol: False),
                   d=(lambda ent: (global_freefun(ent) and not on_hbdry(ent) and
                                   bool(set(ent.connect2elem) & ei_set)),
                      lambda vol: vol.index in ei_set),
                   e=(e_ent, e_vol)),
            Struct(ei=ei_set, ie=ie_set))
            
