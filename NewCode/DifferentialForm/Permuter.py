import numpy as N
from NewCode.Exceptions import AllZero, DeadElement

from NewCode.Utilities import Struct

class Permuter(object):
    """
    Algorithm
    =========

    Local Numbering
    ---------------

    The local numbering has to match the ordering used by the
    DiscretiserElement class. This means that all edge basis functions are
    numbered first, then face and finally volume. Each basis function is
    numbered for all entities before the next basis function is
    considrered. For an hypothetical element with 1 basisfunction per edge and
    six edges, 2 per face and four faces and three per volume, the local
    numbering would be:

    [e1fn1, e2fn1, e3fn1, e4fn1, e5fn1, e6fn1, 
     f1fn1, f2fn1, f3fn1, f4fn1, f1fn2, f2fn2, f3fn2, f4fn2,
     v1fn1, v1fn2, v1fn3
    ]

    Global Numbering
    ----------------

    Numbering for each geometric entity type is handled seperately. I.e. all
    the edge dofs are numbered from 0 thru no_edgedofs, face functions are
    numbered from no_edgedofs+1 thru no_edgedofs + no_facedofs etc.

    Each entity is fully numbered, in order of basis-function. I.e., for M
    global edges with 3 dofs per edge, the dofs are numbered

    [e1fn1, e1fn2, e1fn3, e2fn1, e2fn2, e2fn3, ..., eMfn1, eMfn2, eMfn3]
    
    """
    # Defines the order of the geometric entities, as well which are known
    def __init__(self, parentElement, geomEntities):
        for ent_name in geomEntities:
            assert(ent_name in parentElement.entityNames)
        self.geomEntities = geomEntities

        # Number of DOFs for each geometric entity type
        total_dofs = lambda ent_name,ent_obj: \
            ent_obj.noFree*parentElement.noDOFs[ent_name].perEntity
        self.noDOFs = dict((ent_name, total_dofs(ent_name, ent_obj))
                           for ent_name, ent_obj in geomEntities.items())
        ent_names = tuple(ent_name for ent_name in parentElement.entityNames
                          if ent_name in parentElement.noDOFs)
        self.noDOFsPerEntity = dict((en, parentElement.noDOFs[en].perEntity)
                                     for en in ent_names)
        ent_offsets = [0] + N.cumsum(tuple(
            self.noDOFs[ent_name] for ent_name in ent_names)).tolist()
        # Global DOF Offset for each geometric entity type
        self.offsetDOFs = dict((ent_name, ent_offs) for
                               ent_name, ent_offs in zip(ent_names, ent_offsets))
        # The total number of DOFs in the system
        self.totalDOFs = ent_offsets[-1]
        # Geometrical entities that are actually present, in the correct order
        self.presentEntities = ent_names

    def globalEntityPermutationTable(self, ent_name):
        dofs_per = self.noDOFsPerEntity[ent_name]
        ent_offset = self.offsetDOFs[ent_name]
        constrained = N.ones(dofs_per, N.int32)*-1
        r = N.arange(dofs_per)
        return N.array([r+fn*dofs_per+ent_offset if fn >= 0 else constrained
                        for fn in self.geomEntities[ent_name][:].freeNo],
                       N.int32)

    def partialGlobalEntityPermutationTable(self, ent_name, ent_nos):
        dofs_per = self.noDOFsPerEntity[ent_name]
        ent_offset = self.offsetDOFs[ent_name]
        constrained = N.ones(dofs_per, N.int32)*-1
        r = N.arange(dofs_per)
        return N.array([r+fn*dofs_per+ent_offset if fn >= 0 else constrained
                        for fn in self.geomEntities[ent_name][ent_nos].freeNo],
                       N.int32)
    
    def permuteElementPerEntity(self, element):
        def number_dofs(ent_name, ent_obj):
            dofs_per = element.noDOFs[ent_name].perEntity
            freedom_nos = ent_obj[getattr(element, ent_name+'nos')].freeNo
            return N.hstack(freedom_nos*dofs_per+i
                            for i in range(dofs_per)) 

        return dict((ent_name, number_dofs(ent_name, ent_obj))
                    for ent_name, ent_obj in self.geomEntities.iteritems())

    def permuteElementWithConstrained(self, element):
        ent_dofs = self.permuteElementPerEntity(element)
        dofs = [] #; dof_offs = 0
        for ent_name in self.presentEntities:
            these_dofs = ent_dofs[ent_name]
            dof_offs = self.offsetDOFs[ent_name]
            dofs.append(N.where(these_dofs >= 0,these_dofs+dof_offs,these_dofs))
        return dofs

    def permuteElement(self, element, remove_constrained=True):
        """
        Return a tuple with two arrays (local_dofnos, global_dofnos)

        Constrained DOFs are removed from the two dofno arrays.
        """
        if element.index in element.dead_elements: raise DeadElement
        g_perm = N.hstack(self.permuteElementWithConstrained(element))
        l_perm = N.arange(element.noDOFs.element)
        if remove_constrained:
            l_perm = l_perm[g_perm >= 0]
            if len(l_perm) == 0: raise AllZero
            g_perm = g_perm[g_perm >= 0]
        return (l_perm, g_perm)

    def permuteElementEntities(self, element, remove_constrained=True):
        def number_dofs(ent_name, ent_obj):
            dofs_per = element.noDOFs[ent_name].perEntity
            freedom_nos = ent_obj[getattr(element, ent_name+'nos')].freeNo
            r = N.arange(dofs_per)
            try: return N.vstack(fn*dofs_per+r for fn in freedom_nos)
            except TypeError: return N.array([freedom_nos*dofs_per+r])

        g_ent_perms = dict((ent_name, number_dofs(ent_name, ent_obj))
                            for ent_name, ent_obj in self.geomEntities.iteritems())
        ent_perms = Struct()
        for ent_name, g_ent_perm in g_ent_perms.iteritems():
            ent_obj = self.geomEntities[ent_name]
            lf = ent_obj[getattr(element, ent_name+'nos')].isFree
            if ent_name is not 'vol':
                local_free_ent = N.arange(len(lf))[lf == True]
            else: local_free_ent = 0
            if remove_constrained:
                ent_perms[ent_name] = [
                    local_free_ent,
                    # If the first edge function is constrained, the rest of the
                    # functions on that edge are too. This is assumption is needed
                    # because a 2-D truth-table in [] results in the array being
                    # flattened
                    g_ent_perm[g_ent_perm[:,0] >= 0] + self.offsetDOFs[ent_name]]
            else:
                g_ent_perm[g_ent_perm[:,0] >= 0] += self.offsetDOFs[ent_name]
                ent_perms[ent_name] = g_ent_perm
        
        return ent_perms

    def newDOFs(self, dtype=N.float64):
        return N.zeros(dtype=dtype, shape=self.totalDOFs)
