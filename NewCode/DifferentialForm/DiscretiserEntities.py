import numpy as N
from NewCode import Mesh, ProxyList
from NewCode.Utilities import Struct
from NewCode.DifferentialForm import allfree

geom_entity_names = ('node', 'edge', 'face', 'vol')


class DiscretiserEntity(object):
    proxy_attrs = ('freeNo',)

    @property
    def freeNo(self):
        """
        The 'freedom index' of this entity. If this entity is constrained, -1.

        """
        return self._freeNo[self.index]

    @property
    def isFree(self):
        """
        Returns True if the entity is free, False if it is constrained.
        """
        return self._freeNo[self.index] > -1

    def __init__(self, *names, **kwargs):
        """
        The mesh kwarg is used only to store the mesh identity

        The mesh kwarg is used as a consistency check when several mesh-related objects
        are working togheter to ensure they are all related to the same mesh.
        """
        super(DiscretiserEntity, self).__init__(*names, **kwargs)
        self.mesh = kwargs['mesh']
        try: self.freefun = kwargs['freefun']
        except: self.freefun = None

class BaseDiscretiserEntityList(object):

    def __init__(self, *names, **kwargs):
        super(BaseDiscretiserEntityList, self).__init__(*names, **kwargs)
        self.mesh = self.entity.mesh
        freefun = self.entity.freefun
        # -1 represents a constrained entity
        self.entity._freeNo = N.zeros(len(self), int) -1
        _freeNo = self.entity._freeNo
        freeNo = 0
        for index, entity in enumerate(self):
            if freefun(entity):
                _freeNo[index] = freeNo
                freeNo += 1
        self.noFree = freeNo

class DiscretiserEntityListNoBoundary(BaseDiscretiserEntityList,
                                      ProxyList.ProxyList): pass

class DiscretiserEntityList(BaseDiscretiserEntityList,
                            Mesh.BoundaryProxyEntities): pass

class Edge(DiscretiserEntity, Mesh.Edge):
    entity_type = 'edge'
    proxy_attrs = Mesh.Edge.proxy_attrs + DiscretiserEntity.proxy_attrs

class Face(DiscretiserEntity, Mesh.Face):
    entity_type = 'face'
    proxy_attrs = Mesh.Face.proxy_attrs + DiscretiserEntity.proxy_attrs

class FakeVolEntity(object):
    def __init__(self, freefun=None):
        self.entity = Struct(freefun=freefun)
        
class make_geomEntities(Struct):
    discretiserEntities = dict(edge=Edge,face=Face)
    def __init__(self, mesh, basisSet, freeFun, vol_freeFun=None, el_in_disc=None):
        Struct.__init__(self)
        self.el_in_disc = el_in_disc if el_in_disc else allfree
        del self['el_in_disc']
        self.update(dict((ent_name, DiscretiserEntityList(
            self.discretiserEntities[ent_name](
            attrs=getattr(mesh, ent_name+'s').list_repr(partial=True),
            mesh=mesh, freefun=freeFun))) 
                         for ent_name in basisSet.fns if ent_name != 'vol'))
        # Note that volumes have to be handled by the discretiser object itself
        if 'vol' in basisSet.fns: self['vol'] = FakeVolEntity(vol_freeFun)
        
    

