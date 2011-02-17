import numpy as N
from NewCode import SubDimMesh, ProxyList
from NewCode.DifferentialForm import DiscretiserEntities
geom_entity_names = ('node', 'edge', 'face')

class SubDimDiscretiserEntity(DiscretiserEntities.DiscretiserEntity):
    pass

class SubDimDiscretiserEntityList(DiscretiserEntities.DiscretiserEntityList):
    pass

class Edge(SubDimDiscretiserEntity, SubDimMesh.Edge):
    proxy_attrs = SubDimMesh.Edge.proxy_attrs + SubDimDiscretiserEntity.proxy_attrs

class Face(SubDimDiscretiserEntity, SubDimMesh.FaceElement):
    proxy_attrs = SubDimMesh.FaceElement.proxy_attrs + SubDimDiscretiserEntity.proxy_attrs

