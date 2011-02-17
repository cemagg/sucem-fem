import numpy as N
from NewCode import BrickSubDimMesh, ProxyList
from NewCode.DifferentialForm import BrickDiscretiserEntities
geom_entity_names = ('node', 'edge', 'face')

class SubDimDiscretiserEntity(BrickDiscretiserEntities.DiscretiserEntity):
    pass

class SubDimDiscretiserEntityList(BrickDiscretiserEntities.DiscretiserEntityList):
    pass

class Edge(SubDimDiscretiserEntity, BrickSubDimMesh.Edge):
    proxy_attrs = BrickSubDimMesh.Edge.proxy_attrs + SubDimDiscretiserEntity.proxy_attrs

class Face(SubDimDiscretiserEntity, BrickSubDimMesh.FaceElement):
    proxy_attrs = BrickSubDimMesh.FaceElement.proxy_attrs \
                  + SubDimDiscretiserEntity.proxy_attrs
