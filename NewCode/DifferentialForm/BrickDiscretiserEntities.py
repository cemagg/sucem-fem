import numpy as N
from NewCode import ProxyList
from NewCode.Meshes  import BrickMesh
from NewCode import Mesh
from NewCode.DifferentialForm import DiscretiserEntities 
geom_entity_names = DiscretiserEntities.geom_entity_names

class DiscretiserEntityList(DiscretiserEntities.BaseDiscretiserEntityList,
                            Mesh.BoundaryProxyEntities): pass

class DiscretiserEntity(DiscretiserEntities.DiscretiserEntity): pass

class Edge(DiscretiserEntity, BrickMesh.Edge):
    entity_type = 'edge'
    proxy_attrs = BrickMesh.Edge.proxy_attrs + \
                  DiscretiserEntity.proxy_attrs
    
class Face(DiscretiserEntity, BrickMesh.Face):
    entity_type = 'face'
    proxy_attrs = BrickMesh.Face.proxy_attrs + \
                  DiscretiserEntity.proxy_attrs

class make_geomEntities(DiscretiserEntities.make_geomEntities):
    discretiserEntities = dict(edge=Edge, face=Face)
