import numpy as N
from NewCode import ProxyList
from NewCode.Meshes  import PyramMesh
from NewCode import Mesh
from NewCode.DifferentialForm import DiscretiserEntities as DE

class DiscretiserEntityList(DE.DiscretiserEntityListNoBoundary): pass
class DiscretiserEntity(DE.DiscretiserEntity): pass

class Edge(DiscretiserEntity, PyramMesh.Edge):
    proxy_attrs = PyramMesh.Edge.proxy_attrs + \
                  DiscretiserEntity.proxy_attrs
    
class ApexFace(DiscretiserEntity, PyramMesh.ApexFace):
    proxy_attrs = PyramMesh.ApexFace.proxy_attrs + \
                  DiscretiserEntity.proxy_attrs

class BaseFace(DiscretiserEntity, PyramMesh.BaseFace):
    proxy_attrs = PyramMesh.BaseFace.proxy_attrs + \
                  DiscretiserEntity.proxy_attrs

class make_geomEntities(DE.make_geomEntities):
    discretiserEntities = dict(
        edge=Edge, apexface=ApexFace, baseface=BaseFace, vol=DE.FakeVolEntity)
