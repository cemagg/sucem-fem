# Authors:
# Evan Lezar <mail@evanlezar.com>

from dolfin import *

mesh = Mesh ( "sample_mesh.xml" )
V = FunctionSpace(mesh, "CG", 1)
interior_edges = MeshFunction ( 'uint', mesh, "sample_mesh_facet_region.xml" )

u = TrialFunction(V)
v = TestFunction(V)

E = Expression ( ("1.0", "0.0") )

for interior_edge in [100, 101, 102, 103]:
    n = V.cell().n
    l_p = dot( n('+'), E('+'))*dS(interior_edge)
    print interior_edge, assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )

for interior_edge in [100, 101, 102, 103]:
    n = V.cell().n
    l_p = dot( n('+'), E('-'))*dS(interior_edge)
    print interior_edge, assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )

plot( mesh, interactive=True)
