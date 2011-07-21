# Authors:
# Evan Lezar <mail@evanlezar.com>

from dolfin import *

def use_saved_mesh():
    print 'Saved Mesh:'
    sides = ['left', 'top', 'right', 'bottom']
    
    mesh = Mesh ( "sample_mesh.xml" )
    V = FunctionSpace(mesh, "CG", 1)
    interior_edges = MeshFunction ( 'uint', mesh, "sample_mesh_facet_region.xml" )
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    E = Expression ( ("1.0", "0.0") )
    
    for interior_edge in [100, 101, 102, 103]:
        n = V.cell().n
        l_p = dot( n('+'), E('+'))*dS(interior_edge)
        print sides[interior_edge-100], assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )
    print '----'
    for interior_edge in [100, 101, 102, 103]:
        n = V.cell().n
        l_p = dot( n('+'), E('-'))*dS(interior_edge)
        print sides[interior_edge-100], assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )
    
    plot( mesh, interactive=True)

def use_generated_mesh():
    """construct a 1x1 mesh that has inside it a 0.5x0.5 square mesh.
    
    Integrate the dot product of the facet normal with a unit field in
    the x direction along each one of the sides of the interior square
    (0.5x0.5).

    The integral of the two vertical sides should equal the
    length of each side.
    """
    print 'Generated Mesh:'
    class Left ( SubDomain ):
        def inside (self, x, on_boundary ):
            if x[0] == 0.25 and ( x[1] >= 0.25 and x[1] <= 0.75 ):
                return True
            return False
    
    class Right ( SubDomain ):
        def inside (self, x, on_boundary ):
            if x[0] == 0.75 and ( x[1] >= 0.25 and x[1] <= 0.75 ):
                return True
            return False
    
    class Top ( SubDomain ):
        def inside (self, x, on_boundary ):
            if x[1] == 0.75 and ( x[0] >= 0.25 and x[0] <= 0.75 ):
                return True
            return False
    
    class Bottom ( SubDomain ):
        def inside (self, x, on_boundary ):
            if x[1] == 0.25 and ( x[0] >= 0.25 and x[0] <= 0.75 ):
                return True
            return False
    
    sides = ['left', 'top', 'right', 'bottom']
    
    mesh = UnitSquare ( 8, 8 )
    V = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    # a uniform x-directed field
    E = Expression ( ("1.0", "0.0") )
    
    interior_edges = MeshFunction ( 'uint', mesh, 1 )
    interior_edges.set_all ( 0 )
    
    left_side = Left ( )
    left_side.mark( interior_edges, 100 )
    top_side = Top ( )
    top_side.mark( interior_edges, 101 )
    right_side = Right ( )
    right_side.mark( interior_edges, 102 )
    bottom_side = Bottom ( )
    bottom_side.mark( interior_edges, 103 )
    
    for interior_edge in [100, 101, 102, 103]:
        n = V.cell().n
        l_p = dot( n('+'), E('+'))*dS(interior_edge)
        print sides[interior_edge-100], assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )
    
    print '----'
    for interior_edge in [100, 101, 102, 103]:
        n = V.cell().n
        l_p = dot( n('+'), E('-'))*dS(interior_edge)
        print sides[interior_edge-100], assemble ( l_p, interior_facet_domains=interior_edges, mesh=mesh )
    
#    plot( mesh, interactive=True)


if __name__ == "__main__":
    use_generated_mesh()
    print '----'
    use_saved_mesh()
