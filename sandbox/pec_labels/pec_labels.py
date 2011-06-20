__author__ = "Evan Lezar"
__date__ = "17 June 2011"

from dolfin_utils_meshconvert import convert2xml

from dolfin import *
import numpy as np

class Domain(SubDomain):
    def inside (self, x, on_boundary):
        return on_boundary

def test_mesh_function():
    
    mesh = UnitSquare ( 2, 2 )
    
    boundary = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
    boundary.set_all ( 0 )
    bn = Domain()
    
    bn.mark(boundary, 2 )
    
    file = File("boundary_function.xml")
    file << boundary
        
def process_gmsh ( ifilename = "untitled.msh"):
    
    # Open files
    ifile = open(ifilename, "r")

    # Scan file for cell type
    dim_count = {1: 0, 2: 0, 4: 0}
    vertices_used = {1: [], 2: [], 4: []}
    # Arrays used to store gmsh tags for elements
    tags = {1: [], 2: [], 4: []}
    line = ifile.readline()
    while line:

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        # Read dimension
        if line.find("$Elements") == 0:

            line = ifile.readline()
            num_cells  = int(line)
            num_cells_counted = 0
            if num_cells == 0:
                _error("No cells found in gmsh file.")
            line = ifile.readline()

            # Now iterate through elements to find largest dimension.  Gmsh
            # format might include elements of lower dimensions in the element list.
            # We also need to count number of elements of correct dimensions.
            # Also determine which vertices are not used.
            while line.find("$EndElements") == -1:
                element = line.split()
                elem_type = int(element[1])
                num_tags = int(element[2])
                node_num_list = [int(node)-1 for node in element[3 + num_tags:]]
                
                vertices_used[elem_type].extend ( node_num_list )
                if num_tags > 0:
                    tags[elem_type].append(tuple(int(tag) for tag in element[3:3+num_tags]))
                dim_count[elem_type] += 1
                
                line = ifile.readline()
            break;
        else:
            # Read next line
            line = ifile.readline()
    
    
    
    print dim_count
    print vertices_used
    print tags
    
    ifile.close()
    
    return dim_count, vertices_used, tags

def load_and_construct ():
    mesh = Mesh ( "untitled.xml" )
    mesh.order()
    D = mesh.topology().dim()
    
    dim_count, vertices_used, tags = process_gmsh () 
    
    mesh.init(D-1, 0)
    facets_as_nodes = mesh.topology()(D-1,0)().reshape(mesh.num_facets(), D)
    
    num_facets_to_mark = len(tags[D-1])
    index_and_value_of_facet_to_mark = []
    
    for i in range(num_facets_to_mark):
        nodes = np.sort(np.array(vertices_used[D-1][2*i:(2*i+D)]))
        value  = tags[D-1][i][0]
        
        found = False
        for j in range ( mesh.num_facets() ):
            if np.array_equal(facets_as_nodes[j,:], nodes):
                index = j
                found = True;
                break;
            
        if not found:
            raise Exception ( "The facet (%d) to mark was not found: %s" % (i, nodes) )
        
        index_and_value_of_facet_to_mark.append ( (index, value) )
        
    facet_mark_function = MeshFunction ( 'uint', mesh, D-1 )
    facet_mark_function.set_all( 0 )
    
#    facet_mark_function_values = np.zeros( facet_mark_function.size(), dtype=np.uint32 )
    
    for i_v in index_and_value_of_facet_to_mark:
        facet_mark_function[i_v[0]] = i_v[1]
        
#    facet_mark_function.set ( facet_mark_function_values.tolist() )
    
    print facet_mark_function.str(True)
    
    plot ( facet_mark_function, interactive = True )
    
#    plot ( mesh, interactive=True)
    
    print "Done"



if __name__ == "__main__":
    convert2xml ( "untitled.msh", "untitled.xml" )
    load_and_construct()
    