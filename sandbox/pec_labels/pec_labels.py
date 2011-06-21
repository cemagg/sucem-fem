__author__ = "Evan Lezar"
__date__ = "17 June 2011"

from dolfin import *
from dolfin_utils import meshconvert
import numpy as np
import os
        
def process_gmsh_elements ( ifilename, dim ):
    """
    Process the specified gmsh file and extract information on elements of dimension dim
    
    @param ifilename: The name of the gmsh file to read
    @param dim: The dimension of the elements to read
    """
    # Open files
    ifile = open(ifilename, "r")

    # Scan file for cell type
    dim_count = 0
    vertices_used = []
    # Arrays used to store gmsh tags for elements
    tags = []
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

            # Count number of elements for each dimension.
            # Also determine which vertices are not used.
            while line.find("$EndElements") == -1:
                element = line.split()
                elem_type = int(element[1])
                
                if elem_type == dim:
                    num_tags = int(element[2])
                    node_num_list = [int(node)-1 for node in element[3 + num_tags:]]
                    
                    vertices_used.extend ( node_num_list )
                    if num_tags > 0:
                        tags.append(tuple(int(tag) for tag in element[3:3+num_tags]))
                    dim_count += 1
                
                line = ifile.readline()
            break;
        else:
            # Read next line
            line = ifile.readline()
    
    ifile.close()
    
    return dim_count, vertices_used, tags

def convert_and_create_facet_mesh_function ( ifilename, ofilename ):
    # First convert the gmsh mesh
    meshconvert.convert2xml ( ifilename, ofilename )
    
    # Now load the created mesh and initialise the required connectivity information
    mesh = Mesh ( ofilename )
    mesh.order()
    D = mesh.topology().dim()
    mesh.init(D-1, 0)
    
    # read the data from the gmsh file once again
    dim_count, vertices_used, tags = process_gmsh_elements( ifilename, D-1 )
    # Get the facet-node connectivity information (reshape as a row of node indices per facet)
    facets_as_nodes = mesh.topology()(D-1,0)().reshape ( mesh.num_facets(), D )
    
    # Create and initialise the mesh function
    facet_mark_function = MeshFunction ( 'uint', mesh, D-1 )
    facet_mark_function.set_all( 0 )
    
    # set the relevant values of the mesh function
    facets_to_check = range( mesh.num_facets() )
    for i in range(len(tags)):
        nodes = np.sort(np.array(vertices_used[2*i:(2*i+D)]))
        value  = tags[i][0]
        
        if value != 0:
            found = False
            for j in range(len(facets_to_check)):
                index = facets_to_check[j]
                if np.array_equal(facets_as_nodes[index,:], nodes):
                    found = True;
                    facets_to_check.pop(j)
                    # set the value of the mesh function
                    facet_mark_function[index] = value
                    break;
                
            if not found:
                raise Exception ( "The facet (%d) was not found to mark: %s" % (i, nodes) )
        
    # save the mesh function to file
    fname = os.path.splitext(ofilename)[0]
    mesh_function_file = File("%s_%s.xml" % (fname, "facet_regions"))
    
    mesh_function_file << facet_mark_function

if __name__ == "__main__":
    convert_and_create_facet_mesh_function ( "untitled.msh", "untitled.xml" )
    