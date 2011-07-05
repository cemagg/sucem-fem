# Authors:
# Evan Lezar <mail@evanlezar.com>
__date__ = "17 June 2011"

from dolfin import *
from dolfin_utils import meshconvert
import numpy as np
import os

def gmsh_to_dolfin_mesh ( ifilename, handler ):
    """Convert between .gmsh v2.0 format (http://www.geuz.org/gmsh/) and .xml,
    parser implemented as a state machine:

        0 = read 'MeshFormat'
        1 = read  mesh format data
        2 = read 'EndMeshFormat'
        3 = read 'Nodes'
        4 = read  number of vertices
        5 = read  vertices
        6 = read 'EndNodes'
        7 = read 'Elements'
        8 = read  number of cells
        9 = read  cells
        10 = done

    Afterwards, extract physical region numbers if they are defined in
    the mesh file as a mesh function.

    """

    print "Converting from Gmsh format (.msh, .gmsh) to DOLFIN XML format"
    
    # The dimension of the gmsh element types supported here as well as the dolfin cell types for each dimension 
    gmsh_dim = {1: 1, 2: 2, 4: 3}
    gmsh_cell_type = {1: "interval", 2: "triangle", 3: "tetrahedron" }
    # the gmsh element types supported for conversion
    supported_gmsh_element_types = [1, 2, 4]
    
    # Open files
    ifile = open(ifilename, "r")

    # Scan file for cell type
    cell_type = None
    highest_dim = 0
    line = ifile.readline()
    while line:

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        # Read dimension
        if line.find("$Elements") == 0:

            line = ifile.readline()
            num_elements = int(line)
            num_cells_counted = 0
            if num_elements == 0:
                _error("No elements found in gmsh file.")
            line = ifile.readline()

            # Now iterate through elements to find largest dimension.  Gmsh
            # format might include elements of lower dimensions in the element list.
            # We also need to count number of elements of correct dimensions.
            # Also determine which vertices are not used.
            dim_count = {1: 0, 2: 0, 3: 0}
            vertices_used = {1: [], 2: [], 3: []}
            # Array used to store gmsh tags for 1D (type 1/line), 2D (type 2/triangular) elements and 3D (type 4/tet) elements
            tags_for_dim = {1: [], 2: [], 3: []}
            
            while line.find("$EndElements") == -1:
                element = line.split()
                elem_type = int(element[1])
                num_tags = int(element[2])
                
                if elem_type in supported_gmsh_element_types:
                    dim = gmsh_dim[elem_type]
                    if highest_dim < dim:
                        highest_dim = dim
                    
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_used[dim].extend(node_num_list)
                    if num_tags > 0:
                        tags_for_dim[dim].append(tuple(int(tag) for tag in element[3:3+num_tags]))    
                    dim_count[dim] += 1
                else:
                    #TODO: output a warning here. "gmsh element type %d not supported" % elem_type
                    pass

                line = ifile.readline()
        else:
            # Read next line
            line = ifile.readline()

    # Check that we got the cell type and set num_cells_counted
    if highest_dim == 0:
        _error("Unable to find cells of supported type.")
    
    num_cells_counted = dim_count[highest_dim]
    vertex_set = set(vertices_used[highest_dim])
    vertices_used[highest_dim] = None

    vertex_dict = {}
    for n,v in enumerate(vertex_set):
        vertex_dict[v] = n

    # Step to beginning of file
    ifile.seek(0)
    
    # Set mesh type
    handler.set_mesh_type(gmsh_cell_type[highest_dim], highest_dim)

    # Initialise node list (gmsh does not export all vertexes in order)
    nodelist = {}

    # Current state
    state = 0

    # Write data
    num_vertices_read = 0
    num_cells_read = 0

    # Now handle the facet markings
    if len(tags_for_dim[highest_dim-1]) > 0:
        # first construct the mesh
        from dolfin import MeshEditor, Mesh
        mesh = Mesh()
        me = MeshEditor ()
        me.open( mesh, highest_dim, highest_dim )
    else:
        me = None
    

    while state != 10:
        
        # Read next line
        line = ifile.readline()
        if not line: break

        # Skip comments
        if line[0] == '#':
            continue

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        if state == 0:
            if line == "$MeshFormat":
                state = 1
        elif state == 1:
            (version, file_type, data_size) = line.split()
            state = 2
        elif state == 2:
            if line == "$EndMeshFormat":
                state = 3
        elif state == 3:
            if line == "$Nodes":
                state = 4
        elif state == 4:
            num_vertices = len(vertex_dict)
            handler.start_vertices(num_vertices)
            if me is not None:
                me.init_vertices ( num_vertices )
            state = 5
        elif state == 5:
            (node_no, x, y, z) = line.split()
            node_no = int(node_no)
            x,y,z = [float(xx) for xx in (x,y,z)]
            if vertex_dict.has_key(node_no):
                node_no = vertex_dict[node_no]
            else:
                continue
            nodelist[int(node_no)] = num_vertices_read
            handler.add_vertex(num_vertices_read, [x, y, z])
            if me is not None:
                if highest_dim == 1:
                    me.add_vertex( num_vertices_read, x)
                elif highest_dim == 2:
                    me.add_vertex( num_vertices_read, x, y)
                elif highest_dim == 3:
                    me.add_vertex( num_vertices_read, x, y, z)
                    
            num_vertices_read +=1

            if num_vertices == num_vertices_read:
                handler.end_vertices()
                state = 6
        elif state == 6:
            if line == "$EndNodes":
                state = 7
        elif state == 7:
            if line == "$Elements":
                state = 8
        elif state == 8:
            handler.start_cells(num_cells_counted)
            if me is not None:
                me.init_cells( num_cells_counted )
                
            state = 9
        elif state == 9:
            element = line.split()
            elem_type = int(element[1])
            num_tags  = int(element[2])
            if elem_type in supported_gmsh_element_types:
                dim = gmsh_dim[elem_type]
            else:
                dim = 0
            if dim == highest_dim:
                node_num_list = [vertex_dict[int(node)] for node in element[3 + num_tags:]]
                for node in node_num_list:
                    if not node in nodelist:
                        _error("Vertex %d of %s %d not previously defined." %
                              (node, gmsh_cell_type[dim], num_cells_read))
                cell_nodes = [nodelist[n] for n in node_num_list]
                handler.add_cell(num_cells_read, cell_nodes)
                
                if me is not None:
                    me.add_cell( num_cells_read, *cell_nodes )
                    
                num_cells_read +=1

            if num_cells_counted == num_cells_read:
                handler.end_cells()
                if me is not None:
                    me.close()
                state = 10
        elif state == 10:
            break
    
    
            
    # Write mesh function based on the Physical Regions defined by
    # gmsh, but only if they are not all zero. All zero physical
    # regions indicate that no physical regions were defined.
    if highest_dim not in [1,2,3]:
        _error("Gmsh tags not supported for dimension %i. Probably a bug" % dim)
        
    tags = tags_for_dim[highest_dim]
    physical_regions = tuple(tag[0] for tag in tags)
    if not all(tag == 0 for tag in physical_regions):
        handler.start_meshfunction("physical_region", dim, num_cells_counted)
        for i, physical_region in enumerate(physical_regions):
            handler.add_entity_meshfunction(i, physical_region)
        handler.end_meshfunction()
    
    # Now process the facet markers
    tags = tags_for_dim[highest_dim-1]
    if len(tags) > 0:
        
        print tags
        print vertices_used[highest_dim-1]
        
        physical_regions = tuple(tag[0] for tag in tags)
        if not all(tag == 0 for tag in physical_regions):
            mesh.init(highest_dim-1,0)
        
            # Get the facet-node connectivity information (reshape as a row of node indices per facet)
            facets_as_nodes = mesh.topology()(highest_dim-1,0)().reshape ( mesh.num_facets(), highest_dim )
            
#            from dolfin import MeshFunction
#            # Create and initialise the mesh function
#            facet_mark_function = MeshFunction ( 'uint', mesh, highest_dim-1 )
#            facet_mark_function.set_all( 0 )
            handler.start_meshfunction("facet_region", highest_dim-1, mesh.num_facets() )
            
            facets_to_check = range( mesh.num_facets() )
            
            data = [int(0*k) for k in range(len(facets_to_check)) ]
            
            for i, physical_region in enumerate(physical_regions):
                nodes = [n-1 for n in vertices_used[highest_dim-1][2*i:(2*i+highest_dim)]]
                nodes.sort()
                
                if physical_region != 0:
                    found = False
                    for j in range(len(facets_to_check)):
                        index = facets_to_check[j]
                        if all ( facets_as_nodes[index,k] == nodes[k] for k in range(len(nodes)) ):
                            found = True;
                            facets_to_check.pop(j)
                            # set the value of the mesh function
#                            facet_mark_function[index] = physical_region
                            data[index] = physical_region                
                            break;
                        
                    if not found:
                        raise Exception ( "The facet (%d) was not found to mark: %s" % (i, nodes) )
            
#            fname = os.path.splitext('tmp.xml')[0]
#            mesh_function_file = File("%s_%s.xml" % (fname, "facet_region"))
#            mesh_function_file << facet_mark_function
            
            for index, physical_region in enumerate ( data ):
                handler.add_entity_meshfunction(index, physical_region)
            handler.end_meshfunction()    
            
            mf = MeshFunction ( 'uint', mesh, 'tmp_facet_region.xml' )
            plot ( mf, interactive=True )
            
    # Check that we got all data
    if state == 10:
        print "Conversion done"
    else:
       _error("Missing data, unable to convert \n\ Did you use version 2.0 of the gmsh file format?")

    # Close files
    ifile.close()

        
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
    
    File ( ofilename ) << mesh
    
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


def to_mesh ( ifile ):
    me = meshconvert.XmlHandler( 'tmp.xml' )
    
    gmsh_to_dolfin_mesh( ifile, me )
    
    me.close()

if __name__ == "__main__":
    to_mesh( "untitled.msh" )
#    convert_and_create_facet_mesh_function ( "untitled.msh", "untitled.xml" )
    