!!!
!!! MODULE input_mesh :
!!! 
!!! Module that reads the femfeko native mesh format meshfile. The resulting 
!!! mesh info is stored in module geometry
!!!

MODULE input_mesh
  USE nrtype
  USE problem_info
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: input_mesh_read
CONTAINS

!!!
!!! SUBROUTINE input_mesh_read(meshfilename)
!!!
!!! The high level mesh reading driving routine that is intended to be called 
!!! by other code.
!!!
!!! input:
!!!
!!! meshfilename, which specifies the filename of the mesh to be read
!!!
!!! output:
!!!
!!! Calls module private routines read_nodes, ...
!!! These routines store the part of mesh input they handle in the appropriate
!!! structures in module geometry.
!!!

  SUBROUTINE input_mesh_read(meshfilename)
    USE unit_numbers
    
    IMPLICIT NONE
    CHARACTER(FILENAMELENGTH), INTENT(in) :: meshfilename
    
    OPEN(unit=meshin, status='OLD', file=meshfilename)
    CALL read_nodes(meshin)
    ! Blocks can be in any order in the meshfile
    REWIND(meshin)            
    CALL read_tets(meshin)

  END SUBROUTINE input_mesh_read
  

  SUBROUTINE read_nodes(meshfile)
    USE blockinput_tools
    USE geometry
    IMPLICIT NONE

!!!
!!! Interface variables
!!!
    INTEGER(i4b), INTENT(in) :: meshfile
!!!
!!! Internal variables
!!!
    INTEGER(i4b) :: ios, node_i, node_no
    REAL(dp) :: vert_x, vert_y, vert_z
    CHARACTER(optionlength) :: block_type
    CHARACTER(8) :: temp
    LOGICAL(lgt) :: nodes_found 


    nodes_found = .FALSE.       ! Nodes block not yet found 
    ! Find the nodes block in the mesh file
    find_block: DO
       block_type = blockinput_nextblock(meshfile, ios)
       IF (ios < 0) EXIT        ! End of file
       in_nodeblock: IF (block_type == 'nodes') THEN
          nodes_found = .TRUE.
          READ(meshfile, *) num_nodes ! Number of nodes
          ! Allocate memory for the nodes
          CALL allocate_geometry(alloc_vertices=num_nodes)
          node_loop: DO node_i = 1, num_nodes
             READ(meshfile, *) node_no, vert_x, vert_y, vert_z
             IF (node_no > num_nodes) THEN
                stop 'Node number in mesh file higher than num_nodes!'
             END IF
             ! Assign the node co-ordinates to variable in geometry module
             vertices(node_no)%coord = (/vert_x, vert_y, vert_z/)
          END DO node_loop

          READ(meshfile, '(A)', IOstat=ios) temp
          ! Check that we are indeed at the end of the block
          IF (temp.NE.'ENDBLOCK' .OR. ios < 0) THEN 
             STOP 'Incorrect number of nodes specified in BLOCK nodes'
          END IF

       END IF in_nodeblock
    END DO find_block
    IF (.NOT.nodes_found) THEN
       STOP 'No nodes block found in mesh file!'
    END IF

  END SUBROUTINE read_nodes

  SUBROUTINE read_tets(meshfile)
    USE blockinput_tools
    USE geometry
    IMPLICIT NONE

!!!
!!! Interface variables
!!!
    INTEGER(i4b), INTENT(in) :: meshfile
!!!
!!! Internal variables
!!!
    INTEGER(i4b) :: ios, tet_i, tet_no, tet_material, node_i
    INTEGER(i4b) :: node_a, node_b, node_c, node_d
    CHARACTER(optionlength) :: block_type
    CHARACTER(8) :: temp
    LOGICAL(lgt) :: tets_found 


    tets_found = .FALSE.       ! Nodes block not yet found 
    ! Find the nodes block in the mesh file
    find_block: DO
       block_type = blockinput_nextblock(meshfile, ios)
       IF (ios < 0) EXIT        ! End of file
       in_tetblock: IF (block_type == 'tets') THEN
          tets_found = .TRUE.
          READ(meshfile, *) num_elements ! Number of nodes
          ! Allocate memory for the nodes
          CALL allocate_geometry(alloc_elements=num_elements)
          tet_loop: DO tet_i = 1, num_elements
             READ(meshfile, *) tet_no, tet_material, node_a, node_b, node_c, node_d
             IF (tet_no > num_elements) THEN
                STOP 'Tetraheder number in mesh file higher than num_nodes!'
             END IF
             ! Assign the element nodes to variable in geometry module
             elements(tet_no)%nodes = (/node_a, node_b, node_c, node_d/)
             ! Check that non-existant node numbers weren't specified
             DO node_i = 1,4
                IF (elements(tet_no)%nodes(node_i) > num_nodes .OR.   &
                     elements(tet_no)%nodes(node_i) < 0) THEN
                   STOP 'Tetraheder specified with invalid node number'
                END IF
             END DO
             elements(tet_no)%material = tet_material ! material type
          END DO tet_loop
          
          READ(meshfile, '(A)', IOstat=ios) temp
          ! Check that we are indeed at the end of the block
          IF (temp.NE.'ENDBLOCK' .OR. ios < 0) THEN 
             STOP 'Incorrect number of tets specified in BLOCK tets'
          END IF

       END IF in_tetblock
    END DO find_block
    IF (.NOT.tets_found) THEN
       STOP 'No tets block found in mesh file!'
    END IF

    MAX_EDGES = num_elements*6   ! Maximum number of edges & faces.
    MAX_FACES = num_elements*4   ! Note: worst case - 
                                 ! 6 edges and 4 faces per tet.               
    
  END SUBROUTINE read_tets
    
END MODULE input_mesh


    
