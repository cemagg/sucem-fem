MODULE geomwrap
  USE nrtype
  USE geometry
  IMPLICIT NONE

  INTEGER(i4b), PARAMETER :: SPACEDIM=3
  INTEGER(i4b), PARAMETER :: ELNODES=4
  INTEGER(i4b), PARAMETER :: ELEDGES=6
  INTEGER(i4b), PARAMETER :: ELFACES=4
  INTEGER(i4b), PARAMETER :: FACENODES=3
  INTEGER(i4b), PARAMETER :: EDGENODES=2

  REAL(DP), DIMENSION(:,:), ALLOCATABLE ::  vertex_coords
  INTEGER(i4b), DIMENSION(:), ALLOCATABLE :: node_element_ptr
  INTEGER(i4b), DIMENSION(:), ALLOCATABLE :: node_elements
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: element_nodes
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: element_edges
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: element_faces
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: face_nodes
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: edge_nodes
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: element_connect2elem
  INTEGER(i4b), DIMENSION(:,:), ALLOCATABLE :: element_connect2face

CONTAINS

  SUBROUTINE init_geom()
    IMPLICIT NONE
    
    INTEGER(i4b) :: i, no_elements

    no_elements = SIZE(elements)

!!!
!!! Node elements connectivity pointer
!!!
    IF (ALLOCATED(node_element_ptr)) DEALLOCATE(node_element_ptr)
    ALLOCATE(node_element_ptr(SIZE(Node_ptr)))
    node_element_ptr = Node_ptr
!!!
!!! Node elements connectivity
!!!
    IF (ALLOCATED(node_elements)) DEALLOCATE(node_elements)
    ALLOCATE(node_elements(SIZE(Node_ind)))
    node_elements = Node_ind
!!!
!!! Vertex co-ordinates
!!!
    IF (ALLOCATED(vertex_coords)) DEALLOCATE(vertex_coords)
    ALLOCATE(vertex_coords(SPACEDIM,SIZE(vertices)))
    DO i=1,SIZE(vertices)
       vertex_coords(:,i) = vertices(i)%coord
    END DO
!!!
!!! Global element nodes
!!!
    IF (ALLOCATED(element_nodes)) DEALLOCATE(element_nodes)
    ALLOCATE(element_nodes(ELNODES,no_elements))
    DO i=1,no_elements
       element_nodes(:,i) = elements(i)%nodes
    END DO
!!!
!!! Global element edges
!!!
    IF (ALLOCATED(element_edges)) DEALLOCATE(element_edges)
    ALLOCATE(element_edges(ELEDGES,no_elements))
    DO i=1,no_elements
       element_edges(:,i) = elements(i)%edges
    END DO
!!!
!!! Global element faces
!!!
    IF (ALLOCATED(element_faces)) DEALLOCATE(element_faces)
    ALLOCATE(element_faces(ELFACES,no_elements))
    DO i=1,no_elements
       element_faces(:,i) = elements(i)%faces
    END DO
!!!
!!! element face 2 element connectivity
!!!
    IF (ALLOCATED(element_connect2elem)) DEALLOCATE(element_connect2elem)
    ALLOCATE(element_connect2elem(ELFACES, no_elements))
    DO i=1,no_elements
       element_connect2elem(:,i) = elements(i)%connect2elem
    END DO
!!!
!!! element face 2 face connectivity
!!!
    IF (ALLOCATED(element_connect2face)) DEALLOCATE(element_connect2face)
    ALLOCATE(element_connect2face(ELFACES, no_elements))
    DO i=1,no_elements
       element_connect2face(:,i) = elements(i)%connect2face
    END DO

!!!
!!! Global face nodes
!!!
    IF (ALLOCATED(face_nodes)) DEALLOCATE(face_nodes)
    ALLOCATE(face_nodes(FACENODES,SIZE(faces)))
    DO i=1,SIZE(faces)
       face_nodes(:,i) = faces(i)%nodes
    END DO
!!!
!!! Global edge nodes
!!!
    IF (ALLOCATED(edge_nodes)) DEALLOCATE(edge_nodes)
    ALLOCATE(edge_nodes(EDGENODES,SIZE(edges)))
    DO i=1,SIZE(edges)
       edge_nodes(:,i) = edges(i)%nodes
    END DO
!!! 

  END SUBROUTINE init_geom

END MODULE geomwrap
