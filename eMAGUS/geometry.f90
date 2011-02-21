! Last changed:
! 21 May 2003: Scat/tot field internal boundary search added.
! 14 Dec 2002: Additional iterative solver option
! added in FM card.
! Last changed 21 Feb 02. DBD:
!  Element, edge and face type augmented with mixed_order flag.

MODULE geometry 
  USE math_tools, ONLY: DET_DIM3, DET_DIM4, CROSS_PRODUCT, FIND_IN_LIST
  USE nrtype 
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  USE scattering_analysis_data, ONLY: sph_radius ! This variable was really used 
                                ! to hack around FEMAP's problems. This is a 
                                ! little messy, since logically this belongs 
                                ! with the scattering data. This variable will
                                ! probably go unused when more general scattering
                                ! geometry handling can be used. NM
  IMPLICIT NONE
!*******************************************************************************
! This module contains geometry-related data structures and routines.
! Created: 2001-09-15. MMB.
!*******************************************************************************
  SAVE

  ! Constants of the problem geometry:
  INTEGER(I4B) num_nodes     ! the number of nodes
  INTEGER(I4B) num_elements  ! the number of elements
  INTEGER(I4B) num_edges     ! the number of edges
  INTEGER(I4B) num_faces     ! the number of faces
  INTEGER(I4B) num_BCs       ! Number of BC's

  ! High level mesh information, should replace num_nodes, num_elements, etc.
  TYPE meshinfo_type
     INTEGER(I4B) num_nodes     ! the number of nodes
     INTEGER(I4B) num_elements  ! the number of elements
     INTEGER(I4B) num_edges     ! the number of edges
     INTEGER(I4B) num_PEC_edges ! Number of PEC edges
     INTEGER(I4B) num_faces     ! the number of faces
     INTEGER(I4B) num_PEC_faces ! Number of PEC faces
     INTEGER(I4B) num_BCs       ! Number of BC's
  END TYPE meshinfo_type

!!!
!!! TYPE element_info_type, currently used to pass information about the 
!!! element being integrated over in the S_(EDGE/FACE)_(EDGE/FACE) subroutines
!!! in module S_and_T_matrix.
!!!
  TYPE element_info_type
     INTEGER(i4b) :: num                   ! Global element number
     INTEGER(i4b) :: order                 ! Elemental order
     LOGICAL(LGT) :: mixed_order           ! Mixed order if true
  END TYPE element_info_type
!!!
    
  ! Vertex data:
  TYPE VERTEX
    REAL(DP), DIMENSION(3) :: coord                   !(x,y,z)
    INTEGER(I4B) BC_type                              ! To flag BC type.
                                                      ! 0 : Free
                                                      ! 1 : Port
                                                      ! 2 : PEC
                                                      ! 3 : PMC
                                                      ! 4 : Periodic
    LOGICAL(LGT) free                                 ! Free/not free flag.
  END TYPE VERTEX
  TYPE(VERTEX), DIMENSION(:), ALLOCATABLE, TARGET :: vertices

  ! Element data:
  TYPE ELEMENT
    INTEGER(I4B), DIMENSION(4) :: nodes               ! Global nodes
    INTEGER(I4B) order                                ! Hierarchal order
    LOGICAL(LGT) mixed                                ! Mixed or complete order
    INTEGER(I4B) material                             ! Material index
    INTEGER(I4B), DIMENSION(6) :: edges               ! Global edges
    INTEGER(I4B), DIMENSION(4) :: faces               ! Global faces
    INTEGER(I4B), DIMENSION(4) :: connect2elem        ! Edge element data structures
    INTEGER(I4B), DIMENSION(4) :: connect2face        ! See routine EDGEMAKE
    INTEGER(I4B), DIMENSION(4) :: hash                ! hash function for 
                                                      ! conectivity
    REAL(SP) :: residual                              ! volume residual integral
! Added DBD 16 July 2003
    REAL(SP) :: sigma_x                               ! PML parameters for the element.
    REAL(SP) :: sigma_y
    REAL(SP) :: sigma_z
    REAL(SP) :: Se_x
    REAL(SP) :: Se_y
    REAL(SP) :: Se_z
    REAL(SP) :: Te_x
    REAL(SP) :: Te_y
    REAL(SP) :: Te_z
! End added DBD 16 July 2003
!!! Added DBD 13 July 2005
!!!
!!! Global mid-point nodes. This is an experimental data-structure
!!! for curvilinear test elements, to avoid disrupting other parts of the code.
!!! Mid-point nodes are numbered as in Fig 5.3, p. 177, J. Jin	  
!!! "The FEM in EM", 2nd edn, Wiley 2002, but with nodes 9 and 10 interchanged: i.e. 
!!! node  5 on edge 1 between 1 and 2;
!!! node  6 on edge 2 between 1 and 3; 
!!! node  7 on edge 3 between 1 and 4;
!!! node  8 on edge 4 between 2 and 3; 
!!! node  9 on edge 5 between 2 and 4;
!!! node 10 on edge 6 between 3 and 4.
 
    INTEGER(I4B), DIMENSION(6) :: mid_point_nodes     
!!! End added DBD 13 July 2005


  END TYPE ELEMENT
  TYPE(ELEMENT), DIMENSION(:), ALLOCATABLE, TARGET :: elements

  ! Edge data:
  TYPE EDGE
    INTEGER(I4B), DIMENSION(2) :: nodes               ! Global nodes
    INTEGER(I4B) order                                ! Hierarchal order
    LOGICAL(LGT) mixed                                ! Mixed or complete order
    LOGICAL(LGT) coax_aperture
    INTEGER(I4B) coaxnumber
    LOGICAL(LGT) free
    LOGICAL(LGT) scat_tot_boundary
    LOGICAL(LGT) PEC
    LOGICAL(LGT) CBAA_aperture
    LOGICAL(LGT) Dirichlet
    LOGICAL(LGT) port
    INTEGER(I4B) portnumber
! DBD addition 25 March 2003
    LOGICAL(LGT) ABC
    INTEGER(I4B) ABCnumber
!    REAL(SP) ABC_Yc
! End DBD addition 25 March 2003
  END TYPE EDGE
  TYPE(EDGE), DIMENSION(:), ALLOCATABLE, TARGET :: edges

  ! Face data:
  TYPE FACE
    INTEGER(I4B), DIMENSION(3) :: nodes               ! Global nodes
    INTEGER(I4B) order                                ! Hierarchal order
    LOGICAL(LGT) mixed                                ! Mixed or complete order
    LOGICAL(LGT) coax_aperture
    INTEGER(I4B) coaxnumber
    LOGICAL(LGT) free
    LOGICAL(LGT) scat_tot_boundary
    LOGICAL(LGT) curvilinear
    LOGICAL(LGT) PEC
    LOGICAL(LGT) CBAA_aperture
    LOGICAL(LGT) Dirichlet
    LOGICAL(LGT) port
    INTEGER(I4B) portnumber
    REAL(SP) :: residual                              ! surface residual integral
! DBD addition 25 March 2003
    LOGICAL(LGT) ABC
    INTEGER(I4B) ABCNumber
!    REAL(SP) ABC_Yc
! End DBD addition 25 March 2003
  END TYPE FACE
  TYPE(FACE), DIMENSION(:), ALLOCATABLE, TARGET :: faces

  ! Boundary conditions data:
  TYPE BC
    REAL(SP), DIMENSION(4,3) :: corner                !(x,y,z) coords of each corner.
                                                      ! As for type VERTEX.
    INTEGER(I4B) type                                 ! To flag BC type.
    INTEGER(I4B) port_num                             ! Unique identifiying label assigned to this port (if type=1)
! DBD addition 25 March 2003
    INTEGER(I4B) ABC_num                              ! Unique identifiying label assigned to this ABC (if type=5)
! End DBD addition 25 March 2003
  END TYPE BC
  TYPE(BC), DIMENSION(:), ALLOCATABLE :: BCs

  ! PML data:
  TYPE PerfectlyMatchedLayer
    INTEGER(I4B) m
!    INTEGER(I4B) material_label
!    INTEGER(I4B) material_label_offset
	REAL(SP) sigma_max
	REAL(SP) thickness
	LOGICAL(LGT) absorb_x
	LOGICAL(LGT) absorb_y
	LOGICAL(LGT) absorb_z	
	LOGICAL(LGT) full
	REAL(SP) x_min
	REAL(SP) y_min
	REAL(SP) z_min
	REAL(SP) x_plus
	REAL(SP) y_plus
	REAL(SP) z_plus
  END TYPE PerfectlyMatchedLayer
  TYPE(PerfectlyMatchedLayer) :: PML
  LOGICAL(LGT) PML_PRESENT 

 


  INTEGER(I4B), PARAMETER :: HOMOG_MEDIUM = 1       ! Relevant only for TD scattering analysis using scattered/total field 
                                                    ! approach. 

  
  REAL(SP) avg_edge_length                          ! Average mesh edge length

  ! Connectivity data (previously MODULE connect_geometry, moved MMB 18 Sept 2001):
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: edgeconnectelem, edgeconnectedge, &
                                               faceconnectface, faceconnectelem

  ! These are pointers and array indices that are used 
  ! for fast searching through elements / nodes 
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: Node_ptr, Node_ind 

!*******************************************************************************
! Internal routines:
! (A routine is either listed as an interface, or in the private statement.)

! DBD addition 25 March 2003

!PRIVATE
!*******************************************************************************

CONTAINS
!*******************************************************************************

! DBD addition 25 March 2003
SUBROUTINE ABC_BOUNDARY_SEARCH
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE   
!*******************************************************************************
! Cycles through all EXTERIOR faces (unconnected faces) and flags them an ABC
! if appropriate. Routine works for time-domain open region problems
! terminated with explicit ABC cards (AB).
! Written DBD 25 March 2003.
! Extended DBD 18 Dec 2003 to include frequency domain open region problems
! terminated in an open and otherwise unspecified boundary.
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem, ABC_counter
  INTEGER(I4B), DIMENSION(3) :: tempedges 
  LOGICAL(LGT) :: in_plane
  LOGICAL(LGT), ALLOCATABLE :: ABC_found(:)     ! A flag to check that an ABC
                                                ! does indeed lie (at least partially) on the mesh. 

  IF (TD_ANALYSIS) THEN
    ALLOCATE(ABC_found(num_ABCs))
  ELSE                   
    ALLOCATE(ABC_found(1)) ! Only one ABC otherwise.
  END IF
  ABC_found = .FALSE. ! Array allocation


  IF (TD_ANALYSIS) THEN
    ELEMENT_LOOP1: DO ielem = 1,num_elements
      FACE_LOOP1: DO iface = 1,4 ! Face-wise search to find edges in the quad.
        FACE_CONNECT_TEST1: IF (elements(ielem)%connect2elem(iface).EQ.0) THEN     
          tempedges = LOCAL_FACEEDGES(iface)
          ABC_LOOP: DO ABC_counter = 1,num_ABCs
            CALL FACE_IN_QUADRILATERAL( & 
               BCs(ABCs(ABC_counter)%BC_num)%corner(1,1:3), &
               BCs(ABCs(ABC_counter)%BC_num)%corner(2,1:3), &
               BCs(ABCs(ABC_counter)%BC_num)%corner(3,1:3), &
               BCs(ABCs(ABC_counter)%BC_num)%corner(4,1:3), &
               elements(ielem)%faces(iface),in_plane)
            IF (in_plane) THEN
              ! Flag the three edges and the face as an ABC. 
	          edges(elements(ielem)%edges(tempedges))%ABC = .true.
!              edges(elements(ielem)%edges(tempedges))%ABC_Yc = ABCs(ABC_counter)%Yc
              edges(elements(ielem)%edges(tempedges))%ABCNumber = ABCs(ABC_counter)%label
              faces(elements(ielem)%faces(iface))%ABC = .true.
 !             faces(elements(ielem)%faces(iface))%ABC_Yc = ABCs(ABC_counter)%Yc
              faces(elements(ielem)%faces(iface))%ABCNumber = ABCs(ABC_counter)%label
	          ABC_found(ABC_counter) = .TRUE. ! Will usually be overwritten by subsequent elements.
            END IF
          END DO ABC_LOOP
        END IF FACE_CONNECT_TEST1
      END DO FACE_LOOP1
    END DO ELEMENT_LOOP1
  ELSE IF(FD_SCAT_ANALYSIS) THEN
    ELEMENT_LOOP2: DO ielem = 1,num_elements
      FACE_LOOP2: DO iface = 1,4 ! Face-wise search to find outer boundary faces and edges.
	    FACE_CONNECT_TEST2: IF (elements(ielem)%connect2elem(iface).EQ.0) THEN     
          tempedges = LOCAL_FACEEDGES(iface)
          ! Flag the three edges and the face as an ABC. 
          edges(elements(ielem)%edges(tempedges))%ABC = .true.
          faces(elements(ielem)%faces(iface))%ABC = .true.
	      ABC_found(1) = .TRUE. 
        END IF FACE_CONNECT_TEST2
      END DO FACE_LOOP2
    END DO ELEMENT_LOOP2
  ELSE
    STOP 'ABC not supported for other types of analysis at present'
  END IF
   
  ! Check that at least one element (face) was found per ABC.
  DO ABC_counter = 1,num_ABCs
    IF(.NOT.ABC_found(ABC_counter)) THEN
      CALL ERROR_FEMFEKO(1,4700,ABC_counter)
    END IF
  END DO
  DEALLOCATE(ABC_found)
    
END SUBROUTINE ABC_BOUNDARY_SEARCH
! End DBD addition 25 March 2003
!*******************************************************************************




SUBROUTINE ALLOCATE_GEOMETRY(alloc_vertices,alloc_edges,alloc_faces,alloc_elements)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine handles allocation of the arrays in module geometry.
! Can allocate space for nodes, edges, faces, elements.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN), OPTIONAL :: alloc_vertices
  INTEGER(I4B), INTENT(IN), OPTIONAL :: alloc_edges
  INTEGER(I4B), INTENT(IN), OPTIONAL :: alloc_faces
  INTEGER(I4B), INTENT(IN), OPTIONAL :: alloc_elements

  ! Allocate storage for vertex data:
  IF (PRESENT(alloc_vertices)) THEN
    IF (alloc_vertices.GT.0) THEN
      IF (ALLOCATED(vertices)) DEALLOCATE(vertices)
      ALLOCATE(vertices(alloc_vertices))
    ELSE
      STOP 'IE: Negative number of vertices in ALLOCATE_GEOMETRY.'
    END IF
  END IF

  ! Allocate storage for edge data:
  IF (PRESENT(alloc_edges)) THEN
    IF (alloc_edges.GT.0) THEN
      IF (ALLOCATED(edges)) DEALLOCATE(edges)
      ALLOCATE(edges(alloc_edges))
    ELSE
      STOP 'IE: Negative number of edges in ALLOCATE_GEOMETRY.'
    END IF
  END IF

  ! Allocate storage for face data:
  IF (PRESENT(alloc_faces)) THEN
    IF (alloc_faces.GT.0) THEN
      IF (ALLOCATED(faces)) DEALLOCATE(faces)
      ALLOCATE(faces(alloc_faces))
    ELSE
      STOP 'IE: Negative number of faces in ALLOCATE_GEOMETRY.'
    END IF
  END IF

  ! Allocate storage for element data:
  IF (PRESENT(alloc_elements)) THEN
    IF (alloc_elements.GT.0) THEN
      IF (ALLOCATED(elements)) DEALLOCATE(elements)
      ALLOCATE(elements(alloc_elements))
    ELSE
      STOP 'IE: Negative number of elements in ALLOCATE_GEOMETRY.'
    END IF
  END IF

END SUBROUTINE ALLOCATE_GEOMETRY
!*******************************************************************************


SUBROUTINE CBAA_EDGESORT
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine sorts the edges so that all edges lying in the CBAA aperture
! (z=0 plane) are numbered first, ythis is necessary to insure that the CBAA
! aperture dofs are numbered first.
! MMB 2001-10-04
!*******************************************************************************
  INTEGER(I4B) :: ielem,iedge,jedge,ap_edcount
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: ap_edges

  ! Cycle through all edges and record those that lie on the CBAA aperture:
  ALLOCATE(ap_edges(num_edges))
  ap_edges = .FALSE. ! array assignment
  ap_edcount = 0
  DO iedge = 1,num_edges
    IF ((ABS(vertices(edges(iedge)%nodes(1))%coord(3)).LT.EPS).AND. &
        (ABS(vertices(edges(iedge)%nodes(2))%coord(3)).LT.EPS)) THEN ! edge in aperture
      ap_edcount = ap_edcount + 1
      ap_edges(iedge) = .TRUE.
    END IF
  END DO

  ! Now, starting with edge 1, swop every non-aperture edge with 
  ! the next aperture edge found:
  EDGE_LOOP: DO iedge = 1,ap_edcount
    AP_CHECK1: IF (.NOT.ap_edges(iedge)) THEN
      FIND_LOOP: DO jedge = iedge+1,num_edges
        AP_CHECK2: IF (ap_edges(jedge)) THEN

          CALL EDGE_SWOP(iedge,jedge)
          ap_edges(jedge) = .FALSE.
          EXIT FIND_LOOP

        END IF AP_CHECK2
      END DO FIND_LOOP
    END IF AP_CHECK1
  END DO EDGE_LOOP

  ! Deallocate temporary storage:
  DEALLOCATE(ap_edges)

CONTAINS

  SUBROUTINE EDGE_SWOP(edge1,edge2)
    IMPLICIT NONE
!*******************************************************************************
! This routine swops two edges. Affected data structures:
! * edges
! * elements
!*******************************************************************************
INTEGER(I4B), INTENT(IN) :: edge1,edge2

    INTEGER(I4B), DIMENSION(4) :: edgenodes
    TYPE(EDGE) :: tempedge
    INTEGER(I4B) :: listlen,inode,pos1,pos2,ielem,elemval,jelem,thisedge
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ellist

    ! Swop the edge data:
    tempedge%nodes     = edges(edge1)%nodes
    edges(edge1)%nodes = edges(edge2)%nodes
    edges(edge2)%nodes = tempedge%nodes

    ! Construct a list of elements to be searched:
    edgenodes(1:2) = edges(edge1)%nodes
    edgenodes(3:4) = edges(edge2)%nodes
    listlen = 0
    DO inode = 1,4
      listlen = listlen - Node_ptr(edgenodes(inode)) + Node_ptr(edgenodes(inode)+1)
    END DO
    ALLOCATE(ellist(listlen))
    pos1 = 1
    DO inode = 1,4
      pos2 = (pos1-1) - Node_ptr(edgenodes(inode)) + Node_ptr(edgenodes(inode)+1)
      ellist(pos1:pos2) = &
        Node_ind(Node_ptr(edgenodes(inode))+1:Node_ptr(edgenodes(inode)+1))
      pos1 = pos2 + 1
    END DO

    ! Replace double occurances of element numbers in <ellist>, with -1:
    DO ielem = 1,listlen
      elemval = ellist(ielem)
      DO jelem = ielem+1,listlen
        IF (ellist(jelem).EQ.elemval) ellist(jelem) = -1
      END DO
    END DO

    ! Swop the edge numbers in the elements listed in <ellist>:
    SWOP_LOOP: DO ielem = 1,listlen
      elemval = ellist(ielem)
      IF (elemval.EQ.-1) CYCLE SWOP_LOOP
      EDGE_RENUM: DO thisedge = 1,6
        IF (elements(elemval)%edges(thisedge).EQ.edge1) THEN
          elements(elemval)%edges(thisedge) = edge2
          CYCLE EDGE_RENUM
        END IF
        IF (elements(elemval)%edges(thisedge).EQ.edge2) THEN
          elements(elemval)%edges(thisedge) = edge1
          CYCLE EDGE_RENUM
        END IF
      END DO EDGE_RENUM
    END DO SWOP_LOOP

    DEALLOCATE(ellist)

  END SUBROUTINE EDGE_SWOP

END SUBROUTINE CBAA_EDGESORT
!*******************************************************************************


! DBD addition 26 July 2005
SUBROUTINE CURVILINEAR_BOUNDARY_SEARCH_SPH
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! This routine finds faces lying on the boundary between the internal 
! inhomegeneous SPHERICAL materials (material labels not equal to HOMOG_MEDIUM) 
! and external homogeneous material label HOMOG_MEDIUM and flags them appropriately 
! for treatment by curvilinear elements. 
! IT ONLY WORKS FOR SPHERES, CENTERED ON THE ORIGIN, AND THIS IS NOT CHECKED!!
!
! Written by DBD, 26 July 2005, based on SCAT_TOTAL_BOUNDARY_SEARCH_SPH.
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem,inodes
 ! INTEGER(I4B), DIMENSION(3) :: tempedges ! Add in later if edges also need to be flagged, 
                                           ! as well as two other commented out lines below.
  INTEGER(I4B), DIMENSION(3) :: tempnodes
  REAL(SP), DIMENSION(3) :: rad ! radius of nodes. 
  REAL(SP), DIMENSION(3) :: diff   
  LOGICAL(LGT) :: boundary_found ! A flag to ensure that at least one element has been found.

  ! Initialize
  faces%curvilinear = .FALSE.
  boundary_found    = .FALSE.

  ELEMENT_LOOP: DO ielem = 1,num_elements
    IF (elements(ielem)%material.NE.HOMOG_MEDIUM) THEN
	  ! This is in the inhomogeneous target region. Test if it lies on the (fictitious) boundary between
	  ! the scattered/total field region, by seeing if all the nodes on a face lie on the same radius. 
      DO iface = 1,4
!        tempedges = LOCAL_FACEEDGES(iface)
        tempnodes = GLOBAL_FACENODES(ielem,iface)
		DO inodes = 1,3
          rad(inodes) = SQRT(DOT_PRODUCT(vertices(tempnodes(inodes))%coord,vertices(tempnodes(inodes))%coord))
		  diff(inodes) = ABS(rad(inodes)-sph_radius)
        END DO
        IF (SUM(diff).LT.0.001_SP*sph_radius) THEN
! write(fileout,*) 'Element ',ielem,'  Face  ',iface,'  flagged as curvilinear'
          faces(elements(ielem)%faces(iface))%curvilinear = .TRUE.
          WRITE(*,*) "curvilinear face found in sub CURVILINEAR_BOUNDARY_SEARCH_SPH"
          boundary_found    = .TRUE.
		END IF
	  END DO
    END IF    
  END DO ELEMENT_LOOP
  IF(.NOT.boundary_found) THEN
    WRITE(FILEOUT,*) 'IE: ERROR IN CURVILINEAR_BOUNDARY_SEARCH_SPH, NO ELEMENTS FLAGGED AS CURVILINEAR. CHECK SPHERICAL RADIUS.' 
  END IF

END SUBROUTINE CURVILINEAR_BOUNDARY_SEARCH_SPH

! End DBD addition 30 July 2004
!*******************************************************************************





SUBROUTINE DEALLOCATE_GEOMETRY(dealloc_vertices,dealloc_edges,dealloc_faces, &
                               dealloc_elements,dealloc_connect)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine deallocates the geometry data defined inside module geometry.
!*******************************************************************************
  LOGICAL(LGT), INTENT(IN), OPTIONAL :: dealloc_vertices
  LOGICAL(LGT), INTENT(IN), OPTIONAL :: dealloc_edges
  LOGICAL(LGT), INTENT(IN), OPTIONAL :: dealloc_faces
  LOGICAL(LGT), INTENT(IN), OPTIONAL :: dealloc_elements
  LOGICAL(LGT), INTENT(IN), OPTIONAL :: dealloc_connect

  ! Deallocate vertex data:
  IF (PRESENT(dealloc_vertices)) THEN
    IF (ALLOCATED(vertices).AND.dealloc_vertices) DEALLOCATE(vertices)
  END IF

  ! Deallocate edge data:
  IF (PRESENT(dealloc_edges)) THEN
    IF (ALLOCATED(edges).AND.dealloc_edges) DEALLOCATE(edges)
  END IF

  ! Deallocate face data:
  IF (PRESENT(dealloc_faces)) THEN
    IF (ALLOCATED(faces).AND.dealloc_faces) DEALLOCATE(faces)
  END IF

  ! Deallocate element data:
  IF (PRESENT(dealloc_elements)) THEN
    IF (ALLOCATED(elements).AND.dealloc_elements) DEALLOCATE(elements)
  END IF

  ! Deallocate connectivity data:
  IF (PRESENT(dealloc_connect)) THEN
    IF (dealloc_connect) THEN
      IF (ALLOCATED(edgeconnectelem)) DEALLOCATE(edgeconnectelem)
      IF (ALLOCATED(edgeconnectedge)) DEALLOCATE(edgeconnectedge)
      IF (ALLOCATED(faceconnectface)) DEALLOCATE(faceconnectface)
      IF (ALLOCATED(faceconnectelem)) DEALLOCATE(faceconnectelem)
    END IF
  END IF

END SUBROUTINE DEALLOCATE_GEOMETRY
!*******************************************************************************


! Added DBD 13 July 05 
FUNCTION EDGE_CENTRE(ielem,jedge)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Return coordinates of centre of element ielem's edge jedge..
! DBD 13 July 2005.
!*******************************************************************************
  REAL(SP), DIMENSION(3) :: EDGE_CENTRE
  INTEGER(I4B), INTENT(IN)  :: ielem,jedge
  INTEGER(I4B), DIMENSION(2) :: edgenodes
  edgenodes = GLOBAL_EDGENODES(ielem,jedge)
  EDGE_CENTRE(1) = SUM(vertices(edgenodes(1:2))%coord(1))/2.0_SP 
  EDGE_CENTRE(2) = SUM(vertices(edgenodes(1:2))%coord(2))/2.0_SP 
  EDGE_CENTRE(3) = SUM(vertices(edgenodes(1:2))%coord(3))/2.0_SP 
 
END FUNCTION EDGE_CENTRE
!*******************************************************************************
! End added DBD 13 July 05




SUBROUTINE EDGE_FACE_OUTPUT
  USE boundary_conditions
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE   

!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine writes out edge and face information, if requested.
! Note that although it might logically appear to belong within the output module,
! this creates a circular reference since it also needs access to geometry.
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! DB Davidson.
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 2 April 2003 DBD.
!
!*******************************************************************************
   INTEGER(I4B) i_edge,i_elem,i_face
   IF (OUTPUT_ELEMENT_DATA) THEN
     WRITE(FILEOUT,'(//,20X,A)') 'LIST OF EDGES'
     WRITE(FILEOUT,'(3X,3(A12),5(A7))') 'Edge no.', 'Node 1', 'Node 2', & 
	                              'PEC', 'Port', 'Port #', 'ABC', 'ABC #'
     DO i_edge = 1,NUM_edges
       WRITE(FILEOUT,'(1X,3(4X,I8),6X,L2,2(6X,L2,2X,I4,1X))') i_edge, edges(i_edge)%nodes(1:2), &
            edges(i_edge)%PEC, &
	        edges(i_edge)%port,edges(i_edge)%portnumber,&
			edges(i_edge)%ABC,edges(i_edge)%ABCnumber
     END DO

     WRITE(FILEOUT,'(//,20X,A)') 'LIST OF FACES'
     WRITE(FILEOUT,'(3X,4(A12),5(A7))') 'Face no.', 'Node 1', 'Node 2', 'Node 3',&
	                              'PEC', 'Port', 'Port #', 'ABC', 'ABC #'
     DO i_face = 1,NUM_faces
       WRITE(FILEOUT,'(1X,4(4X,I8),6X,L2,2(6X,L2,2X,I4,1X))') i_face, faces(i_face)%nodes(1:3), & 
	        faces(i_face)%PEC, &
	        faces(i_face)%port,faces(i_face)%portnumber,&
			faces(i_face)%ABC,faces(i_face)%ABCnumber
     END DO

     WRITE(FILEOUT,'(//,20X,A)') 'ELEMENT FACE INTERCONNECTIVITY DATA'
     WRITE(FILEOUT,'(3X,4(A12))') 'Element', 'Face', 'Connected to elements', 'Connected to faces of these elements' 
     DO i_elem =1,num_elements
	   DO i_face = 1,4
	     WRITE(FILEOUT,'(1X,2(4X,I8),6X,2(1X,4(1X,I8)))') i_elem, i_face,elements(i_elem)%connect2elem(1:4),&
		 elements(i_elem)%connect2face(1:4)
       END DO
     END DO


   END IF
END SUBROUTINE EDGE_FACE_OUTPUT

!*******************************************************************************


FUNCTION EDGE_VECTOR(first_node,second_node,i)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! The edge vector from the first to second node of element i.
!*******************************************************************************
   REAL(SP), DIMENSION(3) :: EDGE_VECTOR
   INTEGER(I4B), INTENT(IN) :: first_node,second_node,i

   EDGE_VECTOR = vertices(elements(i)%nodes(second_node))%coord - & 
                 vertices(elements(i)%nodes(first_node))%coord 

END FUNCTION EDGE_VECTOR

! Added DBD 21 May 2003
!*******************************************************************************
FUNCTION EDGE_UNIT_VECTOR(first_node,second_node,i)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! The edge vector of unit length from the first to second node of element i.
!*******************************************************************************
   REAL(SP), DIMENSION(3) :: EDGE_UNIT_VECTOR
   INTEGER(I4B), INTENT(IN) :: first_node,second_node,i

   EDGE_UNIT_VECTOR = vertices(elements(i)%nodes(second_node))%coord - & 
                      vertices(elements(i)%nodes(first_node))%coord 
   EDGE_UNIT_VECTOR = EDGE_UNIT_VECTOR/SQRT(DOT_PRODUCT(EDGE_UNIT_VECTOR,EDGE_UNIT_VECTOR))

END FUNCTION EDGE_UNIT_VECTOR
!*******************************************************************************
! End added DBD 21 May 2003

! Added DBD 06 June 2003
FUNCTION ELEMENT_CENTRE(ii)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Return coordinates of centre of element ii
! DBD 06 June 2003
!*******************************************************************************
  REAL(SP), DIMENSION(3) :: ELEMENT_CENTRE
  INTEGER(I4B), INTENT(IN)  :: ii
  
  ELEMENT_CENTRE(1) = SUM(vertices(elements(ii)%nodes(1:4))%coord(1))/4.0_SP 
  ELEMENT_CENTRE(2) = SUM(vertices(elements(ii)%nodes(1:4))%coord(2))/4.0_SP 
  ELEMENT_CENTRE(3) = SUM(vertices(elements(ii)%nodes(1:4))%coord(3))/4.0_SP

 
END FUNCTION ELEMENT_CENTRE
!*******************************************************************************
! Added DBD 06 June 2003

FUNCTION ELEMENT_NORMAL_VECTOR(elem,local_face)
  USE math_tools, ONLY: FIND_IN_LIST,VECTOR_LENGTH
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This function returns the outward, unit normal vector to local face 
!<local_face> of element <elem>. It uses the fact that the simplex coord
! gradient vectors are inward-normal to the opposing faces to the associated
! nodes.
!
! 2002-03-09: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face
  REAL(SP), DIMENSION(3) :: ELEMENT_NORMAL_VECTOR

  REAL(SP), DIMENSION(4,3) :: grad_lam
  INTEGER(I4B)  :: gfortran_bug_tmp(3)
  INTEGER(I4B), DIMENSION(4) :: node_flags
  INTEGER(I4B) :: face_node
  
  ! Find the simplex coordinate gradients. They are normal to the faces:
  grad_lam = GRADIENT_LAMBDA(elem,.TRUE.)

  ! Find the node opposite the face <local_face>:
  node_flags = 0
  gfortran_bug_tmp = LOCAL_FACENODES(local_face)
  node_flags(gfortran_bug_tmp) = 1
  CALL FIND_IN_LIST(node_flags,4,0,face_node)
  IF (face_node.EQ.0) STOP 'IE: invalid node in ELEMENT_NORMAL_VECTOR.'

  ! The gradient vector must be reversed to point outward and must also be 
  ! normalised, yielding the final result:
  ELEMENT_NORMAL_VECTOR = -(1.0/VECTOR_LENGTH(grad_lam(face_node,1:3))) * &
                          grad_lam(face_node,1:3)

END FUNCTION ELEMENT_NORMAL_VECTOR
!*******************************************************************************


FUNCTION FACE_AREA(element_num,face_num)
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!******************************************************************************
! A function to compute the actual (unsigned) area of face_num of element_num.
! This uses the identity A_i = \nabla \lambda_i *3V [Lee and Mittra,eqn. 15]
! where A_i is the area normal to triangular face comprising nodes {j,k,l} 
! of element {i,j,k,l}. 
! Note that the S&P numbering convention is slightly different, 
! (see for eg. documentation for GLOBAL_FACENODES) and is taken into account
! here.
! DBD 1997. Revised Jan-Feb 1999.
!******************************************************************************
  INTEGER(I4B), INTENT(IN)  :: element_num,face_num
  REAL(SP), DIMENSION(4,3) :: grad_lam
  REAL(SP), DIMENSION(3) :: VECT_AREA
  REAL(SP) FACE_AREA
  INTEGER(I4B) :: simplex_coord

  SELECT CASE (face_num)
  CASE (1)
    simplex_coord = 4
  CASE(2)
    simplex_coord = 3
  CASE(3)
    simplex_coord = 2
  CASE(4)
    simplex_coord = 1
  CASE DEFAULT    
    STOP 'IE: Call with invalid face number in function FACE_AREA.'
  END SELECT
 
  grad_lam = GRADIENT_LAMBDA(element_num,.FALSE.)
  VECT_AREA = 0.5_SP*grad_lam(simplex_coord,1:3)
  ! Rather than divide by 6V and then multiply by 3V, only the factor of 1/2
  ! is included. 
  FACE_AREA = SQRT(DOT_PRODUCT(VECT_AREA,VECT_AREA))
  ! This is by definition a positively signed quantity.
  IF(DEBUG_AREA) THEN ! Output debugging info
     WRITE(FILEOUT,'(//A,A)') 'Element area'
     WRITE(FILEOUT,'(A,I8,1X,A,I1,1X,A,1X,G10.4)') & 
                                  'Element number: ',element_num, & 
                                  'Face number: ',face_num, &
                                  'Area: ', FACE_AREA
  END IF
END FUNCTION FACE_AREA
!*******************************************************************************


SUBROUTINE FACE_AREA_VECTORS(element_num,area)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This subroutine returns the area vectors for the four faces of the 
! tetrahedron.
! See Lee and Mittra, "Edge-elements for 3-D inhomegoneously-filled cavities", 
! IEEE MTT, Sep 192. 1767-1773 for theory. This routine 
! essentially implements egn (15) thereof.
! The gradient of the barycentric funtons (\nabla /lambda_{i} in eqn 15)
! are computed used eqn. (5.9), p.298, from Silvester and Ferrari, 
! "Finite Elements for Electrical Engineers", 3dr edition, 
! Cambridge Univ Press 1996.
! Note the S&F numbers the nodes from 1 to 4; Lee and Mittra from 0 to 3.
! This routine uses the former convention.
!
! Output data is the matrix area. Note that there is a 3V in eqn.(15) [L&M]
! which (almost) cancels the 6V in eqn.(5.9) [S&F,p298], leaving a 1/2.
! Note also that the latter eqn. give the OUTWARD normals;  
! [L&M] use INWARD normals. This is corrected in the code by
! multiplying by -0.5
! Since the area is always used in either a dot- or cross-product
! with another area term of the same element, the direction of the normal 
! is  actually irrelevant.
!
! 18 March 1997 by DBD. 
! 27 May 1997. Most functionality replaced by calls to function 
! GRADIENT_LAMBDA.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: element_num
  REAL(SP), DIMENSION(4,3), INTENT(OUT) :: area

  area = GRADIENT_LAMBDA(element_num,.false.)
  area = -0.5_SP * area

END SUBROUTINE FACE_AREA_VECTORS
!*******************************************************************************


! Added DBD 27 Jan 05
FUNCTION FACE_CENTRE(ielem,jface)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Return coordinates of centre of element ielem's face jface.
! DBD 06 June 2003
!*******************************************************************************
  REAL(SP), DIMENSION(3) :: FACE_CENTRE
  INTEGER(I4B), INTENT(IN)  :: ielem,jface
  INTEGER(I4B), DIMENSION(3) :: facenodes
  facenodes = GLOBAL_FACENODES(ielem,jface)

  FACE_CENTRE(1) = SUM(vertices(facenodes(1:3))%coord(1))/3.0_SP 
  FACE_CENTRE(2) = SUM(vertices(facenodes(1:3))%coord(2))/3.0_SP 
  FACE_CENTRE(3) = SUM(vertices(facenodes(1:3))%coord(3))/3.0_SP 
 

END FUNCTION FACE_CENTRE
!*******************************************************************************
! End added DBD 27 Jan 05


SUBROUTINE FACE_IN_CIRCLE(normal,center,radius,face_val,in_circ)
  USE math_tools, ONLY: VECTOR_LENGTH
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! This subroutine recieves the circle's plane-normal, center coordinates and
! radius. The flag <in_circ> is set to indicate whether the global face <face_val>
! is located within the circle. 
! 2002-02-06: Created. MMB
!*******************************************************************************
  REAL(SP), DIMENSION(3), INTENT(IN) :: normal,center
  REAL(SP), INTENT(IN) :: radius
  INTEGER(I4B), INTENT(IN) :: face_val
  LOGICAL(LGT), INTENT(OUT) :: in_circ

  REAL(SP), DIMENSION(4) :: pl_eq          ! equation for plane: ax+by+cz+d=0
  LOGICAL(LGT), DIMENSION(3) :: node_flags ! flags for all face nodes' aperture status
  REAL(SP), DIMENSION(3) :: node_coord     ! coordinates of a face node
  INTEGER(I4B) :: inode                    ! counters

  ! Set up equation to describe the circle plane:
  pl_eq(1:3) = normal             ! 0 = ax + by + cz
  pl_eq(4) = - SUM(normal*center) ! + d (pl_eq = [a b c d])

  ! Assign the flag values:
  node_flags = .FALSE.   ! intialise

  DO inode = 1,3 ! cycle through the 3 face nodes 
    node_coord = vertices(faces(face_val)%nodes(inode))%coord(1:3)
	
	! Check whether node is in the circle plane:
    IF (ABS(SUM(pl_eq(1:3)*node_coord)+pl_eq(4)).LT.EPS) THEN

      ! Check whether node is within the radius distance of the center:
      IF ((VECTOR_LENGTH(node_coord-center)/radius-1.0).LT.EPS) THEN
        node_flags(inode) = .TRUE.
      END IF
    END IF
  END DO
  in_circ = node_flags(1).AND.node_flags(2).AND.node_flags(3)

END SUBROUTINE FACE_IN_CIRCLE
!*******************************************************************************


SUBROUTINE FACE_IN_QUADRILATERAL(r1,r2,r3,r4,face_val,in_quad)
  USE math_tools, ONLY: CROSS_PRODUCT
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE   
!*******************************************************************************
! This subroutine recieves the quadrilateral's nodes (r1-r4) and then sets the 
! flag <in_quad> to indicate whether the global face <face_val> is located 
! within the quadrilateral. 
! Note: All interior quad. angles must be < 180 deg.
!*******************************************************************************
  REAL(SP), DIMENSION(3), INTENT(IN) :: r1,r2,r3,r4     ! Quad. nodes in right hand order
  INTEGER(I4B), INTENT(IN) :: face_val                  ! Global face number
  LOGICAL(LGT), INTENT(OUT) :: in_quad                  ! results
  LOGICAL(LGT), DIMENSION(3) :: node_flags              ! calculated (if needed in future)
  REAL(SP), DIMENSION(3) :: pl_norm,basis1,basis2       ! basis vectors for plane's geometry
  REAL(SP), DIMENSION(3,3) :: BC_vertmat1,BC_vertmat2,  &
                               BC_simpmat1,BC_simpmat2  ! Matrices to check if point is
                                                        ! inside the 2 triangles making up the
                                                        ! quadrilateral BC.
  REAL(SP), DIMENSION(4) :: pl_eq                       ! equation for plane: ax+by+cz+d=0
  REAL(SP), DIMENSION(3) :: r_vec                       ! position vector in the plane, (basis1,basis2,1.0),
                                                        ! for simplex evaluation purposes.     
  REAL(SP), DIMENSION(3) :: s_vec1,s_vec2               ! 2D simplex co-ordinates
  INTEGER(I4B) :: inode                                 ! counters
  INTEGER(I4B) :: info                                  ! LAPACK
  INTEGER(I4B), DIMENSION(3) :: ipiv_local              ! LAPACK
  REAL(SP), DIMENSION(3) :: work_local                  ! LAPACK
  
  ! Establish an orthogonal basis to describe the plane:
  pl_norm = CROSS_PRODUCT(r2-r1,r4-r1)
  pl_norm = pl_norm/SQRT(SUM(pl_norm**2))
  pl_eq(1:3) = pl_norm                                ! 0 = ax + by + cz
  pl_eq(4) = - SUM(pl_norm*r1)                        ! + d (pl_eq = [a b c d])
  basis1 = (r1-r2)/SQRT(SUM((r1-r2)**2))
  basis2 = CROSS_PRODUCT(pl_norm,basis1)
 
  ! Use simplex co-ordinates to check whether a point falls within the two triangles
  ! that the quadrilateral comprises of.
 
  ! Setup the vertices matrices:
  BC_vertmat1(1,1) = SUM(r1*basis1)
  BC_vertmat1(1,2) = SUM(r2*basis1)
  BC_vertmat1(1,3) = SUM(r3*basis1)
  BC_vertmat1(2,1) = SUM(r1*basis2)
  BC_vertmat1(2,2) = SUM(r2*basis2)
  BC_vertmat1(2,3) = SUM(r3*basis2)
  BC_vertmat1(3,1:3) = 1.0
  BC_vertmat2(1,1) = SUM(r3*basis1)
  BC_vertmat2(1,2) = SUM(r4*basis1)
  BC_vertmat2(1,3) = SUM(r1*basis1)
  BC_vertmat2(2,1) = SUM(r3*basis2)
  BC_vertmat2(2,2) = SUM(r4*basis2)
  BC_vertmat2(2,3) = SUM(r1*basis2)
  BC_vertmat2(3,1:3) = 1.0

  ! Now invert the vertices matrices for simplex evaluation of the node locations:
  ! (Uses LAPACK routine SGETRF and SGETRI to invert.)
  BC_simpmat1 = BC_vertmat1 
  BC_simpmat2 = BC_vertmat2 
  ! BC_simpmat1:
  CALL SGETRF ((3),(3),BC_simpmat1,(3),ipiv_local,info)  
  ! Check error conditions - shouldn't occur.
  IF(info.NE.0) CALL ERROR_FEMFEKO(1,4032)
  CALL SGETRI ((3),BC_simpmat1,(3),ipiv_local,work_local,(3),info)
  IF(info.NE.0) CALL ERROR_FEMFEKO(1,4033)
  ! BC_simpmat2:
  CALL SGETRF ((3),(3),BC_simpmat2,(3),ipiv_local,info)  
  ! Check error conditions - shouldn't occur.
  IF(info.NE.0) CALL ERROR_FEMFEKO(1,4032)
  CALL SGETRI ((3),BC_simpmat2,(3),ipiv_local,work_local,(3),info)
  ! Check error conditions - also shouldn't occur.
  IF(info.NE.0) CALL ERROR_FEMFEKO(1,4033)
  ! Note that the matrices have now been overwritten with their inverses.

  ! Assign the flag values:
  node_flags = .FALSE.   ! intialise
  DO inode = 1,3         ! cycle through the 2 nodes 
    ! Check whether node is in the quad plane:
    IF (ABS(SUM(pl_eq(1:3)*vertices(faces(face_val)%nodes(inode))%coord(1:3))+pl_eq(4)).LT.EPS) THEN
      r_vec(1) = SUM(vertices(faces(face_val)%nodes(inode))%coord(1:3)*basis1)
      r_vec(2) = SUM(vertices(faces(face_val)%nodes(inode))%coord(1:3)*basis2)
      r_vec(3) = 1.0
      s_vec1 = MATMUL(BC_simpmat1,r_vec)
      s_vec2 = MATMUL(BC_simpmat2,r_vec)
      ! Check whether in the quad:
      IF ( (SUM(ABS(s_vec1))-1.0.LT.EPS) .OR. (SUM(ABS(s_vec2))-1.0.LT.EPS) ) THEN
        node_flags(inode) = .TRUE.
      END IF
    END IF
  END DO
  in_quad = node_flags(1).AND.node_flags(2).AND.node_flags(3)

END SUBROUTINE FACE_IN_QUADRILATERAL
!*******************************************************************************


SUBROUTINE FAST_CONNECT_ELEMENT_TO_FACE
  USE math_tools, ONLY: FIND_IN_LIST
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! The edge and face connections are found here, to be used in the construction 
! of the sparse matrix indices col_ind and row_ind.
! i.e: set up 
!  - faceconnectelem
! 
! Note: This routine requires that NODEELEMENT_INDEXLIST_MAKE was called
!       before hand. 
!
! Reworked MMB. 2001-10-10. Made this an order(N) process.
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem,this_node,this_elem,listpos,elempos

  ALLOCATE(faceconnectelem(num_faces,2))
  faceconnectelem = 0 

  FACE_LOOP: DO iface = 1,num_faces
    elempos = 0
    this_node = faces(iface)%nodes(1) ! Any one of the three nodes will
                                      ! be connected to both elements
                                      ! that share the face.
    ELEMENT_LOOP: DO ielem = Node_ptr(this_node)+1,Node_ptr(this_node+1)
      this_elem = Node_ind(ielem)
      CALL FIND_IN_LIST(elements(this_elem)%faces(1:4),4,iface,listpos)
      IF (listpos.GT.0) THEN ! this elemement is connected to face <iface>
        elempos = elempos + 1
        faceconnectelem(iface,elempos) = this_elem
      END IF
      IF (elempos.EQ.2) EXIT ELEMENT_LOOP
    END DO ELEMENT_LOOP
  END DO FACE_LOOP

END SUBROUTINE FAST_CONNECT_ELEMENT_TO_FACE
!*******************************************************************************


SUBROUTINE FAST_EDGE_CONNECT
   USE nrtype
   USE unit_numbers
   USE problem_info
   USE output_error, ONLY: ERROR_FEMFEKO
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine handles building edge interconnectivity data.
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! FJCM: Modification to DBDs EDGE_CONNECT
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version: 15-05-01
! 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the list edge nodes.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is face interconnectivity data in the (integer) array
!
!   edgeconnectelem
!
! which lists the elements that the edge is common to.
! 
! and the array 
!
!   edgeconnectedge
!
! which lists (by edge) all the other EDGES that are "connected".
!
! NB. Note that the search assumes that global edge numbers have been 
! assigned in increasing order in accordance with the element numbers. 
! 
! FJCM NOTE: col_start and row_subs not calculated (as in original edge_connect)
!            this must be done if FSS are to be considered (seems the only place)
!            where these variable arrays are used. 
!******************************************************************************* 
   INTEGER(I4B) ielem, jedge, jelem, kedge, ledge, row, knodes, column, counter 
   ! counters 
   INTEGER(I4B) previous_col, previous_elem, previous_edge! more counters 
   INTEGER(I4B) length ! array length
   INTEGER(I4B) connect_edge, connect_elem
   INTEGER(I4B) index, indx_j, t_col
   INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: temp
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: tempvec1
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: tempvec2
   LOGICAL(LGT) new_edge, jelem_already_linked
   INTEGER(I4B) ii,jj, key
   ! Sparsity info.
   INTEGER(I4B), DIMENSION(2) :: kedgenodes
 
   length = 2 ! An intial guess that will be corrected as needed.
   allocate(edgeconnectelem(num_edges,length))
 
   ! Initialize all connectivity data first (to zero, indicates no connection)
   edgeconnectelem =  0          ! Note array assignment
   ! Build up connectivity data. 
   OUTER_ELEMENT_LOOP: DO ielem=1,num_elements 
     OUTER_EDGE_LOOP: DO kedge = 1,6
       ! Locate relevant row of "edgeconnectelem"
       row = elements(ielem)%edges(kedge)
       ! 1st column is the "originating" element - lowest element # connected
       ! to this edge
       SEARCH_EDGE: IF (edgeconnectelem(row,1).EQ.0) THEN
         ! This edge has not previously been searched for.
         edgeconnectelem(row,1)  = ielem
         ! OPTIMIZER PRODUCES AN ERROR AT NEXT LINE (which disappears
         ! if the print statement is executed.)
         ! print *,edgeconnectelem(row,:),row,ielem
         column = 2

!        Note that here we don't search through all elements but only
!        through those elements with kfacenodes(knodes) as per node-element index 

         kedgenodes = GLOBAL_EDGENODES(ielem,kedge)
       
         INNER_NODEEDGE_LOOP: DO knodes = 1,2
           INNER_ELEMENT_LOOP:  DO indx_j = Node_ptr(kedgenodes(knodes))+1,Node_ptr(kedgenodes(knodes)+1)          
             jelem = Node_ind(indx_j)
             
             ! Test if this element has already been tested for common edges
             ! This is also a self test!!!
             jelem_already_linked = .false.
             JDONE_TEST1: DO t_col = 1, column-1  
               IF (edgeconnectelem(row,t_col).EQ.jelem) THEN
                 jelem_already_linked = .true.
                 EXIT JDONE_TEST1
               END IF  
             END DO JDONE_TEST1

             JDONE_TEST2: IF (.NOT.jelem_already_linked) THEN          
               INNER_EDGE_LOOP: DO ledge = 1,6
                 EDGE_TEST: IF (row.EQ.elements(jelem)%edges(ledge)) THEN
                 ! A common edge has been found
                 IF (column.GT.length) THEN ! Increase the memory allocation
                   length = length + 1
                   allocate(temp(num_edges,length))
                   temp = 0 ! Array allocation: ensure last column(s) zeroed.
                   temp(:,1:length-1) = edgeconnectelem(:,:)
                   deallocate(edgeconnectelem)
                   allocate(edgeconnectelem(num_edges,length))
                   edgeconnectelem = temp
                   deallocate(temp)
                 END IF
                 edgeconnectelem(row,column)  = jelem 
                 column = column+1 ! More elements may still be connected
                 EXIT INNER_EDGE_LOOP ! Edge can only appear once in an element
               END IF EDGE_TEST
             END DO INNER_EDGE_LOOP
             END IF JDONE_TEST2
           END DO INNER_ELEMENT_LOOP
         END DO INNER_NODEEDGE_LOOP
       ELSE ! SEARCH_EDGE
         CONTINUE
         ! This edge has already been searched as part of a lower # element
         ! Do not repeat search; go on to next edge.
       END IF SEARCH_EDGE
     END DO OUTER_EDGE_LOOP
   END DO OUTER_ELEMENT_LOOP


   length = 6 ! An intial guess that will be corrected as needed.     
              ! At least 6 essential.
   allocate(edgeconnectedge(num_edges,length))
   edgeconnectedge = 0 ! array
   IF (length.LT.6) THEN 
     STOP 'IE: Internal error in CONNECT. Increase length to at least 6 (edgeconnectedge)'
   END IF
   ! Build up a table of all the other EDGES which an edge could 
   ! possibly be connected to for generating non-zero matrix entries.
   DO row = 1,num_edges
     ! First, the edges comprising this element:
     edgeconnectedge(row,1:6) = elements(edgeconnectelem(row,1))%edges(1:6)
     counter = 7 
     DO column = 2,SIZE(edgeconnectelem,2)
       connect_elem = edgeconnectelem(row,column)
       IF (connect_elem.EQ.0) THEN
         EXIT ! No more connected elements to search through 
       END IF
       DO jedge = 1,6
         connect_edge = elements(connect_elem)%edges(jedge)
         new_edge = .true.
         ! Search ALL the previous elements to see if this edge has 
         ! previously appeared.
         DO previous_col = 1,column-1
           DO kedge = 1,6
             previous_elem = edgeconnectelem(row,previous_col)
             previous_edge = elements(previous_elem)%edges(kedge)
             IF (connect_edge.EQ.previous_edge) THEN
               new_edge = .false.
               EXIT ! This edge has previously appeared in another element
             END IF
           END DO
         END DO
         NEW_EDGE_TEST: IF (new_edge) THEN
           ! This edge has not previously appeared in another element
           IF (counter.GT.length) THEN ! Increase the memory allocation
             ! Note that counter may be more than 1 larger than length...
             allocate(temp(num_edges,counter))
             temp = 0 ! Array allocation: ensure last column(s) zeroed.
             temp(:,1:length) = edgeconnectedge(:,:)
             deallocate(edgeconnectedge)
             allocate(edgeconnectedge(num_edges,counter))
             edgeconnectedge = temp
             deallocate(temp)
             length = counter
           END IF
           edgeconnectedge(row,counter) = connect_edge               
           counter = counter +1
         END IF NEW_EDGE_TEST
       END DO
     END DO
   END DO

END SUBROUTINE FAST_EDGE_CONNECT     
!*******************************************************************************
       

SUBROUTINE FAST_EDGEMAKE
! NB NOTE: Not implemented for CBAA as was edgemake! 
! Better to use this routine to create adges and a permutation
! to sort aperture edges first, as required by CBAA
   USE nrtype
   USE problem_info
   USE unit_numbers
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! NOTE: Routine based on subroutine EDGEMAKE but using node-element list
!       to get O(N) for large probloms. 

! This subroutine makes the list of edges associated with the elements. 
! The edges are always made with the unit vector pointing from lower
! to higher node number; the elements need to have been thus sorted first. 
! See Lee and Mittra, "Edge-elements for 3-D inhomegoneously-filled cavities", 
! IEEE MTT, Sep 1992, pp. 1767-1773 for definitions.
!
! The local edge numbering scheme numbers the nodes as:
!  Edge in L& M  Edge in this code    Local edge number
!
!  e01            e12                 1
!  e02            e13                 2
!  e03            e14                 3
!  e12            e23                 4
!  e13            e24                 5
!  e23            e34                 6
! 
! The edge numbering scheme in this code is consistent with Savage and 
! Peterson [IEEE MTT June 96,see main code for full reference.]
!
! Global edge numbers are assigned from 1 upwards; within an element,
! glocal nodes are incremented in the same pattern as the local edges. 
! (i.e. if the next element shared edges e12 and e13, the next element
! would have global edge numbers 7,8,9 and 10 corresponding to local edges
! e14, e23, e24 and e34).
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! FJCM based on FACEMAKE of DBD
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! First version: 20.05.2001 -> FJCM
! Last editef:   20.05.2001 -> FJCM
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the sorted data structure "elements", passed 
! via MODULE geometry.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is a list of the global edge numbers, as well
! as a list of which global nodes a global edge connects.
!
! The data structures changed are:
!
!   elements%edges: 
!   %edges stores the global edges corresponding to the local edges, numbered
!   as above, per element.
! 
! and 
!
!   Tmp_edges%nodes:
!   %nodes  stores the global nodes that the edge connects, per edge.
!
! The number of edges (numedges) is also found by this routine.
!
! Note that this routine - AS DO SEVERAL OTHERS! - relies on the 
! nodes for each element being sorted into ascending order; this ensures
! the edge vectors are consistent from element to element.
!
!*******************************************************************************
    INTEGER(I4B) ielem,jelem,iedge,jedge,kedge,knodes,indx_j
    INTEGER(I4B), DIMENSION(2) ::  temp_edge
    LOGICAL(LGT) new_edge
    TYPE TMP_EDGE
      INTEGER(I4B), DIMENSION(2) :: nodes ! Global nodes
    END TYPE TMP_EDGE
    TYPE(TMP_EDGE), DIMENSION(:), ALLOCATABLE :: Tmp_edges

    ! Temporary edge data structure for building element connection data:
    ! Needed to minimize memory usage and to allocate the exact space
    ! needed for <edges>.)
    ALLOCATE (Tmp_edges(MAX_EDGES))        

    ! Find element edges:
    num_edges = 0

    ! Initialize all element edges, indicating edge not numbered yet 
    DO ielem=1,num_elements
      elements(ielem)%edges(:) = 0       
    END DO

    ELEMENT_LOOP: DO ielem=1,num_elements   ! Global element counter
      OUTER_EDGE_LOOP: DO iedge = 1,6       ! Local edge counter
        temp_edge = GLOBAL_EDGENODES(ielem,iedge) 
                      
        ! Search only already defined edges which are connected
        ! to ielem through node-element index list
        
        new_edge = .true.   
        INNER_NODEEDGE_LOOP: DO knodes = 1,2
          INNER_ELEMENT_LOOP: DO indx_j = Node_ptr(temp_edge(knodes))+1,Node_ptr(temp_edge(knodes)+1)          
          jelem = Node_ind(indx_j)
          
          ! Search all existing edges on connected elements         
          INNER_ELEMENT_EDGELOOP: DO jedge = 1,6 
            kedge = elements(jelem)%edges(jedge)
            IF (kedge.NE.0) THEN
            ! Search all existing edges on connected elements
              IF ( (temp_edge(1).EQ.Tmp_edges(kedge)%nodes(1)) .AND. &
                 (temp_edge(2).EQ.Tmp_edges(kedge)%nodes(2)) ) THEN 
                ! Existing edge found
                elements(ielem)%edges(iedge) = kedge
                new_edge = .false.
                EXIT INNER_NODEEDGE_LOOP ! No further search necessary
              END IF
            END IF   
          END DO INNER_ELEMENT_EDGELOOP
          END DO INNER_ELEMENT_LOOP
        END DO INNER_NODEEDGE_LOOP

        IF (new_edge) THEN 
          ! The search didn't find an existing edge, 
          ! so this is a new edge.
          num_edges=num_edges+1  
          elements(ielem)%edges(iedge) = num_edges
          Tmp_edges(num_edges)%nodes(1) = temp_edge(1)
          Tmp_edges(num_edges)%nodes(2) = temp_edge(2)
        END IF
      END DO OUTER_EDGE_LOOP
    END DO ELEMENT_LOOP

    CALL allocate_geometry(alloc_edges=num_edges)

    ! Transfer the data from TMP array:
    DO iedge=1,num_edges
      edges(iedge)%nodes(1:2) = TMP_edges(iedge)%nodes(1:2)
    END DO

    ! Deallocate temporary storage:
    DEALLOCATE(TMP_edges)
    
END SUBROUTINE FAST_EDGEMAKE
!*******************************************************************************


SUBROUTINE FAST_FACE_CONNECT
   USE nrtype
   USE unit_numbers
   USE problem_info
   USE output_error, ONLY: ERROR_FEMFEKO
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine builds both edge and face connectivity. It does it in a fast
! simple, but fast way by fisrt creating a node-element list and index. This 
! fast face connect routine should be of the order O(N). 
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! FJC Meyer
! Routine adapted from Metis mesh2graph routine
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 2 April 2001
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the list of nodes per element.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is face interconnectivity data in the (integer) arrays 
!   connect2elem   and connect2faces
! The following eg. shows how data is stored for a element i:
! connect2elem(i,1:4) 0 1 2 0 
! connect2face(i,1:4) 0 4 3 0
! The above eg. indicates that faces 1 and 4 are not connected;
! face 2 is connected to element 1's face 4;
! face 3 is connected to element 2's face 3.
!
! (At least) two different face numbering schemes are used in the literature.
! This code handles one:
!
! Savage and Peterson's convention:
! Face  Local node numbers
! 1         1 2 3
! 2         1 2 4 
! 3         1 3 4
! 4         2 3 4 
!*******************************************************************************
   INTEGER(I4B) ielem,jelem,kface,lface ! counters
   INTEGER(I4B) i,k, knodes,indx_j      ! more counters    
   INTEGER(I4B), DIMENSION(1:3) :: lfacenodes, kfacenodes
  
   ! Initialize all connectivity data first (to zero, indicates no connection)
   DO ielem=1,num_elements
     elements(ielem)%connect2elem(:) = 0 ! Note array assignment
     elements(ielem)%connect2face(:) = 0 !       ditto
   END DO

   OUTER_ELEMENT_LOOP: DO ielem=1,num_elements 
     OUTER_FACE_LOOP: DO kface = 1,4     
       IF (elements(ielem)%connect2elem(kface).NE.0) CYCLE OUTER_FACE_LOOP
       ! If a connection has already been made, no further search for this
       ! face is needed...   
       kfacenodes = GLOBAL_FACENODES(ielem,kface)
       ! Note that here we don't search through all elements but only
       ! through those elements with kfacenodes(knodes) as per node-element index 
       ! created above.        
       INNER_NODEFACE_LOOP: DO knodes = 1,3
         INNER_ELEMENT_LOOP: DO indx_j = Node_ptr(kfacenodes(knodes))+1,Node_ptr(kfacenodes(knodes)+1)         
          jelem = Node_ind(indx_j)
          ! No connection to self
          SELFFACE_TEST: IF (jelem.NE.ielem) THEN 
           INNER_FACE_LOOP: DO lface = 1,4     
             lfacenodes = GLOBAL_FACENODES(jelem,lface)
             
             NODE_TEST: IF ( (kfacenodes(1).eq.lfacenodes(1) ).AND.  & 
                             (kfacenodes(2).eq.lfacenodes(2) ).AND.  & 
                             (kfacenodes(3).eq.lfacenodes(3) ) ) THEN
               ! Note that this test relies on the nodes being sorted into
               ! ascending order for each element. 
               elements(ielem)%connect2elem(kface) = jelem 
               elements(ielem)%connect2face(kface) = lface
               elements(jelem)%connect2elem(lface) = ielem
               elements(jelem)%connect2face(lface) = kface  
               EXIT INNER_NODEFACE_LOOP
               ! Face "kface" of element "ielem" is  connected to element "jelem"
               ! and it is connected to that element's face "lface";
               ! and vice-versa. Then skip the rest of this face's search; 
               ! a face can only be connected to one other face -
               ! which has just been found.
             END IF NODE_TEST            
           END DO INNER_FACE_LOOP
          END IF SELFFACE_TEST
         END DO INNER_ELEMENT_LOOP
       END DO INNER_NODEFACE_LOOP
     END DO OUTER_FACE_LOOP
   END DO OUTER_ELEMENT_LOOP
     
END SUBROUTINE FAST_FACE_CONNECT   
!*******************************************************************************


SUBROUTINE FAST_FACEMAKE
   USE nrtype
   USE problem_info
   USE unit_numbers   
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! NOTE: Routine based on subroutine FACEMAKE but using node-element list
!       to get O(N) for large probloms. 
! 
! This subroutine makes the list of faces associated with the elements. 
! It is assumed - as for the EDGEMAKE routine - that the elements 
! have been sorted first with local node 1 corresponding to lowest global 
! node and the rest sorted in ascending order. This is ESSENTIAL for the
! searching algorithm for new faces to work properly.
! 
! The face numbering scheme in this code is consistent with Savage and 
! Peterson [IEEE MTT June 96,see main code for full reference.] (This is 
! actually established in routine GLOBAL_FACENODES). The local numbering 
! scheme is:
!
! Face   Node
!   1   1 2 3
!   2   1 2 4
!   3   1 3 4
!   4   2 3 4
!
! Global face numbers are assigned from 1 upwards; within an element,
! global faces are incremented in the same pattern as the local faces. 
! (i.e. if the second element shared faces 1 and 2 with the first element, 
! the next element would have global face numbers 5 and 6 
! corresponding to local faces 3 and 4 of element 2).
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! FJCM based on FACEMAKE of DBD
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! First version: 20.05.2001 -> FJCM
! Last editef:   20.05.2001 -> FJCM
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the sorted data structures "elements", passed via MODULE
! geometry.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the global face information.
! The data structures are:
!
!   - elements%faces
!   %faces stores the global face corresponding to the local face numbered
!   as above, per element.
!   - Tmp_faces%nodes
!   %nodes stores the global nodes that the face contains, per face.
! Again, the searching algorithms are quite possibly sub-optimal for 
! very large meshes.
!
!*******************************************************************************
    INTEGER(I4B) ielem,jelem,iface,jface,kface,knodes,indx_j
    INTEGER(I4B), DIMENSION(3) ::  temp_face
    LOGICAL(LGT) new_face
    TYPE TMP_FACE
     INTEGER(I4B), DIMENSION(3) :: nodes ! Global nodes
    END TYPE TMP_FACE
    TYPE(TMP_FACE), DIMENSION(:), ALLOCATABLE :: Tmp_faces
 
    ! Allocate temporary storage for face data to minimize memory usage and
    ! to allocato only the exact number of <faces>:
    ALLOCATE (Tmp_faces(MAX_FACES))   

    ! Initialize all element faces, indicating face not numbered yet 
    DO ielem=1,num_elements
      elements(ielem)%faces(:) = 0       
    END DO  

    ! Find element faces
    num_faces = 0

    ELEMENT_LOOP: DO ielem=1,num_elements   ! Global element counter
      OUTER_FACE_LOOP: DO iface = 1,4       ! Local face counter
        temp_face = GLOBAL_FACENODES(ielem,iface) 

        ! temporary face now contains the global nodes for this element...
        FIRST_ELEM_TEST: IF (ielem.EQ.1) THEN   
          ! First element requires special treatment:
          ! all faces are new in this case.
          ! All faces for this element have vector sense +1.
          num_faces=num_faces+1  
          elements(ielem)%faces(iface) = num_faces 
          Tmp_faces(num_faces)%nodes = temp_face ! Array assignment
        ELSE 
          ! Second and subsequent elements...
          new_face = .true.

          ! Search only already defined faces which are connected
          ! to ielem through node-element index list
           
          INNER_NODEFACE_LOOP: DO knodes = 1,3
            INNER_ELEMENT_LOOP: DO indx_j = Node_ptr(temp_face(knodes))+1,Node_ptr(temp_face(knodes)+1)        
            jelem = Node_ind(indx_j)
            
            ! Search all existing faces on connected elements           
            INNER_ELEMENT_FACELOOP: DO jface = 1,4 
              kface = elements(jelem)%faces(jface)
              IF (kface.NE.0) THEN
                IF ((temp_face(1).EQ.Tmp_faces(kface)%nodes(1)) .AND. &
                   (temp_face(2).EQ.Tmp_faces(kface)%nodes(2)) .AND. & 
                   (temp_face(3).EQ.Tmp_faces(kface)%nodes(3))       &                           
                  ) THEN 
                  ! Existing face found
                  elements(ielem)%faces(iface) = kface
                  new_face = .false.
                  EXIT INNER_NODEFACE_LOOP ! No further search necessary
                END IF
              END IF
            END DO INNER_ELEMENT_FACELOOP  
            END DO INNER_ELEMENT_LOOP
          END DO INNER_NODEFACE_LOOP
          IF (new_face) THEN 
            ! The search didn't find an existing face, 
            ! so this is a new face.
            num_faces=num_faces+1  ! increment number of faces
            elements(ielem)%faces(iface) = num_faces         
            Tmp_faces(num_faces)%nodes = temp_face ! Array assignment
         END IF
        END IF FIRST_ELEM_TEST
      END DO OUTER_FACE_LOOP
    END DO ELEMENT_LOOP

    ! Allocate the exact no of faces
    CALL allocate_geometry(alloc_faces=num_faces) 

    ! Transfer the data from TMP array
    DO iface = 1,num_faces
      faces(iface)%nodes(1:3) = TMP_faces(iface)%nodes(1:3)
    END DO
  
    ! Deallocate the temporary storage:
    DEALLOCATE(TMP_faces)

END SUBROUTINE FAST_FACEMAKE
!*******************************************************************************


SUBROUTINE FREE_ASSIGN_SEARCH
  USE nrtype
  USE eigen_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! Assigns the free flag for all faces and edges according to their PEC flag:
! free = .NOT.PEC
!*******************************************************************************
  INTEGER(I4B) :: ielem,iedge,iface
  INTEGER(I4B), DIMENSION(2) :: temp_edge_nodes

  ! Assign the edges' free property:
  DO iedge = 1,num_edges
    IF (edges(iedge)%PEC) THEN
      edges(iedge)%free = .FALSE.

! DBD change 4 Dec 02 and again 13 Dec
      IF (PREDICT_SPURIOUS_EIGENMODES) THEN
        vertices(edges(iedge)%nodes(1))%free = .FALSE.
        vertices(edges(iedge)%nodes(2))%free = .FALSE.
      END IF

    END IF       
  END DO
! End DBD change 4 Dec 02

  ! Assigns the faces' free property:
  DO iface = 1,num_faces
    IF (faces(iface)%PEC) THEN
      faces(iface)%free = .FALSE.
    END IF       
  END DO

END SUBROUTINE FREE_ASSIGN_SEARCH
!*******************************************************************************

SUBROUTINE FREE_DIR_ASSIGN_SEARCH
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Assigns the free flag for all faces and edges according to their Dirchlet flag:
! free = .NOT.PEC and also .NOT.Dirchlet. Not used for eigenvalue computations.
! DBD 5 Aug 04.
!*******************************************************************************
  INTEGER(I4B) :: ielem,iedge,iface
  INTEGER(I4B), DIMENSION(2) :: temp_edge_nodes

  ! Assign the edges' free property:
  DO iedge = 1,num_edges
    IF (edges(iedge)%PEC.OR.edges(iedge)%Dirichlet) THEN
      edges(iedge)%free = .FALSE.
    END IF       
  END DO

  ! Assigns the faces' free property:
  DO iface = 1,num_faces
    IF (faces(iface)%PEC.OR.faces(iface)%Dirichlet) THEN
      faces(iface)%free = .FALSE.
    END IF       
  END DO

END SUBROUTINE FREE_DIR_ASSIGN_SEARCH
!*******************************************************************************

SUBROUTINE geometry_init
  USE nrtype
  IMPLICIT NONE

  ! Added DBD 3 June and 11 July 03
  PML_PRESENT  = .FALSE.
  PML%absorb_x = .FALSE.
  PML%absorb_y = .FALSE.
  PML%absorb_z = .FALSE.
  PML%full     = .TRUE.
  ! End DBD 3 June and 11 July 03

  num_nodes     = 0
  num_elements  = 0
  num_edges     = 0
  num_faces     = 0
  num_BCs       = 0

  PML_PRESENT = .FALSE.

  avg_edge_length = 0

END SUBROUTINE geometry_init

FUNCTION get_meshinfo()
!!!
!!! Output variable
!!!    
  TYPE(meshinfo_type) :: get_meshinfo
!!!
!!! Local Variables
!!!
  INTEGER(i4b) :: i

  get_meshinfo%num_nodes    = num_nodes    
  get_meshinfo%num_elements = num_elements 
  get_meshinfo%num_PEC_edges = 0 ! Will count them below
  get_meshinfo%num_edges    = num_edges    
  get_meshinfo%num_faces    = num_faces    
  get_meshinfo%num_PEC_faces = 0 ! Will count them below
  get_meshinfo%num_BCs      = num_BCs      

  DO i=1, num_faces             ! Count PEC faces
     IF (faces(i)%PEC) THEN 
        get_meshinfo%num_PEC_faces = get_meshinfo%num_PEC_faces + 1
     END IF
  END DO

  DO i=1, num_edges             ! Count PEC faces
     IF (edges(i)%PEC) THEN 
        get_meshinfo%num_PEC_edges = get_meshinfo%num_PEC_edges + 1
     END IF
  END DO
  

END FUNCTION get_meshinfo

FUNCTION GLOBAL_EDGENODES(elemnum,edgenum)
  USE nrtype
  IMPLICIT NONE
!******************************************************************************
! A function to return the two GLOBAL nodes associated with edges 1-6.
! Note that the edges 1-6 correspond to e12,e13,e14,e23,e24 and e34 respectively.
! DBD 1997. Revised Jan-Feb 1999.
!******************************************************************************
  INTEGER(I4B), INTENT(IN) :: edgenum,elemnum 
  INTEGER(I4B), DIMENSION(2) :: GLOBAL_EDGENODES

  SELECT CASE(edgenum)
  CASE(1) ! e12
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(1)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(2)
  CASE(2) ! e13
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(1)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(3)
  CASE(3) ! e14
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(1)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(4)
  CASE(4) ! e23
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(2)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(3)
  CASE(5) ! e24
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(2)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(4)
  CASE(6) ! e34
    GLOBAL_EDGENODES(1) = elements(elemnum)%nodes(3)
    GLOBAL_EDGENODES(2) = elements(elemnum)%nodes(4)
  CASE DEFAULT
    STOP 'IE: Edge no. out of range in function GLOBAL_EDGENODES.'
  END SELECT

! Better, valid code, using array constructors that would not compile 
! on Salford FTN95 compiler: 
!   SELECT CASE(edgenum)
!   CASE(1) ! e12
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(1), elements(elemnum)%nodes(2)/)
!   CASE(2) ! e13
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(1), elements(elemnum)%nodes(3)/) 
!   CASE(3) ! e14
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(1), elements(elemnum)%nodes(4)/)
!   CASE(4) ! e23
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(2), elements(elemnum)%nodes(3)/)
!   CASE(5) ! e24
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(2), elements(elemnum)%nodes(4)/)
!   CASE(6) ! e34
!     GLOBAL_EDGENODES = (/elements(elemnum)%nodes(3), elements(elemnum)%nodes(4)/)
!   CASE DEFAULT
!     STOP 'Error in function GLOBAL_EDGENODES'
!   END SELECT

END FUNCTION GLOBAL_EDGENODES
!*******************************************************************************


FUNCTION GLOBAL_FACENODES(elemnum,facenum)
   USE nrtype
   IMPLICIT NONE
!******************************************************************************
! A function to return the three GLOBAL nodes associated with faces 1-4, 
! using the convention (after Savage and Peterson)
! Face   Node
!   1   1 2 3
!   2   1 2 4
!   3   1 3 4
!   4   2 3 4
! DBD Jan-Feb 1999.
!******************************************************************************
   INTEGER(I4B), INTENT(IN) :: facenum,elemnum
   INTEGER(I4B), DIMENSION(3) :: GLOBAL_FACENODES
  
   SELECT CASE(facenum)
   CASE(1)
     GLOBAL_FACENODES(1) = elements(elemnum)%nodes(1)
     GLOBAL_FACENODES(2) = elements(elemnum)%nodes(2) 
     GLOBAL_FACENODES(3) = elements(elemnum)%nodes(3)
   CASE(2)
     GLOBAL_FACENODES(1) = elements(elemnum)%nodes(1)
     GLOBAL_FACENODES(2) = elements(elemnum)%nodes(2) 
     GLOBAL_FACENODES(3) = elements(elemnum)%nodes(4)
   CASE(3)
     GLOBAL_FACENODES(1) = elements(elemnum)%nodes(1)
     GLOBAL_FACENODES(2) = elements(elemnum)%nodes(3) 
     GLOBAL_FACENODES(3) = elements(elemnum)%nodes(4)
   CASE(4)
     GLOBAL_FACENODES(1) = elements(elemnum)%nodes(2)
     GLOBAL_FACENODES(2) = elements(elemnum)%nodes(3) 
     GLOBAL_FACENODES(3) = elements(elemnum)%nodes(4)
   CASE DEFAULT
     STOP 'IE: Local facenum. out of range in function GLOBAL_FACENODES.'
    END SELECT
 
END FUNCTION GLOBAL_FACENODES
!*******************************************************************************


FUNCTION GRADIENT_LAMBDA(element_num,normalize)
   USE math_tools, ONLY: DET_DIM3
   USE nrtype
   USE problem_info
   USE unit_numbers
   IMPLICIT NONE
!*******************************************************************************
! A function to return the gradient of the simplex coordinates
! for the element number specified. 
! The gradient should usually be normalized. However, FE codes often 
! require subsequent multiplication by the volume, so the option is 
! provided NOT to carry this out, permitting the volume to be exactly
! cancelled.
! DBD 1997. Revised Jan-Feb 1999.
! MMB 18 Sept 2001, changed to the gradients of all 4 simplex coords instead of only 1.
!*******************************************************************************
   INTEGER(I4B), INTENT(IN)   :: element_num
   LOGICAL(LGT), INTENT(IN)   :: normalize
   REAL(SP), DIMENSION(3,3)   ::  x, y, z
   REAL(SP), DIMENSION(4,3)     :: GRADIENT_LAMBDA
   REAL(SP), DIMENSION(4,3) :: nabla_lambda
   REAL(SP) x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4
 
   ! Get nodal values for nodes 1 to 4 of element number element_num
   x1 = vertices(elements(element_num)%nodes(1))%coord(1)
   x2 = vertices(elements(element_num)%nodes(2))%coord(1)
   x3 = vertices(elements(element_num)%nodes(3))%coord(1)
   x4 = vertices(elements(element_num)%nodes(4))%coord(1)
   y1 = vertices(elements(element_num)%nodes(1))%coord(2)
   y2 = vertices(elements(element_num)%nodes(2))%coord(2)
   y3 = vertices(elements(element_num)%nodes(3))%coord(2)
   y4 = vertices(elements(element_num)%nodes(4))%coord(2)
   z1 = vertices(elements(element_num)%nodes(1))%coord(3)
   z2 = vertices(elements(element_num)%nodes(2))%coord(3)
   z3 = vertices(elements(element_num)%nodes(3))%coord(3)
   z4 = vertices(elements(element_num)%nodes(4))%coord(3)

   ! Work out the gradient of the simplex coordinates
   ! This is given by:
   !                            ^  ^   ^ 
   !   __               +1   | 0 x  y  z  | 
   !   \/ lambda_1  =   --   | 1 x2 y2 z2 |
   !                    6V   | 1 x3 y3 z3 | 
   !                         | 1 x4 y4 z4 |
   !
   ! and similar expressions for the others (the signs alternate):
   !                            ^  ^   ^ 
   !   __               -1   | 0 x  y  z  | 
   !   \/ lambda_2  =   --   | 1 x1 y1 z1 |
   !                    6V   | 1 x3 y3 z3 | 
   !                         | 1 x4 y4 z4 |
   ! 
   ! Note there that the volume here may be SIGNED!
 
   x(1,:)=(/1.0_SP, y2,z2/) 
   x(2,:)=(/1.0_SP, y3,z3/) 
   x(3,:)=(/1.0_SP, y4,z4/) 

   y(1,:)=(/1.0_SP, x2,z2/) 
   y(2,:)=(/1.0_SP, x3,z3/) 
   y(3,:)=(/1.0_SP, x4,z4/) 

   z(1,:)=(/1.0_SP, x2,y2/) 
   z(2,:)=(/1.0_SP, x3,y3/) 
   z(3,:)=(/1.0_SP, x4,y4/) 
!!$   PRINT*, '1:'
!!$   PRINT*, 'x:'
!!$   PRINT*, x(1,:)
!!$   PRINT*, x(2,:)
!!$   PRINT*, x(3,:)
!!$   PRINT*, 'y:'
!!$   PRINT*, y(1,:)
!!$   PRINT*, y(2,:)
!!$   PRINT*, y(3,:)
!!$   PRINT*, 'z:'
!!$   PRINT*, z(1,:)
!!$   PRINT*, z(2,:)
!!$   PRINT*, z(3,:)

   nabla_lambda(1,1) = - DET_DIM3(x)
   nabla_lambda(1,2) = + DET_DIM3(y)
   nabla_lambda(1,3) = - DET_DIM3(z)

   x(1,:)=(/1.0_SP, y1,z1/)  ! Others unchanged
   y(1,:)=(/1.0_SP, x1,z1/) 
   z(1,:)=(/1.0_SP, x1,y1/) 
 
!!$   PRINT*, '2:'
!!$   PRINT*, 'x:'
!!$   PRINT*, x(1,:)
!!$   PRINT*, x(2,:)
!!$   PRINT*, x(3,:)
!!$   PRINT*, 'y:'
!!$   PRINT*, y(1,:)
!!$   PRINT*, y(2,:)
!!$   PRINT*, y(3,:)
!!$   PRINT*, 'z:'
!!$   PRINT*, z(1,:)
!!$   PRINT*, z(2,:)
!!$   PRINT*, z(3,:)

   nabla_lambda(2,1) = + DET_DIM3(x)
   nabla_lambda(2,2) = - DET_DIM3(y)
   nabla_lambda(2,3) = + DET_DIM3(z)

   x(2,:)=(/1.0_SP, y2,z2/) ! Others unchanged
   y(2,:)=(/1.0_SP, x2,z2/) 
   z(2,:)=(/1.0_SP, x2,y2/) 
!!$
!!$   PRINT*, '3:'
!!$   PRINT*, 'x:'
!!$   PRINT*, x(1,:)
!!$   PRINT*, x(2,:)
!!$   PRINT*, x(3,:)
!!$   PRINT*, 'y:'
!!$   PRINT*, y(1,:)
!!$   PRINT*, y(2,:)
!!$   PRINT*, y(3,:)
!!$   PRINT*, 'z:'
!!$   PRINT*, z(1,:)
!!$   PRINT*, z(2,:)
!!$   PRINT*, z(3,:)

   nabla_lambda(3,1) = - DET_DIM3(x)
   nabla_lambda(3,2) = + DET_DIM3(y)
   nabla_lambda(3,3) = - DET_DIM3(z)

   x(3,:)=(/1.0_SP, y3,z3/) ! Others unchanged
   y(3,:)=(/1.0_SP, x3,z3/) 
   z(3,:)=(/1.0_SP, x3,y3/) 

!!$   PRINT*, '4:'
!!$   PRINT*, 'x:'
!!$   PRINT*, x(1,:)
!!$   PRINT*, x(2,:)
!!$   PRINT*, x(3,:)
!!$   PRINT*, 'y:'
!!$   PRINT*, y(1,:)
!!$   PRINT*, y(2,:)
!!$   PRINT*, y(3,:)
!!$   PRINT*, 'z:'
!!$   PRINT*, z(1,:)
!!$   PRINT*, z(2,:)
!!$   PRINT*, z(3,:)

   nabla_lambda(4,1) = + DET_DIM3(x)
   nabla_lambda(4,2) = - DET_DIM3(y)
   nabla_lambda(4,3) = + DET_DIM3(z)

   ! Normalize, if required.
   IF (normalize) THEN
     nabla_lambda = nabla_lambda/(6.0_SP*SIGNED_VOLUME(element_num))
   END IF

   GRADIENT_LAMBDA = nabla_lambda

END FUNCTION GRADIENT_LAMBDA
!******************************************************************************


SUBROUTINE INIT_GEOMETRY_ATTRIBUTES
  USE eigen_analysis_data
  IMPLICIT NONE
!******************************************************************************
! Initializes the properties of the edges and faces to their default values.
!******************************************************************************

  ! (Note array assignments.)
  edges%free          = .TRUE.    ! General
  edges%PEC           = .FALSE.   ! General
  edges%CBAA_aperture = .FALSE.   ! CBAA only
  edges%coax_aperture = .FALSE.
  edges%port          = .FALSE.   ! WG only
  edges%portnumber    = 0         ! WG only (0 implies not a port.) 
! Added DBD 2 April 2003
  edges%ABC           = .FALSE.    ! TD only
  edges%ABCnumber     = 0          ! TD only (0 implies not an ABC.) 
! End added DBD 2 April 2003
! Added DBD 5 Aug 2004
  edges%Dirichlet      = .FALSE.   ! Scattering only
! End added DBD 5 Aug 2004

  faces%free          = .TRUE.    ! General
  faces%PEC           = .FALSE.   ! General
  faces%CBAA_aperture = .FALSE.   ! CBAA only
  faces%coax_aperture = .FALSE.
  faces%port          = .FALSE.   ! WG only
  faces%portnumber    = 0         ! WG only (0 implies not a port.) 
! Added DBD 2 April 2003
  faces%ABC            = .FALSE.  ! TD only
  faces%ABCnumber      = 0        ! TD only (0 implies not an ABC.) 
! End added DBD 2 April 2003
! Added DBD 5 Aug 2004
  faces%Dirichlet      = .FALSE.   ! Scattering only
! End added DBD 5 Aug 2004
! Added DBD 4 Aug 2005
  faces%curvilinear    = .FALSE.   ! Scattering only
! End added DBD 4 Aug 2005


  IF ((REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS).AND.&
    PREDICT_SPURIOUS_EIGENMODES) THEN
    vertices%free = .TRUE. 
  END IF

END SUBROUTINE INIT_GEOMETRY_ATTRIBUTES
!******************************************************************************

FUNCTION JACOBIAN_POLY2(element_num,u,v,w,inverse)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! A function to return the Jacobian of a 2nd order polynomial curvilinear element, element number
! element_num, with local coordinates (u,v,w). The numbering scheme of mid-point nodes is defined 
! as in the header of this module.
! Via the optional flag inverse, the inverse can be requested instead of the Jabobian. 
! Written from scratch, DBD 13 July 2005.
!
!*******************************************************************************
  SAVE ! Much of this need not be re-computed for the same element
  INTEGER(I4B), INTENT(IN)   :: element_num
  REAL(SP), INTENT(IN)   :: u,v,w
  LOGICAL(LGT), INTENT(IN)  :: inverse
  
  REAL(SP), DIMENSION(3,3) :: JACOBIAN_POLY2,JJ
  REAL(SP), DIMENSION(10) :: dN_by_du, dN_by_dv, dN_by_dw
  LOGICAL(LGT):: first_call = .true. ! This is SAVEd by the overall SAVE statement.
  INTEGER(I4B)   :: last_element_num
  INTEGER(I4B) iedge,inode,iface

  REAL(SP), DIMENSION(10)   ::  x, y, z
  INTEGER(I4B) :: info               ! Info from LAPACK factorization.
  INTEGER(I4B) :: lwork              ! LAPACK dimensioning data
  PARAMETER(lwork=3)
  REAL(SP), DIMENSION(lwork) :: temp_work ! LAPACK temporary workspace.
  INTEGER(I4B), DIMENSION(3) :: temp_ipiv ! Pivot indices from factorization. 
!  REAL(SP) :: phi_1, phi_2, theta_1, theta_2
  
  IF(element_num.NE.last_element_num.OR.first_call) THEN
    CALL MID_EDGE_NODES(element_num,x,y,z)
  END IF

  dN_by_du = (/ 1.0_SP-4.0_SP*(1.0_SP-u-v-w) , 4.0_SP*u-1.0_SP , 0.0_SP , 0.0_SP , 4.0_SP*(1.0_SP-2.0_SP*u-v-w) , & 
                -4.0_SP*v, -4.0_SP*w, 4.0_SP*v, 4.0_SP*w, 0.0_SP /)
  dN_by_dv = (/ 1.0_SP-4.0_SP*(1.0_SP-u-v-w) , 0.0_SP, 4.0_SP*v-1.0_SP , 0.0_SP , -4.0_SP*u , & 
                4.0_SP*(1.0_SP-u-2.0_SP*v-w), -4.0_SP*w, 4.0_SP*u, 0.0_SP, 4.0_SP*w /)
  dN_by_dw = (/ 1.0_SP-4.0_SP*(1.0_SP-u-v-w) , 0.0_SP, 0.0_SP, 4.0_SP*w-1.0_SP , -4.0_SP*u , & 
                -4.0_SP*v, 4.0_SP*(1.0_SP-u-v-2.0_SP*w), 0.0_SP, 4.0_SP*u, 4.0_SP*v /)

  JJ = 0.0 ! Initialize on each call.
  DO inode = 1,10
    JJ(1,1) = JJ(1,1)+dN_by_du(inode)*x(inode)
    JJ(1,2) = JJ(1,2)+dN_by_du(inode)*y(inode)
    JJ(1,3) = JJ(1,3)+dN_by_du(inode)*z(inode)
    JJ(2,1) = JJ(2,1)+dN_by_dv(inode)*x(inode)
    JJ(2,2) = JJ(2,2)+dN_by_dv(inode)*y(inode)
    JJ(2,3) = JJ(2,3)+dN_by_dv(inode)*z(inode)
    JJ(3,1) = JJ(3,1)+dN_by_dw(inode)*x(inode)
    JJ(3,2) = JJ(3,2)+dN_by_dw(inode)*y(inode)
    JJ(3,3) = JJ(3,3)+dN_by_dw(inode)*z(inode)
  END DO
  IF (.NOT.inverse) THEN
    CONTINUE
  ELSE IF(inverse) THEN
    ! Invert matrix. Uses LAPACK routine SGETRF to factorize, and then 
    ! SGETRI to invert. Note - THESE ARE SINGLE PRECISION ROUTINES!! 
    ! See note above!!
    CALL SGETRF ((3),(3),JJ,(3),temp_ipiv,info)  
    ! Check error conditions - shouldn't occur.
    IF(Info.NE.0) THEN 
    WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRF in subroutine JACOBIAN_POLY2.'
    END IF
    CALL SGETRI ((3),JJ,(3),temp_ipiv,temp_work,lwork,info)
    ! Check error conditions - also shouldn't occur.
    ! J has now been overwritten with its inverse
    IF(Info.NE.0) THEN 
      WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRI in subroutine JACOBIAN_POLY2.'
    END IF
  END IF
  JACOBIAN_POLY2 = JJ ! Returns either Jacobian or inverse thereof.
  last_element_num = element_num
  first_call = .FALSE.
END FUNCTION JACOBIAN_POLY2


SUBROUTINE MID_EDGE_NODES(element_num,x,y,z)
  USE nrtype
  USE scattering_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! A function to return all the nodes associated with a 2nd order polynomial curvilinear element, element number
! element_num. The numbering scheme of mid-point nodes is defined as in the header of this module.
!
! The mid-point nodes are presently either generated by linear interpolation of the vertex coordinates 
! (i.e. effectively providing a rectilinear element, used for tested the curvilinear schemes) 
! or by using the radius of a scatterer as specified in the input file.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN)   :: element_num
  REAL(SP), DIMENSION(10), INTENT(OUT)   ::  x, y, z ! Nodal cofficients.
  INTEGER(I4B) iedge,inode,iface
  REAL(SP) :: r1, r2, r_mid
  REAL(SP) :: phi_mid, theta_mid
  REAL(SP), DIMENSION(3) :: temp_coords
  INTEGER(I4B), DIMENSION(3) :: faceedges 
  INTEGER(I4B), DIMENSION(2) :: edgenodes

  ! Get the (x,y,z) coordinates of, firstly, the vertices...
  DO inode = 1,4 ! These are the  vertex coordindates
    x(inode) = vertices(elements(element_num)%nodes(inode))%coord(1)
    y(inode) = vertices(elements(element_num)%nodes(inode))%coord(2)
    z(inode) = vertices(elements(element_num)%nodes(inode))%coord(3)
  END DO
  ! Now, compute the mid-edge nodes by linear interpolation:
  DO iedge = 1,6
    inode = iedge+4
    temp_coords = EDGE_CENTRE(element_num,iedge)
    x(inode) = temp_coords(1)
    y(inode) = temp_coords(2)
    z(inode) = temp_coords(3)
  END DO
  
  ! If required, improve the mid-edge node if the scatterer is known to be spherical.
  IF (CURVILINEAR.AND.PEC_SPHERE_SCAT) THEN
    DO iface = 1,4
      faceedges = LOCAL_FACEEDGES(iface)
      IF(faces(elements(element_num)%faces(iface))%curvilinear) THEN
        EDGE_LOOP: DO iedge = 1,3
          edgenodes = GLOBAL_EDGENODES(element_num,faceedges(iedge))          
          r1 = SQRT((vertices(edgenodes(1))%coord(1))**2+(vertices(edgenodes(1))%coord(2))**2+(vertices(edgenodes(1))%coord(3))**2)
          r2 = SQRT((vertices(edgenodes(2))%coord(1))**2+(vertices(edgenodes(2))%coord(2))**2+(vertices(edgenodes(2))%coord(3))**2)
          IF (ABS(r1-sph_radius).GT.0.001_SP*sph_radius.OR.ABS(r2-sph_radius).GT.0.001_SP*sph_radius) THEN ! Double-check...
            STOP 'IE IN ROUTINE JACOBIAN_POLY2; RADIUS OF A CURVILINEAR ELEMENT IS INCORRECT'
          END IF
! Following "exact" algorithm has quadrant ambiguity problems...
!		   phi_1    = ATAN2(vertices(edgenodes(1))%coord(2),vertices(edgenodes(1))%coord(1))
!		   phi_2    = ATAN2(vertices(edgenodes(2))%coord(2),vertices(edgenodes(2))%coord(1))
!		   theta_1  = ACOS(vertices(edgenodes(1))%coord(3)/sph_radius)
!		   theta_2  = ACOS(vertices(edgenodes(2))%coord(3)/sph_radius)
!		   inode    = faceedges(iedge)+4
!    	   x(inode) = SIN((theta_1+theta_2)/2.0_SP)*COS((phi_1+phi_2)/2.0_SP)*sph_radius
!    	   y(inode) = SIN((theta_1+theta_2)/2.0_SP)*SIN((phi_1+phi_2)/2.0_SP)*sph_radius
!    	   z(inode) = COS((theta_1+theta_2)/2.0_SP)*sph_radius

! This algorithm is possibly approximate, but does not have quadrant ambiguity problems. It uses the angles phi and theta to the mid-point 
! coordinates, and then corrects the nodal coordinates for the radius.
		  inode     = faceedges(iedge)+4
          phi_mid   = ATAN2(y(inode),x(inode))
          theta_mid = ACOS(z(inode)/sph_radius)
    	  x(inode) = SIN(theta_mid)*COS(phi_mid)*sph_radius
    	  y(inode) = SIN(theta_mid)*SIN(phi_mid)*sph_radius
    	  z(inode) = COS(theta_mid)*sph_radius
! WRITE(FILEOUT,*) 'Mid point nodes of element', element_num, ' node ', inode, '(x,y,z)=',x(inode),y(inode),z(inode)
            ! Double-check calculation
		  r_mid = SQRT(x(inode)**2 + y(inode)**2 + z(inode)**2)
		  IF (ABS(r_mid-sph_radius).GT.0.001_SP*sph_radius) THEN 
            STOP 'IE IN ROUTINE JACOBIAN_POLY2; RADIUS OF MID-POINT NODE ON CURVILINEAR ELEMENT INCORRECT'
          END IF
        END DO EDGE_LOOP 
	  END IF
	END DO
  END IF
END SUBROUTINE MID_EDGE_NODES

!******************************************************************************


! DBD addition 21 May 2003
SUBROUTINE SCAT_TOTAL_BOUNDARY_SEARCH
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! This routine finds faces lying on the boundary between the internal inhomegeneous materials
! (material labels not equal to HOMOG_MEDIUM) and external homogeneous material label HOMOG_MEDIUM
! and flags them appropriately.
! Note that this should be performed BEFORE PML_ASSIGN_PROPERTIES is called!
! (This routine overwrites some of the external homogeneous material label with PML-specific 
! properties. 
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem
  INTEGER(I4B), DIMENSION(3) :: tempedges 

  ! Initialize
  edges%scat_tot_boundary = .FALSE.
  faces%scat_tot_boundary = .FALSE.

  ELEMENT_LOOP: DO ielem = 1,num_elements

    IF (elements(ielem)%material.NE.HOMOG_MEDIUM) THEN
	  ! This is in the inhomogeneous target region. Test if it lies on the (fictitious) boundary between
	  ! the scattered/total field region. 
      DO iface = 1,4
        tempedges = LOCAL_FACEEDGES(iface)
        IF(elements(ielem)%connect2elem(iface).EQ.0) THEN
		  WRITE(FILEOUT,'(A,I8,2X,A,I2,2X,A,I2)')'Possible meshing error: unconnected internal element ',ielem, 'face ',iface, & 
		                                 'material ',elements(ielem)%material
        ELSE IF (elements(elements(ielem)%connect2elem(iface))%material.EQ.HOMOG_MEDIUM) THEN
		  ! This element is connected to an element in the homogeneous region 
		  ! Flag appropriately.
	      edges(elements(ielem)%edges(tempedges))%scat_tot_boundary = .TRUE.
          faces(elements(ielem)%faces(iface))%scat_tot_boundary = .TRUE.
		END IF
		! Note that it is possible for an element to have MORE than one face on this boundary (eg in a corner)
		! so the search must proceed over all faces. 
	  END DO
    END IF    
  END DO ELEMENT_LOOP
END SUBROUTINE SCAT_TOTAL_BOUNDARY_SEARCH
! End DBD addition 21 May 2003
!*******************************************************************************


! DBD addition 30 July 2004
SUBROUTINE SCAT_TOTAL_BOUNDARY_SEARCH_SPH
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! This routine finds faces lying on the boundary between the internal inhomegeneous SPHERICAL materials
! (material labels not equal to HOMOG_MEDIUM) and external homogeneous material label HOMOG_MEDIUM
! and flags them appropriately.
! IT ONLY WORKS FOR SPHERES, CENTERED ON THE ORIGIN, AND THIS IS NOT CHECKED!!
! It was written as a work-around due to problems caused by FEMAP generating non-conforming meshes
! at the region interface.
!
! Note that this should be performed BEFORE PML_ASSIGN_PROPERTIES is called!
! (This routine overwrites some of the external homogeneous material label with PML-specific 
! properties. 
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem,inodes
  INTEGER(I4B), DIMENSION(3) :: tempedges, tempnodes
  REAL(SP), DIMENSION(3) :: rad ! radius of nodes. 
  REAL(SP), DIMENSION(3) :: diff   

  ! Initialize
  edges%scat_tot_boundary = .FALSE.
  faces%scat_tot_boundary = .FALSE.

  ELEMENT_LOOP: DO ielem = 1,num_elements
    IF (elements(ielem)%material.NE.HOMOG_MEDIUM) THEN
	  ! This is in the inhomogeneous target region. Test if it lies on the (fictitious) boundary between
	  ! the scattered/total field region, by seeing if all the nodes on a face lie on the same radius. 
      DO iface = 1,4
        tempedges = LOCAL_FACEEDGES(iface)
        tempnodes = GLOBAL_FACENODES(ielem,iface)
		DO inodes = 1,3
          rad(inodes) = SQRT(DOT_PRODUCT(vertices(tempnodes(inodes))%coord,vertices(tempnodes(inodes))%coord))
		  diff(inodes) = ABS(rad(inodes)-sph_radius)
        END DO
        IF (SUM(diff).LT.0.001_SP*sph_radius) THEN
! write(fileout,*) 'Element ',ielem,'  Face  ',iface,'  flagged as scat_tot_boundary'
          edges(elements(ielem)%edges(tempedges))%scat_tot_boundary = .TRUE.
          faces(elements(ielem)%faces(iface))%scat_tot_boundary = .TRUE.
		END IF
	  END DO
    END IF    
  END DO ELEMENT_LOOP
END SUBROUTINE SCAT_TOTAL_BOUNDARY_SEARCH_SPH

! End DBD addition 30 July 2004
!*******************************************************************************


! DBD addition 4 August 04
SUBROUTINE SCAT_TOTAL_BOUNDARY_CHECK
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! This routine checks that for elements lying on the boundary between the internal inhomegeneous materials
! (material labels not equal to HOMOG_MEDIUM) and external homogeneous material label HOMOG_MEDIUM, 
! there are no unconnected faces, which indicates a probably meshing error.
! DBD 4 Aug 04. 
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem
  INTEGER(I4B), DIMENSION(3) :: tempedges 

  ELEMENT_LOOP: DO ielem = 1,num_elements

    IF (elements(ielem)%material.NE.HOMOG_MEDIUM) THEN
	  ! This is in the inhomogeneous target region. Test if it lies on the (fictitious) boundary between
	  ! the scattered/total field region. 
      DO iface = 1,4
        tempedges = LOCAL_FACEEDGES(iface)
        IF(elements(ielem)%connect2elem(iface).EQ.0) THEN
		  WRITE(FILEOUT,'(A,I8,2X,A,I2,2X,A,I2)')'WARNING: Possible meshing error on interface : unconnected internal element ',&
		        ielem, 'face ',iface, & 
		        'material ',elements(ielem)%material
		END IF
	  END DO
    END IF    
  END DO ELEMENT_LOOP
END SUBROUTINE SCAT_TOTAL_BOUNDARY_CHECK
! End DBD addition 4 August 04
!*******************************************************************************



FUNCTION LOCAL_EDGENODES(edgenum)
   USE nrtype
   IMPLICIT NONE
!******************************************************************************
! A function to return the two LOCAL nodes associated with edges 1-6.
! DBD 1997. Revised Jan-Feb 1999.
!******************************************************************************
   INTEGER(I4B), INTENT(IN) :: edgenum
   INTEGER(I4B), DIMENSION(2) :: LOCAL_EDGENODES
   SELECT CASE(edgenum)
   CASE(1) ! e12
     LOCAL_EDGENODES = (/1,2/)
   CASE(2) ! e13
     LOCAL_EDGENODES = (/1,3/) 
   CASE(3) ! e14
     LOCAL_EDGENODES = (/1,4/) 
   CASE(4) ! e23
     LOCAL_EDGENODES = (/2,3/) 
   CASE(5) ! e24
     LOCAL_EDGENODES = (/2,4/) 
   CASE(6) ! e34
     LOCAL_EDGENODES = (/3,4/) 
   CASE DEFAULT
     STOP 'IE: Edge no. out of range in function LOCAL_EDGENODES.'
   END SELECT
      
END FUNCTION LOCAL_EDGENODES
!*******************************************************************************


FUNCTION LOCAL_FACEEDGES(facenum)
   USE nrtype
   IMPLICIT NONE
!******************************************************************************
! A function to return the the three LOCAL edges comprising a face
! using the Savage and Peterson face numbering convention.
! Note that the edges 1-6 correspond to e12,e13,e14,e23,e24 and e34.
! DBD 1997. Revised Jan-Feb 1999.
!******************************************************************************
   INTEGER(I4B), INTENT(IN) :: facenum
   INTEGER(I4B), DIMENSION(3) :: LOCAL_FACEEDGES
  
   SELECT CASE(facenum) 
     ! Work out the corresponding local edges.
     CASE(1) ! Face 1; local nodes 1;2;3 -> edges 1,2,4
       LOCAL_FACEEDGES= (/1,2,4/)
     CASE(2) ! Face 2; local nodes 1;2;4 -> edges 1,3,5
       LOCAL_FACEEDGES= (/1,3,5/)
     CASE(3) ! Face 3; local nodes 1;3;4 -> edges 2,3,6
       LOCAL_FACEEDGES= (/2,3,6/)
     CASE(4) ! Face 4; local nodes 2;3;4 -> edges 4,5,6
       LOCAL_FACEEDGES= (/4,5,6/)
     CASE DEFAULT
       STOP 'IE: Local facenum. out of range in function LOCAL_FACEEDGES.'
   END SELECT

END FUNCTION LOCAL_FACEEDGES
!******************************************************************************


FUNCTION LOCAL_FACENODES(facenum)
   USE nrtype
   IMPLICIT NONE
!******************************************************************************
! A function to return the three LOCAL nodes associated with faces 1-4, using
! the convention (after Savage and Peterson): (NOTE: INCOMPATIBLE with
! Lee and Mittra's convention!)
! Face Node
!   1  1 2 3
!   2  1 2 4
!   3  1 3 4
!   4  2 3 4
! DBD Jan-Feb 1999.
!******************************************************************************
   INTEGER(I4B), INTENT(IN) :: facenum
   INTEGER(I4B), DIMENSION(3) :: LOCAL_FACENODES
  
   SELECT CASE(facenum)
   CASE(1)
     LOCAL_FACENODES = (/1,2,3/)
   CASE(2)
     LOCAL_FACENODES = (/1,2,4/)
   CASE(3)
     LOCAL_FACENODES = (/1,3,4/)
   CASE(4)
     LOCAL_FACENODES = (/2,3,4/)
   CASE DEFAULT
     STOP 'IE: Local facenum. out of range in function LOCAL_FACENODES.'
   END SELECT
 
END FUNCTION LOCAL_FACENODES
!******************************************************************************


FUNCTION LOCAL_TRI_EDGENODES(edgenum,facenum)
   USE nrtype
   IMPLICIT NONE
!******************************************************************************
! A function to return the two LOCAL nodes associated with local edges i-iii on 
! the compatible 2D triangular local face of the tet, numbered using the S&P
! convention.
! DBD March 2000.
!******************************************************************************
   INTEGER(I4B), INTENT(IN) :: edgenum,facenum
   INTEGER(I4B), DIMENSION(2) :: LOCAL_TRI_EDGENODES
   SELECT CASE(facenum)
   CASE(1) ! Face 1, local nodes 1,2,3
     SELECT CASE(edgenum)
     CASE(1) ! e12
       LOCAL_TRI_EDGENODES = (/1,2/)
     CASE(2) ! e13
       LOCAL_TRI_EDGENODES = (/1,3/) 
     CASE(3) ! e23
       LOCAL_TRI_EDGENODES = (/2,3/) 
     CASE DEFAULT
       STOP 'IE: Local edgenum. out of range in function LOCAL_TRI_EDGENODES.'
     END SELECT
   CASE(2) ! Face 2, local nodes 1,2,4
     SELECT CASE(edgenum)
     CASE(1) ! e12
       LOCAL_TRI_EDGENODES = (/1,2/)
     CASE(2) ! e14
       LOCAL_TRI_EDGENODES = (/1,4/) 
     CASE(3) ! e24
       LOCAL_TRI_EDGENODES = (/2,4/) 
     CASE DEFAULT
       STOP 'IE: Local edgenum. out of range in function LOCAL_TRI_EDGENODES.'
     END SELECT
   CASE(3) ! Face 3, local nodes 1,3,4
     SELECT CASE(edgenum)
     CASE(1) ! e13
       LOCAL_TRI_EDGENODES = (/1,3/)
     CASE(2) ! e14
       LOCAL_TRI_EDGENODES = (/1,4/) 
     CASE(3) ! e34
       LOCAL_TRI_EDGENODES = (/3,4/) 
     CASE DEFAULT
       STOP 'IE: Local edgenum. out of range in function LOCAL_TRI_EDGENODES.'
     END SELECT
   CASE(4) ! Face 4, local nodes 2,3,4
     SELECT CASE(edgenum)
     CASE(1) ! e23
       LOCAL_TRI_EDGENODES = (/2,3/)
     CASE(2) ! e24
       LOCAL_TRI_EDGENODES = (/2,4/) 
     CASE(3) ! e34
       LOCAL_TRI_EDGENODES = (/3,4/) 
     CASE DEFAULT
       STOP 'IE: Local edgenum. out of range in function LOCAL_TRI_EDGENODES.'
     END SELECT
   CASE DEFAULT
       STOP 'IE: Local facenum. out of range in function LOCAL_TRI_EDGENODES.'
   END SELECT
      
END FUNCTION LOCAL_TRI_EDGENODES
!******************************************************************************


FUNCTION LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face,local)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This function returns the corresponding local tetrahedral dof number, of the
! local, facial dof number <local>, associated with local face <local_face>.
! Max range of local = 1:12 => QT/QN element.
! 2 March 2002 DBD.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: local_face,local
  INTEGER(I4B) :: LOCAL_TRI_TO_LOCAL_TET_INDEX

  INTEGER(I4B), DIMENSION(3) :: faceedges

  ! Find the local edges of the face:
  faceedges = LOCAL_FACEEDGES(local_face)

  SELECT CASE(local)
  CASE(1:3)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = faceedges(local)
  CASE (4:6)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = 6 + faceedges(local-3)
  CASE (7)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = 12 + local_face
  CASE (8)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = 16 + local_face
  CASE (9:11)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = 20 + faceedges(local-8)
  CASE (12)
    LOCAL_TRI_TO_LOCAL_TET_INDEX = 26 + local_face
  CASE DEFAULT
    STOP 'IE: out-of-range local value in LOCAL_TRI_TO_LOCAL_TET_INDEX.'
  END SELECT

END FUNCTION LOCAL_TRI_TO_LOCAL_TET_INDEX
!*******************************************************************************


FUNCTION MAX_ORDER(elem,local_face)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Determine the maximum hierarchal order of any entity associated with the
! requested element <elem>. Optionally the maximum hierarchal order of any
! entity associated with local face <face> is returned instead.
! Created MMB. 2001-10-11.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem
  INTEGER(I4B), INTENT(IN), OPTIONAL :: local_face
  INTEGER(I4B) :: MAX_ORDER
  INTEGER(I4B), DIMENSION(3) :: face_edges

  IF (.NOT.PRESENT(local_face)) THEN
    MAX_ORDER = MAX(elements(elem)%order,                 &
                    faces(elements(elem)%faces(1))%order, &
                    faces(elements(elem)%faces(2))%order, &
                    faces(elements(elem)%faces(3))%order, &
                    faces(elements(elem)%faces(4))%order, &
                    edges(elements(elem)%edges(1))%order, &
                    edges(elements(elem)%edges(2))%order, &
                    edges(elements(elem)%edges(3))%order, &
                    edges(elements(elem)%edges(4))%order, &
                    edges(elements(elem)%edges(5))%order, &
                    edges(elements(elem)%edges(6))%order)
  ELSE
    face_edges = LOCAL_FACEEDGES(local_face)
    MAX_ORDER = MAX(faces(elements(elem)%faces(local_face))%order,    &
                    edges(elements(elem)%edges(face_edges(1)))%order, &
                    edges(elements(elem)%edges(face_edges(2)))%order, &
                    edges(elements(elem)%edges(face_edges(3)))%order)
  END IF

END FUNCTION MAX_ORDER
!*******************************************************************************

FUNCTION MIXED_ORDER(elem,local_face)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Determine the  mixed/complete order of any entity associated with the
! requested element <elem>. Note the complete order overrides mixed
! order. Optionally the mixed/complete order of any entity associated with 
! local face <face> is returned instead.
! Returns <.true.> if the maximum order entity(ies) are of mixed order, and 
! <.false.> if there is a maximum order entity of complete order.
!
! 2002-03-04: Created. DBD.
! 2002-03-17: Corrections - an order 2-mixed entity overrides an order 
!             1-complete entity. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem
  INTEGER(I4B), INTENT(IN), OPTIONAL :: local_face
  LOGICAL(LGT) :: MIXED_ORDER

  INTEGER(I4B), DIMENSION(3) :: face_edges
  INTEGER(I4B) :: max_order_entity,icount
  
  ! Initialise:
  MIXED_ORDER = .TRUE.

  IF (.NOT.PRESENT(local_face)) THEN

    max_order_entity = MAX_ORDER(elem)
	
	! Check the element's mixed status:
	IF (elements(elem)%order.EQ.max_order_entity) THEN
	  MIXED_ORDER = MIXED_ORDER.AND.elements(elem)%mixed
	END IF

	DO icount = 1,4 ! cycle through faces
	  IF (faces(elements(elem)%faces(icount))%order.EQ.max_order_entity) THEN
	    MIXED_ORDER = MIXED_ORDER.AND.faces(elements(elem)%faces(icount))%mixed
	  END IF
	END DO
	
	DO icount = 1,6 ! cycle through edges
	  IF (edges(elements(elem)%edges(icount))%order.EQ.max_order_entity) THEN
	    MIXED_ORDER = MIXED_ORDER.AND.edges(elements(elem)%edges(icount))%mixed
	  END IF
	END DO

  ELSE

    face_edges = LOCAL_FACEEDGES(local_face)
    max_order_entity = MAX_ORDER(elem,local_face)

	! Check the element face's mixed status:
	IF (faces(elements(elem)%faces(local_face))%order.EQ.max_order_entity) THEN
	  MIXED_ORDER = MIXED_ORDER.AND.faces(elements(elem)%faces(local_face))%mixed
	END IF

	DO icount = 1,3 ! cycle through edges
	  IF (edges(elements(elem)%edges(face_edges(icount)))%order.EQ.max_order_entity) THEN
	    MIXED_ORDER = MIXED_ORDER.AND.edges(elements(elem)%edges(face_edges(icount)))%mixed
	  END IF
	END DO

  END IF

END FUNCTION MIXED_ORDER
!*******************************************************************************

SUBROUTINE MESHSORT 
   USE nrtype
   USE problem_info
   USE unit_numbers
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine handles mesh sorting. The nodes associated with an 
! element are sorted into increasing order.
! The sort algorithm used is Insertion Sort; see "Introduction to 
! Algorithms", Cormen et al., MIT Press, 1996, pp2-4.
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! D B Davidson
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! 12 March 1997 by DBD. 
! Minor revision 28 April.
! Minor revision 16 June.
! Revised 24 Jan 1999 by DBD - nodes numbered from 1 now.
! Corrected 28 Jan 1999.
! Documentation corrected 27 Apr 1999 by DBD.
! Updated to use new data structures 7 Mar 2000 by DBD.
! Minor output format change 10 May 2000 by DBD.
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the data structure "elements".
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is "elements" with the list of nodes per element "elements%nodes"
! sorted into ascending order. This is done so that the edge directions are 
! consistent for elements sharing an edge throughout the code.
!*******************************************************************************
    INTEGER(I4B) i, jj, k, key, inode, jnode


    INTEGER(I4B), DIMENSION(4) :: A     ! Note that this must have the same
                                 ! dimension as the numbers of vertices
                                 ! Four for a tetrahedron of order one-half.
    CHARACTER(8) :: elem_type    ! string to output the type of elements (CT/LN etc.)

    DO k=1,num_elements ! Start sort loop

       ! Set up temporary data structure.
       A(1:4)  = elements(k)%nodes(1:4)
  
       ! Now perform insertion sort
       DO jj=2,SIZE(A)
          key = A(jj)  
          ! Insert A(jj) into the sorted sequence A(1..jj-1)
          i = jj-1
          DO WHILE (i.GT.0.AND.A(i).GT.key)
            A(i+1) = A(i)
            i = i-1 
          END DO 
          A(i+1) = key
       END DO ! Sorting finished.

       ! Update elements
       elements(k)%nodes(1:4) = A(1:4)
      
    END DO ! Element loop

    ! LIST OF ELEMENTS
    WRITE(FILEOUT,'(//,20X,A)') 'LIST OF ELEMENTS'
    WRITE(FILEOUT,'(/,7(A14))') 'Element no.','Type', 'Material no.', 'Node 1', 'Node 2', 'Node 3', 'Node 4'
    SELECT CASE (HIERARCHAL_ORDER)
    CASE (1) ! Changed DBD 31 Aug 2004.
	  IF (MIXED_ORDER_FLAG) THEN
        elem_type = '   CT/LN'
	  ELSE
        elem_type = '   LT/LN'
	  END IF 
    CASE (2)
	  IF (MIXED_ORDER_FLAG) THEN
        elem_type = '   LT/QN'
	  ELSE
        elem_type = '   QT/QN'
	  END IF 
          END SELECT ! End changed DBD 31 Aug 2004.
    IF (OUTPUT_ELEMENT_DATA) THEN
      DO i=1,num_elements
         WRITE(FILEOUT,'(2X,I10,2X,4X,A10,5(2X,I10,2X))') i, elem_type, elements(i)%material, &
                                                             elements(i)%nodes(1:4)
      END DO
    ELSE 
      WRITE(FILEOUT,'(1X,A)') 'Suppressed by option OUTPUT_ELEMENT_DATA'
    END IF

  RETURN

END SUBROUTINE MESHSORT
!*******************************************************************************


SUBROUTINE NODEELEMENT_INDEXLIST_CLEAN
  IMPLICIT NONE
!******************************************************************************
! Purpose: Clean up arrays of node-element index list (allocated in subroutine
! MAKE_NODEELEMENT_INDEXLIST).
! Author:  FJCM
!******************************************************************************

  IF (ALLOCATED(Node_ind)) DEALLOCATE(Node_ind)    
  IF (ALLOCATED(Node_ptr)) DEALLOCATE(Node_ptr)

END SUBROUTINE NODEELEMENT_INDEXLIST_CLEAN
!******************************************************************************


SUBROUTINE NODEELEMENT_INDEXLIST_MAKE
   USE nrtype
   IMPLICIT NONE   
!*******************************************************************************
! This subroutine builds a node-element index list, indicating which elements 
! is connected to every node. It must be called before calling fast_faceconnect,
! fast_edgeconnect, fast_edgemake, fast_facemake.  
! Output is:
! * An array Node_ptr(:) pointing to the start of each node in the node-element list
! * An array Node_nind(:) containing the elements associated with each node starting
!   at node one and moving upwards. This is the node-element list.  
!
! Original version 2 April 2001 - FJC Meyer. Routine adapted
!                                 from Metis mesh2graph routine.
!*******************************************************************************
   INTEGER(I4B) ielem,jelem,kface,lface ! counters
   INTEGER(I4B) i,k, knodes,indx_j      ! more counters    
   INTEGER(I4B), DIMENSION(1:3) :: lfacenodes, kfacenodes
   
   ! First check to see if the data structures are already created:
   IF (ALLOCATED(Node_ind).AND.ALLOCATED(Node_ptr)) RETURN

   ! Create a node-element list
   ALLOCATE(Node_ptr(0:num_nodes+1))     
   Node_ptr(:) = 0
   DO i=1,num_elements         ! Loops over all elements, counting how many times
     DO k=1,4                  ! each node is referenced. (Comment NM 2005.04)
       Node_ptr(elements(i)%nodes(k)) = Node_ptr(elements(i)%nodes(k)) + 1
     END DO 
   END DO

   ! CSR format from long array
   DO i=2, num_nodes            ! Generate a cumulative sum of 
     Node_ptr(i) = Node_ptr(i) + Node_ptr(i-1) ! refcounts. ie. Node_ptr(i)
   END DO                       ! = total no. of refs from node 1 through i
   DO i=num_nodes+1, 2, -1   
     Node_ptr(i) = Node_ptr(i-1) 
   END DO
   Node_ptr(1) = 0 

   ! Create a node indexing array
   ALLOCATE(Node_ind(Node_ptr(num_nodes+1)))     
   DO i=1, num_elements
     DO k=1,4
       ! Add one for we have to move one up in the index each time the same
       ! node is worked with
       Node_ptr(elements(i)%nodes(k)) = Node_ptr(elements(i)%nodes(k)) + 1     
       Node_ind(Node_ptr(elements(i)%nodes(k))) = i
     END DO
   END DO
   
   ! Interesting. This restores the original nptr values.  
   DO i=num_nodes+1, 2, -1   
     Node_ptr(i) = Node_ptr(i-1)
   END DO
   Node_ptr(1) = 0

END SUBROUTINE NODEELEMENT_INDEXLIST_MAKE
!******************************************************************************


SUBROUTINE PEC_BOUNDARY_SEARCH
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! Flags the boundary edges and faces as PEC if they are not assigned as 
! CBAA/coax/port.
! Note: this routine must only be called after these other flags have been 
! assigned.
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem
  INTEGER(I4B), DIMENSION(3) :: tempedges 

  ELEMENT_LOOP: DO ielem = 1,num_elements
    FACE_LOOP: DO iface = 1,4
      FACE_CONNECT_TEST: IF (elements(ielem)%connect2elem(iface).EQ.0) THEN     
        tempedges = LOCAL_FACEEDGES(iface)

        ! The rest of the faces are PEC in the cases of CBAA, EIG and GW:
        IF ((.NOT.faces(elements(ielem)%faces(iface))%CBAA_aperture).AND. &
            (.NOT.faces(elements(ielem)%faces(iface))%coax_aperture).AND. &
            (.NOT.faces(elements(ielem)%faces(iface))%port)) THEN
          edges(elements(ielem)%edges(tempedges))%PEC = .TRUE.
          faces(elements(ielem)%faces(iface))%PEC = .TRUE.
        END IF

      END IF FACE_CONNECT_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE PEC_BOUNDARY_SEARCH
!*******************************************************************************


SUBROUTINE PEC_OBJECT_SEARCH
  USE nrtype
  IMPLICIT NONE   
!*******************************************************************************
! Search for edges/faces belonging to internal/surface plates (BC 2 card) and set
! the appropriate flags.
!*******************************************************************************
  INTEGER(I4B) :: icard,row,iedge,iface,ielem
  LOGICAL(LGT) :: in_PEC_plane
  INTEGER(I4B), DIMENSION(3) :: tempedges 

  CARD_LOOP: DO icard = 1,num_BCs ! Cycle through BC 2 cards - quadrilateral PEC surfaces
    IF (BCs(icard)%type.NE.2) CYCLE CARD_LOOP ! not a PEC boundary condition
    ELEMENT_LOOP: DO ielem = 1,num_elements
      FACE_LOOP: DO iface = 1,4 ! Face-wise search to find edges in the quad.
        CALL FACE_IN_QUADRILATERAL(BCs(icard)%corner(1,1:3),BCs(icard)%corner(2,1:3),  &
                                   BCs(icard)%corner(3,1:3),BCs(icard)%corner(4,1:3),  &
                                   elements(ielem)%faces(iface),in_PEC_plane)
        IF (in_PEC_plane) THEN
          tempedges = LOCAL_FACEEDGES(iface)
          DO iedge = 1,3 ! Flag the three edges as PEC
            edges(elements(ielem)%edges(tempedges(iedge)))%PEC = .TRUE.
          END DO
          faces(elements(ielem)%faces(iface))%PEC = .TRUE.
        END IF
      END DO FACE_LOOP
    END DO ELEMENT_LOOP
  END DO CARD_LOOP

END SUBROUTINE PEC_OBJECT_SEARCH
!*******************************************************************************

SUBROUTINE PEC_SCATTERER_SEARCH
  USE nrtype
  use scattering_analysis_data
  IMPLICIT NONE   
!*******************************************************************************
! Search for edges/faces belonging to PEC scatterers and set
! the appropriate flags.  Current implementation is very specific to 
! scattering analysis; it is assumed that the region internal to the homogoneous
! outer is PEC. 
! NB! Routine SCAT_TOTAL_BOUNDARY_SEARCH must be called BEFORE this one,
! for the scattered field formulation.
!
! Note also a subtle difference between the total field and scattered field 
! formulations; in the former, the field is computed EXTERNAL to the PEC; 
! in the latter, the field must be computed EVERYWHERE. 
!
! Written DBD 28 Jul 04.
! Extended and corrected DBD 2 Aug 04.  
! For FD scattering analysis.
!*******************************************************************************
  INTEGER(I4B) :: iedge,iface,ielem
  INTEGER(I4B), DIMENSION(3) :: tempedges 

  IF (SCAT_FIELD) THEN
    ELEMENT_LOOP1: DO ielem = 1,num_elements
      ! Flag all edges and faces on the scattered/total boundary as Dirichlet. 
      FACE_LOOP1: DO iface = 1,4 
	    IF (faces(elements(ielem)%faces(iface))%scat_tot_boundary) THEN
          faces(elements(ielem)%faces(iface))%Dirichlet = .TRUE. ! Ditto face.
	      tempedges = LOCAL_FACEEDGES(iface)
          DO iedge = 1,3 
            edges(elements(ielem)%edges(tempedges(iedge)))%Dirichlet = .TRUE.
          END DO
  		END IF
      END DO FACE_LOOP1	  
    END DO ELEMENT_LOOP1
  ELSE ! Total field formulation.
    ELEMENT_LOOP2: DO ielem = 1,num_elements
      IF (elements(ielem)%material.NE.HOMOG_MEDIUM) THEN
        ! Set all edges and faces corresponding to the scatterer to PEC, INCLUDING those on the boundary of the scatterer. 
        FACE_LOOP2: DO iface = 1,4 
          faces(elements(ielem)%faces(iface))%PEC = .TRUE. ! Ditto face.
          tempedges = LOCAL_FACEEDGES(iface)
          DO iedge = 1,3 ! Flag the three edges as PEC.
            edges(elements(ielem)%edges(tempedges(iedge)))%PEC = .TRUE.
          END DO
  	    END DO FACE_LOOP2
	  END IF
    END DO ELEMENT_LOOP2
  END IF 

END SUBROUTINE PEC_SCATTERER_SEARCH
!*******************************************************************************


SUBROUTINE PML_ASSIGN_PROPERTIES 
  USE boundary_conditions
  USE nrtype
  USE material_properties
  IMPLICIT NONE   
!*******************************************************************************
! Assigns PML propeties to appropriate reqions in mesh. 
! There are seven different zones which are labelled:
! 
! Absorbing in     Label (offset)
! +/- directions   
!  
! x                1 
! y                2
! z                3
! x&y              4
! x&z              5
! y&z              6
! x&y&z            7
!
!*******************************************************************************
  INTEGER(I4B) :: ielem,ii
  REAL(SP), DIMENSION(3) :: elem_centre
  LOGICAL (LGT) S1_found, S2_found
  REAL(SP) d, dist_x, dist_y, dist_z
  INTEGER(I4B) zone_label
  
  ! Logic similar to SUBROUTINE TD_INC_FIELD_SETUP 
  S1_found = .false.
  S2_found = .false.

  ! First, assign opposite points corresponding to S1 and S2.
  DO ii = 1,num_DPoints
    IF(ABC_Box%S1.EQ.DPoints(ii)%name) THEN
      PML%x_min = DPoints(ii)%coords(1)
      PML%y_min = DPoints(ii)%coords(2)
      PML%z_min = DPoints(ii)%coords(3)
      S1_found = .true.
	ELSE IF (ABC_Box%S2.EQ.DPoints(ii)%name) THEN
      PML%x_plus = DPoints(ii)%coords(1)
      PML%y_plus = DPoints(ii)%coords(2)
      PML%z_plus = DPoints(ii)%coords(3)
	  S2_found = .true.
    END IF
  END DO

  ! Check that S1 and S2 were found.
  IF(.NOT.S1_found.OR..NOT.S2_found) THEN
    CALL ERROR_FEMFEKO(1,4091)
  END IF
  
  ELEMENT_LOOP: DO ielem = 1,num_elements
    zone_label = 0 ! Default - not a PML region. 
	d = PML%thickness
    elem_centre = ELEMENT_CENTRE(ielem)
	dist_x = MIN(ABS(elem_centre(1)-PML%x_min ),ABS(elem_centre(1)-PML%x_plus))
    dist_y = MIN(ABS(elem_centre(2)-PML%y_min ),ABS(elem_centre(2)-PML%y_plus))
    dist_z = MIN(ABS(elem_centre(3)-PML%z_min ),ABS(elem_centre(3)-PML%z_plus))
!    if (ielem.eq.33) then
!	continue; ! for debugging
!	end if
    IF(dist_x.LE.d.AND.dist_y.GT.d.AND.dist_z.GT.d) THEN 
	   zone_label = 1
	END IF
    IF(dist_x.GT.d.AND.dist_y.LE.d.AND.dist_z.GT.d) THEN 
	   zone_label = 2
	END IF
    IF(dist_x.GT.d.AND.dist_y.GT.d.AND.dist_z.LE.d) THEN 
       zone_label = 3
	END IF
    IF(dist_x.LE.d.AND.dist_y.LE.d.AND.dist_z.GT.d) THEN 
       zone_label = 4
	END IF
    IF(dist_x.LE.d.AND.dist_y.GT.d.AND.dist_z.LE.d) THEN 
       zone_label = 5
	END IF
    IF(dist_x.GT.d.AND.dist_y.LE.d.AND.dist_z.LE.d) THEN 
       zone_label = 6
	END IF
    IF(dist_x.LE.d.AND.dist_y.LE.d.AND.dist_z.LE.d) THEN 
       zone_label = 7
	END IF
	IF(zone_label.NE.0) THEN
      elements(ielem)%material = MAX_MATERIALS + zone_label
      material_type(elements(ielem)%material) = 3 ! Label as PML.
	ELSE
	  CONTINUE ! No action, this element is not in the PML region. 
	END IF
    ! The above labels are now used when the element in the PML is accessed. 

  END DO ELEMENT_LOOP


END SUBROUTINE PML_ASSIGN_PROPERTIES
!*******************************************************************************



FUNCTION SIGNED_VOLUME(i)
  USE problem_info
  USE math_tools, ONLY: DET_DIM4
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Compute signed volume of element i.
! Otherwise this function is identical to FUNCTION VOLUME.
!*******************************************************************************
  REAL(SP) SIGNED_VOLUME
  INTEGER(I4B), INTENT(IN)  :: i
  REAL(SP), DIMENSION(4,4) :: temp

  !   Set up temporary matrix  
  temp(:,1) = 1
  temp(1,2:4) = vertices(elements(i)%nodes(1))%coord
  temp(2,2:4) = vertices(elements(i)%nodes(2))%coord
  temp(3,2:4) = vertices(elements(i)%nodes(3))%coord
  temp(4,2:4) = vertices(elements(i)%nodes(4))%coord

  SIGNED_VOLUME = DET_DIM4(temp)/6.0_SP  

  IF(DEBUG_VOLUME) THEN ! Output debugging info
    WRITE(FILEOUT,'(//A,A)') 'Element volumes - signed'
    WRITE(FILEOUT,'(A,I8,A,G10.4)') 'Element number: ',i,'Volume: ', SIGNED_VOLUME
  END IF

END FUNCTION SIGNED_VOLUME
!*******************************************************************************


SUBROUTINE SIMPLEX_COEFFICIENTS(element_num,a,b,c,d,coord_mat)    
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the coefficients a_i,b_i,c_i and d_i:
! \lambda = a_i +b_i x + c_i y +d_i z  (i=1,..,4)
! 
! See eqns. (4) & (6), Savage and Peterson, IEEE T-AP, June 96.
!
! Special warning - this routine uses Fortran 77 LAPACK 
! SINGLE PRECISION routines SGETRF and SGETRI; this must be changed MANUALLY if 
! different precision is required! Note also: be careful on especially PC ports 
! - default integers lengths vary, many PC compiler use 2 bytes. 
! I4B is a four-byte integer. Integer lengths must be consistent between
! this routine and (possibly pre-compiled) LAPACK routines. 
! They may be downloaded from:
! http://www.netlib.org/lapack/single
!
! Initial version 16 March 2000. Based on tested code moved from routine
! S_AND_AND_T_MAKE_HIERARCHAL. DBD
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: element_num
  REAL(SP), DIMENSION(4), INTENT(OUT), OPTIONAL :: a,b,c,d  
                                      ! gradient terms: See [eqn(6),S&P] 
  REAL(SP), DIMENSION(4,4), INTENT(OUT), OPTIONAL :: coord_mat
                                      ! coordinate matrix augmentened by ones
  REAL(SP), DIMENSION(4,4) :: temp_mat
  INTEGER(I4B), DIMENSION(4) :: el_nodes
  INTEGER(I4B) :: info               ! Info from LAPACK factorization.
  INTEGER(I4B) :: lwork              ! LAPACK dimensioning data
  PARAMETER(lwork=4)
  REAL(SP), DIMENSION(lwork) :: temp_work ! LAPACK temporary workspace.
  INTEGER(I4B), DIMENSION(4) :: temp_ipiv ! Pivot indices from factorization. 


 ! Build the matrix on RHS of [eqn(6),S&P] 
  el_nodes(1:4) = elements(element_num)%nodes(1:4)
  temp_mat(1,1:4) = vertices(el_nodes(1:4))%coord(1) ! x nodal co-ordinates
  temp_mat(2,1:4) = vertices(el_nodes(1:4))%coord(2) ! y nodal co-ordinates
  temp_mat(3,1:4) = vertices(el_nodes(1:4))%coord(3) ! z nodal co-ordinates
  temp_mat(4,1:4) = SPREAD(1.0_SP,1,4)        ! ones  

  IF (PRESENT(coord_mat)) THEN
    coord_mat = temp_mat
  END IF

  IF ((.NOT.PRESENT(a)).OR.(.NOT.PRESENT(b)).OR. &
      (.NOT.PRESENT(c)).OR.(.NOT.PRESENT(d))) RETURN

  ! Now invert matrix. Uses LAPACK routine SGETRF to factorize, and then 
  ! SGETRI to invert. Note - THESE ARE SINGLE PRECISION ROUTINES!! 
  ! See note above!!
  CALL SGETRF ((4),(4),temp_mat,(4),temp_ipiv,info)  
  ! Check error conditions - shouldn't occur.
  IF(Info.NE.0) THEN 
    CALL ERROR_FEMFEKO(1,4034)
  END IF
  CALL SGETRI ((4),temp_mat,(4),temp_ipiv,temp_work,lwork,info)
  ! Check error conditions - also shouldn't occur.
  IF(Info.NE.0) THEN 
    CALL ERROR_FEMFEKO(1,4035)
  END IF
  ! Note that the matrix has now been overwritten with its inverse

  b(1:4) = temp_mat(1:4,1)     
  c(1:4) = temp_mat(1:4,2) 
  d(1:4) = temp_mat(1:4,3) 
  a(1:4) = temp_mat(1:4,4)

END SUBROUTINE SIMPLEX_COEFFICIENTS
!*******************************************************************************


FUNCTION SIMPLEX_COORDINATES(element_num,x,y,z)
   USE problem_info
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE unit_numbers
   IMPLICIT NONE
!*******************************************************************************
! A function that returns the simplex coordinates at point x,y,z for element
! element_num. Note that x,y, and z need not lie within the element, in 
! which case the simplex coordinates may become negative and individually exceed
! one. DBD June 99, revised 17 Feb 2000.
!*******************************************************************************
! ERROR DIRECTORY
!*******************************************************************************
! 4950: IE: internal error following call to SGETRF.
! 4951: IE: internal error following call to SGETRI. 
! 4952: IE: internal error computing simplex coordinates. 
!*******************************************************************************
   INTEGER(I4B), INTENT(IN)::  element_num         ! Element number
   REAL(SP), INTENT(IN) :: x, y, z                 ! Coordinates of point
   REAL(SP), DIMENSION(4) :: SIMPLEX_COORDINATES   ! Simplex coordinates of
                                                   !  point.
   INTEGER(I4B) inode                              ! Node number
   REAL(SP), DIMENSION(4):: global_coords          ! See eqns (4&6) in ref. 
                                                   ! below.
   INTEGER(I4B), DIMENSION(4) :: el_nodes          ! Nodes of element
   REAL(SP), DIMENSION(4,4) :: coord_mat           ! Matrix in eqn.(6) of ref.
                                                   ! see below
   INTEGER(I4B) :: info               ! Info from LAPACK factorization.
   INTEGER(I4B) :: lwork              ! LAPACK dimensioning data
   PARAMETER(lwork=4)
   REAL(SP), DIMENSION(lwork) :: temp_work ! LAPACK temporary workspace.
   INTEGER(I4B), DIMENSION(4) :: temp_ipiv ! Pivot indices from factorization. 
   LOGICAL(LGT), SAVE :: first_pass  =.true.       ! Flag (initialized)

   ! Find the simplex coordinates of point (x,y,z). 
   ! See [eqn(6),Savage and Peterson, IEEE T-AP June 1996, pp.874--879] 
   el_nodes(1:4) = elements(element_num)%nodes(1:4)
   coord_mat(1,1:4) = vertices(el_nodes(1:4))%coord(1) ! x nodal co-ordinates
   coord_mat(2,1:4) = vertices(el_nodes(1:4))%coord(2) ! y nodal co-ordinates
   coord_mat(3,1:4) = vertices(el_nodes(1:4))%coord(3) ! z nodal co-ordinates
   coord_mat(4,1:4) = SPREAD(1.0_SP,1,4)        ! ones  
   global_coords(1)  = x
   global_coords(2)  = y
   global_coords(3)  = z
   global_coords(4)  = 1.0_SP

   ! Now invert matrix. Uses LAPACK routine SGETRF to factorize, and then 
   ! SGETRI to invert. Note - THESE ARE SINGLE PRECISION ROUTINES!! 
   CALL SGETRF ((4),(4),coord_mat,(4),temp_ipiv,info)  
   IF(Info.NE.0)  CALL ERROR_FEMFEKO(1,4100,int1=info)
   CALL SGETRI ((4),coord_mat,(4),temp_ipiv,temp_work,lwork,info)
   IF(Info.NE.0)  CALL ERROR_FEMFEKO(1,4101,int1=info)

   ! Note that the matrix has now been overwritten with its inverse
   SIMPLEX_COORDINATES = MATMUL(coord_mat,global_coords)
   ! Note that these are "unnormalized" in S&P notation but consistent
   ! with usage in rest of code, since they effectively INCLUDE the volume.
   
   ! Do some error checking. Simplex coordinates should sum to 1.
   IF ((ABS(SUM(SIMPLEX_COORDINATES))-1.0).GT.EPS) THEN
     CALL ERROR_FEMFEKO(1,4102)
   END IF
   
   IF (DEBUG_ELEMENT_E_FIELD) THEN
     IF (first_pass) THEN
       WRITE(FILEOUT,'(//A,A/)') '********** SIMPLEX ',&    
                     'COORDINATES  ****************'
       first_pass = .false.
     END IF
     WRITE(FILEOUT,'(/A,I4,3(A,F6.3))') 'element= ',element_num,& 
       ' x=',x,' y=',y,' z=',z
     WRITE(FILEOUT,*) '  lambda_1', '  lambda_2','  lambda_3', '  lambda_4'
     WRITE(FILEOUT,'(4(F9.5,2X))') SIMPLEX_COORDINATES(1:4)
   END IF

END FUNCTION SIMPLEX_COORDINATES
!*******************************************************************************


FUNCTION TET2TRI_EDGENUMS(iedge,jface)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This function converts an edge numbered in 
! the 3D tetraheron-based edge numbering scheme 
! to one numbered in the 2D triangle-based edge numbering scheme.
!
! The conversion scheme, using the S&P edge and face numbering scheme, is:
! (this is fixed once the edge and face numbering scheme is chosen):
! 
! Tet face  Tet Nodes Tet edges   Tri edges
!    1       1 2 3      1 2 4       1 2 3 
!    2       1 2 4      1 3 5       1 2 3 
!    3       1 3 4      2 3 6       1 2 3 
!    4       2 3 4      4 5 6       1 2 3 
!
! MMB - moved from gw_sys.f90 and made external 6 Feb 2001 (routine by DBD).
!*******************************************************************************
  INTEGER (I4B), INTENT(IN):: iedge,jface
  INTEGER (I4B) :: TET2TRI_EDGENUMS,iedge2D
  SELECT CASE(jface)
  CASE(1) ! Nodes 1,2,3 
    SELECT CASE(iedge) ! 3D: edges 1,2 and 4
    CASE(1:2) 
      iedge2D = iedge 
    CASE(4) 
      iedge2D = iedge-1
    END SELECT
  CASE(2) ! Nodes 1,2,4
    SELECT CASE(iedge) ! 3D: edges 1,3 and 5
    CASE(1) 
      iedge2D = iedge 
    CASE(3) 
      iedge2D = iedge-1
    CASE(5) 
      iedge2D = iedge-2
    END SELECT
  CASE(3) ! Nodes 1,3,4
    SELECT CASE(iedge) ! 3D: edges 2,3 and 6
    CASE(2:3) 
      iedge2D = iedge-1
    CASE(6) 
      iedge2D = iedge-3
    END SELECT
  CASE(4) ! Nodes 2,3,4
    SELECT CASE(iedge) ! 3D: edges 4,5 and 6
    CASE(4:6) 
      iedge2D = iedge-3
    END SELECT
  CASE DEFAULT
    STOP 'IE: In TET2TRI_EDGENUMS, out-of-range face number.'
  END SELECT
  TET2TRI_EDGENUMS = iedge2D
END FUNCTION TET2TRI_EDGENUMS
!*******************************************************************************


FUNCTION T_LENGTH(global_first_node,global_second_node)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! The length of the vector from the (global) first node 1 to global 
! second node.
!*******************************************************************************
   REAL(SP) t_length
   REAL(SP), DIMENSION(3) :: temp
   INTEGER(I4B), INTENT(IN) :: global_first_node,global_second_node
   temp = vertices(global_second_node)%coord  - & 
          vertices(global_first_node)%coord
   t_length = SQRT(DOT_PRODUCT(temp,temp))
END FUNCTION T_LENGTH
!*******************************************************************************


FUNCTION T_LENGTHS(element_num)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! Find ALL the lengths of edges of element number element_num.
! THIS ROUTINE ASSUMES S&P NUMBERING !
!*******************************************************************************
   REAL(SP), DIMENSION(6) :: t_lengths
   INTEGER(I4B), INTENT(IN) :: element_num
   INTEGER(I4B), DIMENSION(4) :: elnodes
   elnodes(1:4) = elements(element_num)%nodes(1:4)
   t_lengths(1) = t_length(elnodes(1),elnodes(2))
   t_lengths(2) = t_length(elnodes(1),elnodes(3))
   t_lengths(3) = t_length(elnodes(1),elnodes(4))
   t_lengths(4) = t_length(elnodes(2),elnodes(3))
   t_lengths(5) = t_length(elnodes(2),elnodes(4))    
   t_lengths(6) = t_length(elnodes(3),elnodes(4))
END FUNCTION T_LENGTHS   
!*******************************************************************************

FUNCTION VOLUME(ii)
  USE problem_info
  USE math_tools, ONLY: DET_DIM4
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Compute (unsigned) volume of element ii   
!*******************************************************************************
  REAL(SP) VOLUME
  INTEGER(I4B), INTENT(IN)  :: ii
  REAL(SP), DIMENSION(4,4) :: temp

  ! Set up temporary matrix  
  temp(:,1) = 1
  temp(1,2:4) = vertices(elements(ii)%nodes(1))%coord
  temp(2,2:4) = vertices(elements(ii)%nodes(2))%coord
  temp(3,2:4) = vertices(elements(ii)%nodes(3))%coord
  temp(4,2:4) = vertices(elements(ii)%nodes(4))%coord

  VOLUME = ABS(DET_DIM4(temp))/6.0_SP  

  IF(DEBUG_VOLUME) THEN ! Output debugging info
    WRITE(FILEOUT,'(//A,A)') 'Element volumes'
    WRITE(FILEOUT,'(A,I8,A,G10.4)') 'Element number: ',ii,'Volume: ', VOLUME
  END IF
  ! Note that the determinant is 6 times the (potentially signed) "volume".
  ! The expressions in the references are for the actual (unsigned) volumes.

END FUNCTION VOLUME
!*******************************************************************************


FUNCTION XYZ_COORDINATES(element_num,lambda)
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE unit_numbers
   IMPLICIT NONE
!*******************************************************************************
! A function that returns the Cartesian coordinates (x,y,z) of the point
! given by simplex coordinates (lambda_1,lambda_2,lambda_3,lambda_4) 
! in element element_num.  See Lab notebook, FEM VolII, pp.59-60.
! DBD 17 March 2000.
!*******************************************************************************
   INTEGER(I4B), INTENT(IN)::  element_num         ! Element number
   REAL(SP), DIMENSION(4), INTENT(IN) :: lambda    ! Simplex coordinates of
                                                   ! point
   REAL(SP), DIMENSION(3) :: XYZ_COORDINATES       ! Cartesian coordinates of
                                                   ! point.
   REAL(SP), DIMENSION(4,4)  :: vert_mat
   REAL(SP), DIMENSION(4)  :: cart

   ! Set up the vertex matrix: (S&P IEEE T-AP June 1996, eq.6)
   vert_mat(1:3,1) = vertices(elements(element_num)%nodes(1))%coord(1:3)
   vert_mat(1:3,2) = vertices(elements(element_num)%nodes(2))%coord(1:3)
   vert_mat(1:3,3) = vertices(elements(element_num)%nodes(3))%coord(1:3)
   vert_mat(1:3,4) = vertices(elements(element_num)%nodes(4))%coord(1:3)
   vert_mat(4,1:4) = 1.0

   cart = MATMUL(vert_mat,lambda)
  
   IF (ABS(cart(4)-1.0).GT.EPS) THEN
	 CALL ERROR_FEMFEKO(1,4103)
   END IF

   XYZ_COORDINATES = cart(1:3)

END FUNCTION XYZ_COORDINATES
!*******************************************************************************


FUNCTION XYZ_TO_ELNUM(xco,yco,zco)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! Finds the element number in which a given x-y-z is located. 17 Jan 2001 MMB
!*******************************************************************************
   REAL(SP), INTENT(IN) :: xco,yco,zco
   INTEGER(I4B) :: XYZ_TO_ELNUM

   REAL(SP), DIMENSION(4) :: cart_vec,simp_vec
   REAL(SP), DIMENSION(4,4) :: cart_mat,simp_mat
   INTEGER(I4B) :: ielem
         
   XYZ_TO_ELNUM = 0
   cart_vec = (/ xco , yco , zco , 1.0 /)

   inside_element: DO ielem = 1,num_elements
     CALL SIMPLEX_COEFFICIENTS(ielem,simp_mat(1:4,4),simp_mat(1:4,1), &
                               simp_mat(1:4,2),simp_mat(1:4,3),cart_mat)
     simp_vec = MATMUL(simp_mat,cart_vec)
     IF (ABS(SUM(ABS(simp_vec))-1.0).LT.EPS) THEN ! inside element(ielem)
       XYZ_TO_ELNUM = ielem
       EXIT inside_element
     END IF
   END DO inside_element

END FUNCTION XYZ_TO_ELNUM
!*******************************************************************************




END MODULE geometry




