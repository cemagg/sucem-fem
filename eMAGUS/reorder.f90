! Bug somewhere around 2b - appears to be in subroutine....
! Last changed 6 Mar 00

SUBROUTINE REORDER
  USE geometry
  USE renumbering
  USE feminterface, ONLY: count_nonzeros, exchange_layeredges
  USE unit_numbers
  USE problem_info
  USE matrix
  IMPLICIT NONE
!*************************************************************************
!  Subroutine description
!*************************************************************************

! This routine attempts to minimise the bandwidth of the system matrices
! by a renumbering technique.  The algorithm is described in : An 
! Introduction to Finite Element Computations by Hinton and Owen 
! and and also in F. Meyer's masters and PhD theses.
!
!*************************************************************************
! AUTHOR
!*************************************************************************

!Riana Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
! 17 May 1999 10:00 RHG
! 17 Nov 14h30 DBD. USE statement modified to conform to MIPS 7.2.1.
! 25 Nov 17h30 DBD. Double de-allocation of first_layer removed.
!  6 Mar 2000 DBD. Updated to use new data structures.
!
!*************************************************************************
! Input
!*************************************************************************
!
! element_edges and edge_nodes, edgeconnectedge and edgeconnectelem
!
!*************************************************************************
! Output
!*************************************************************************
!
! element_edges renumbered and edge_nodes adjusted accordingly;
! also updated edgeconnectedge and edgeconnectelem
!
!*************************************************************************

! INTEGER(I4B) :: edge
INTEGER(I4B) :: layer_size, least, leastpos, search, edge_counter
INTEGER(I4B) :: current_layersize,  new_edge_one, old_edge, k, bandwd  
INTEGER(I4B) :: display_counter, layer_count, kedge, pos
INTEGER(I4B) :: edge_search, node_search, find_max, first_layer_size
INTEGER(I4B) :: temp_pos, layer_elem,conn_search, search_edge, listsize,connpos
INTEGER(I4B) :: tempval_int, i_sort, stopval,maxdim,maxsize

INTEGER(I4B), DIMENSION(3) :: index

REAL(SP) :: first_nodecor, second_nodecor, xmax,xmin,ymax,ymin,zmax,zmin
REAL(SP) :: zlength, ylength, xlength, maxval, width, edge_length, length_sum
REAL(SP) :: tempval_real, maxcor, mincor, value

REAL(SP), DIMENSION(3) :: start_p, end_p

INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: order, conn_list, first_layer, temp_layer
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: current_layer, next_layer, in_all_layers, list
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: edgeconn_list, this_ring, this_ring_and_more,new_list
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: in_prev_ring
CHARACTER(FILENAMELENGTH) ::fname

INTEGER(I4B) :: first_edge, temp_layer_pos, i_search, start, maxconnect
INTEGER(I4B) :: more_listsize, listpos, insert_i, this_ring_size, prev_first_edge

LOGICAL(LGT) DEBUG_OUTPUT, DEBUG_CONNECT

!*******************************************************************
! note that allocations of the type list = 0, where list has dimension
!larger than one, is in fact an array allocation.
!*******************************************************************
fname = 'atest.out'

OPEN (UNIT=TESTFILE,STATUS='UNKNOWN',file=fname)
maxconnect = SIZE(edgeconnectedge,2)
maxdim = num_edges ! this is a rough maximum that should be much larger than 
                   ! any layer could be      

allocate(order(num_edges))
allocate(conn_list(maxdim))
allocate(first_layer(maxdim))
first_layer = 0
allocate(next_layer(num_edges))

 
!*****************************************************************
! A first layer with more than one element is selected in this 
! case (3D reordering)
! In the case of a flat exterior surface, the whole surface is selected
! however, for a curved surface this is problematic.  There might
! not even be one edge that lies entirely in this exterior plane.
! To solve this problem, and to provide a renumbering that is as
! general as possible, some options are examined and the best one 
! then selected.
!*****************************************************************
 
!****************************************************************

!search for maximum and minimum coordinates, to form a cube around
!the mesh.  One of the planes  opposite the maximum length axis,
!is selected and the edges in this plane forms the first layer 

!****************************************************************
 maxsize = size(edgeconnectedge)
 xmax = vertices(1)%coord(1)
 xmin = vertices(1)%coord(1)
 ymax = vertices(1)%coord(2)
 ymin = vertices(1)%coord(2)
 zmax = vertices(1)%coord(3)
 zmin = vertices(1)%coord(3)
 
 DO node_search = 1,num_nodes
   !find x max and min coordinate  
   IF (vertices(node_search)%coord(1).GT.xmax) THEN
     xmax = vertices(node_search)%coord(1)
   END IF
   IF (vertices(node_search)%coord(1).LT.xmin) THEN   
     xmin = vertices(node_search)%coord(1)
   END IF
     
   !find y max and min coordinate  
   IF (vertices(node_search)%coord(2).GT.ymax) THEN
     ymax = vertices(node_search)%coord(2)
   END IF
   IF (vertices(node_search)%coord(2).LT.ymin) THEN   
     ymin = vertices(node_search)%coord(2)
   END IF

   !find z max and min coordinate  
   IF (vertices(node_search)%coord(3).GT.zmax) THEN
     zmax = vertices(node_search)%coord(3)
   END IF
   IF (vertices(node_search)%coord(3).LT.zmin) THEN   
     zmin = vertices(node_search)%coord(3)
   END IF
 END DO

 zlength = zmax - zmin
 ylength = ymax - ymin
 xlength = xmax - xmin
 
 length_array(1) = xlength
 length_array(2) = ylength
 length_array(3) = zlength
 
 
 maxval = xlength
 maxpos =1

 ! do general sort here, since we want the minumum direction as well
 index = 0
 index = (/ (k, k=1,3) /)
 sort : DO stopval = 3, 2, -1
   DO i_sort = 1,stopval-1
     IF (length_array(i_sort+1).LE.length_array(i_sort)) THEN
       tempval_real = length_array(i_sort+1)
       length_array(i_sort+1) = length_array(i_sort)
       length_array(i_sort) = tempval_real
       tempval_int = index(i_sort+1)
       index(i_sort+1) = index(i_sort)
       index(i_sort) = tempval_int
     END IF
   END DO
 END DO sort   
 maxpos = index(3)
 midpos = index(2)
 minpos = index(1) 
 IF (maxpos.LT.1.OR.maxpos.GT.3) THEN 
   STOP 'IE: In REORDER.'
 END IF


 min_array = (/xmin, ymin, zmin/)
 max_array = (/xmax, ymax, zmax /)
 
 first_layer = 0   
 pos = 0
 DO edge_search = 1,num_edges
   first_nodecor  = vertices(edges(edge_search)%nodes(1))%coord(maxpos)
   second_nodecor = vertices(edges(edge_search)%nodes(2))%coord(maxpos)
   IF ((first_nodecor.LE.(min_array(maxpos))).OR.(second_nodecor.LE.(min_array(maxpos)))) THEN
   !IF (first_nodecor.LE.min_array(maxpos)) THEN
     !IF (second_nodecor.LE.min_array(maxpos)) THEN
       pos = pos +1
       first_layer(pos) = edge_search
       
     !END IF
   END IF 
 END DO
WRITE(TESTFILE,*) 'first layer size : ', pos
IF (pos.GT.maxsize) THEN
  WRITE(TESTFILE,*) 'ERROR, MAX LAYER SIZE EXCEEDED FOR FIRST LAYER, PLEASE INCREASE MAX'
  STOP 'IE: MAX LAYER SIZE EXCEEDED FOR FIRST LAYER, contact author'
END IF 

!*********************************************************************
!For a more complex shape, the first layer requirements are different
!********************************************************************* 
!calculate the average edge length
length_sum = 0_SP

DO edge_search = 1,num_edges
  start_p = vertices(edges(edge_search)%nodes(1))%coord
  end_p   = vertices(edges(edge_search)%nodes(2))%coord
   
  edge_length =  SQRT( (start_p(1)-end_p(1))*(start_p(1)-end_p(1)) + &
        & (start_p(2)-end_p(2))*(start_p(2)-end_p(2))  + &
        & (start_p(3)-end_p(3))*(start_p(3)-end_p(3)))
 
  length_sum = length_sum + edge_length 
END DO
 
avg_edge_length = length_sum/num_edges
WRITE(TESTFILE,*) 'ave_length : ', avg_edge_length
 
IF (pos.LE.5) THEN
  WRITE(TESTFILE,*) 'The first layer is still empty, so relax the requirements a bit' 
  !The first layer is still empty, so relax the requirements a bit...
  first_layer = 0   
  pos = 0
  DO edge_search = 1,num_edges
    first_nodecor  = vertices(edges(edge_search)%nodes(1))%coord(maxpos)
    second_nodecor = vertices(edges(edge_search)%nodes(2))%coord(maxpos)
    
    width = avg_edge_length*2_SP
    IF ((first_nodecor.LE.(min_array(maxpos)+width)).AND.(second_nodecor.LE.(min_array(maxpos)+width))) THEN
      pos = pos +1
      IF (pos.GT.maxsize) THEN
        STOP 'IE: MAXIMUM REORDER LAYER SIZE EXCEEDED.'
      END IF      
      first_layer(pos) = edge_search
      
    END IF 
  END DO
END IF  
first_layer_size = pos


!************************************************************* 
!proceed with renumbering algorithm
!*************************************************************
  
  allocate(in_all_layers(num_edges))
  in_all_layers = 0
  CALL count_nonzeros(SIZE(first_layer),first_layer,current_layersize)
  allocate(current_layer(current_layersize))
  current_layer = 0
  current_layer = first_layer(1:current_layersize)
  deallocate(first_layer)
  maxval = SIZE(current_layer)
  CALL count_nonzeros(SIZE(current_layer),current_layer,current_layersize)
  edge_counter = 0
  layer_count = 0

  WALK_LAYERS : DO
     layer_count = layer_count + 1 
     debug_var_int = layer_count
     WRITE(TESTFILE,*) 'the number of this layer is : ', layer_count
     ascending = .TRUE.
     CALL exchange_layeredges(current_layer,next_layer,in_all_layers,& 
          edge_counter,layer_count)
     ! The current layer completed, proceed to the next layer:
     
     CALL count_nonzeros(SIZE(next_layer),next_layer,current_layersize)
     IF (current_layersize.EQ.0) THEN
       EXIT WALK_LAYERS
     END IF  
     deallocate(current_layer)
     allocate(current_layer(current_layersize))
     current_layer = 0
     current_layer = next_layer(1:current_layersize)
     IF (edge_counter.EQ.num_edges) THEN
      ! The renumbering is complete
      EXIT walk_layers
     END IF
  END DO WALK_LAYERS
  
  deallocate(order)
  deallocate(conn_list)
!  deallocate(first_layer) ! Corrected DBD 25 Nov, already de-allocated above.
  deallocate(current_layer)
  deallocate(next_layer)
  deallocate(in_all_layers)
  WRITE(TESTFILE,*) 'layer_count : ', layer_count
  !close(TESTFILE)  ! this will be closed in femfeko

END SUBROUTINE REORDER

SUBROUTINE remove_multiple(listsize, list, newlist)

!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! A list of edge connections was generated and serves as an input to this routine.
! Some of the edges occur multiply in this list and the additional entries must
! be removed.  Indexlist holds ones where an edge occurs that has already been
! found in the list.

!*******************************************************************************
! AUTHOR
!*******************************************************************************
!
!Riana Geschke
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
!
!9 May 00:10 RHG
!17 Nov 14h30 DBD. IMPLICIT NONE added. 
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
!
!The list is the edge connection list.  Listsize is the number of consecutive 
!nonzero entries starting form the first entry in the list.  
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
!
!The list with each edge occuring only once, is called newlist
!
!*******************************************************************************

USE unit_numbers
USE problem_info
IMPLICIT NONE
 
INTEGER(I4B) :: listsize
INTEGER(I4B), DIMENSION(listsize), INTENT(IN)  :: list 
INTEGER(I4B), DIMENSION(listsize), INTENT(OUT) :: newlist

INTEGER(I4B) :: k_search, edge_index, insert_i, pos, edge
INTEGER(I4B), DIMENSION(listsize) :: indexlist, templist

templist = 0
indexlist = 0 
DO edge_index = 1,listsize
  edge = list(edge_index)
  DO k_search = edge_index,listsize
    IF (list(k_search).EQ.edge) THEN
      IF (edge_index.NE.k_search) THEN
        indexlist(k_search) = 1  
      END IF
    END IF
  END DO
END DO
pos = 0
DO insert_i = 1,listsize
  IF (indexlist(insert_i).EQ.0) THEN
    pos = pos + 1
    templist(pos) = list(insert_i)
  END IF
END DO
newlist = templist
END SUBROUTINE remove_multiple

