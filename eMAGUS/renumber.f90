! Last changed 09 Feb 2001 MMB - FJCM updates

SUBROUTINE exchange_edges(new_edge,old_edge)

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine swaps the numbers of two edges.  This requires updating of 
! all the structures listed below.
!
!*************************************************************************
! AUTHOR
!*************************************************************************
!
!Riana Helena Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
!  9 May 1999 00:20 RHG
!  17 Nov 14h45 DBD. USE statement changed to conform to MIPS 7.2.1
!  6 March 2000 DBD. Updated to use new data structures.
!
!*************************************************************************
! Input
!*************************************************************************
!
! The numbers of the two edges that must be swapped
!
!*************************************************************************
! Output
!*************************************************************************
! 
! The following are updated:
!  element_edges 
!  edgeconnectedge
!  edgeconnectelem
!  edge_nodes 
!*************************************************************************


  USE problem_info
  USE geometry
  USE renumbering
  USE unit_numbers
  USE feminterface, ONLY: COUNT_NONZEROS
  USE matrix
  IMPLICIT NONE
  
  INTEGER(I4B), INTENT(IN) :: new_edge, old_edge
  INTEGER(I4B) iedge
  INTEGER(I4B) renum, newlistsize, oldlistsize,  k, listsize, maxconnelems, maxconnedges
  INTEGER(I4B) display_counter, pos, search_counter, insert_counter
  INTEGER(I4B) replace_c, edge_c, found
  INTEGER(I4B),DIMENSION(2) :: tempnodes
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: new_elemconnlist, old_elemconnlist
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: done_list, edge_conn_list, edge_sat
!  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: new_edgeconnlist, old_edgeconnlist
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: temp, old_conn_edges, new_conn_edges

  
  IF (new_edge.GT.num_edges) THEN
    STOP 'The edge_numbers are corrupted.'
  END IF  
  
  
  maxconnelems = SIZE(edgeconnectelem,2)
  maxconnedges = SIZE(edgeconnectedge,2)
  allocate(old_elemconnlist(maxconnelems))
  allocate(new_elemconnlist(maxconnelems)) 
!  allocate(old_edgeconnlist(maxconnelems))
!  allocate(new_edgeconnlist(maxconnedges))
  allocate(done_list(num_edges))
  allocate(temp(maxconnedges))
  new_elemconnlist = edgeconnectelem(new_edge,:)
  old_elemconnlist = edgeconnectelem(old_edge,:)
  CALL count_nonzeros(maxconnelems,new_elemconnlist,newlistsize)
  CALL count_nonzeros(maxconnelems,old_elemconnlist,oldlistsize)
  
  !**********************************************************
  ! element_edges
  !**********************************************************
  done_list = 0
  DO renum = 1,newlistsize
    done_list(new_elemconnlist(renum)) = 1
    DO iedge = 1,6  
      IF (elements(new_elemconnlist(renum))%edges(iedge).EQ.new_edge) THEN
        elements(new_elemconnlist(renum))%edges(iedge) = old_edge
      ELSE IF (elements(new_elemconnlist(renum))%edges(iedge).EQ.old_edge) THEN 
        elements(new_elemconnlist(renum))%edges(iedge) = new_edge
      END IF  
    END DO
  END DO
  DO renum = 1,oldlistsize
    !be careful not to reverse the changes above: 
    IF (done_list(old_elemconnlist(renum)).EQ.0) THEN
      DO iedge = 1,6  
        IF (elements(old_elemconnlist(renum))%edges(iedge).EQ.new_edge) THEN
          elements(old_elemconnlist(renum))%edges(iedge) = old_edge
        ELSE IF (elements(old_elemconnlist(renum))%edges(iedge).EQ.old_edge) THEN 
          elements(old_elemconnlist(renum))%edges(iedge) = new_edge
        END IF  
      END DO
    END IF 
  END DO
  
  !*********************************************************
  ! edgeconnectedge
  !*********************************************************
  
  ! similarly here the entries in the structure that contain
  ! references to the two edge numbers, must all be swopped.
  ! additionally, after this was done, the two rows involved
  ! must be exchanged
 
  allocate(old_conn_edges(maxconnedges))
  allocate(new_conn_edges(maxconnedges)) 
  old_conn_edges = edgeconnectedge(old_edge,:)
  new_conn_edges = edgeconnectedge(new_edge,:)
  CALL count_nonzeros(SIZE(old_conn_edges),old_conn_edges,oldlistsize)
  CALL count_nonzeros(SIZE(new_conn_edges),new_conn_edges,newlistsize)  
  allocate(edge_conn_list(oldlistsize+newlistsize))
  edge_conn_list = 0
  edge_conn_list(1:oldlistsize) = old_conn_edges(1:oldlistsize)
  pos = oldlistsize
  DO insert_counter = 1,newlistsize  
    found = 0
    DO search_counter = 1,oldlistsize
      IF (new_conn_edges(insert_counter).EQ.old_conn_edges(search_counter)) THEN
        found = 1
      END IF
    END DO
    IF (found.EQ.0) THEN
      pos = pos + 1
      edge_conn_list(pos) =  new_conn_edges(insert_counter)
    END IF 
  END DO

  DO replace_c = 1,pos
    DO edge_c = 1,maxconnedges
      IF (edgeconnectedge(edge_conn_list(replace_c),edge_c).EQ.old_edge) THEN
              edgeconnectedge(edge_conn_list(replace_c),edge_c) = new_edge
      ELSE IF (edgeconnectedge(edge_conn_list(replace_c),edge_c).EQ.new_edge) THEN
        edgeconnectedge(edge_conn_list(replace_c),edge_c) = old_edge
      END IF
    END DO
  END DO 
  
  !the swopping of rows:
  temp = edgeconnectedge(old_edge,:)
  edgeconnectedge(old_edge,:) = edgeconnectedge(new_edge,:)
  edgeconnectedge(new_edge,:)  = temp

  deallocate(edge_conn_list)
  deallocate(old_conn_edges)
  deallocate(new_conn_edges) 

  !*********************************************************
  ! edgeconnectelem
  !*********************************************************
  ! simply swop the two rows :  

  ! FJCM
  ! was 
  ! temp = edgeconnectelem(old_edge,:)
  ! is
  temp(1:SIZE(edgeconnectelem, DIM=2)) = edgeconnectelem(old_edge,:)

  edgeconnectelem(old_edge,:) = edgeconnectelem(new_edge,:)
  edgeconnectelem(new_edge,:)  = temp 
  !*********************************************************
  ! edge_nodes
  !*********************************************************
  
  tempnodes = edges(old_edge)%nodes
  edges(old_edge)%nodes = edges(new_edge)%nodes
  edges(new_edge)%nodes = tempnodes

  
  deallocate(new_elemconnlist)
  deallocate(old_elemconnlist)
  deallocate(done_list)
!  deallocate(new_edgeconnlist)
!  deallocate(old_edgeconnlist) 
  deallocate(temp)
  
END SUBROUTINE exchange_edges 




SUBROUTINE exchange_layeredges(current_layer, next_layer,in_this_layer,edge_counter,layer_number)


!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine renumbers an entire layer. It assigns new edge numbers to
! each edge in the layer, and then uses exchange_edges to swop the edge
! numbers ( new and old)
!
!*************************************************************************
! AUTHOR
!*************************************************************************
!
!Riana Helena Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
!   9 May 99 00:25 RHG
!  17 Nov 14h45 DBD. USE statement changed to conform to MIPS 7.2.1
!
!*************************************************************************
! Input
!*************************************************************************
!
! The current layer (all its edges) and the number of edges already
! renumbered (edge_counter). 
!
!*************************************************************************
! Output
!*************************************************************************
!
! The current layer in its renumberd form, as well as the next layer to be
! renumbered.
!
!*************************************************************************

USE feminterface, ONLY: COUNT_NONZEROS
USE geometry
USE renumbering
USE unit_numbers
USE problem_info

IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: current_layer
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: next_layer
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: in_this_layer
  INTEGER(I4B), INTENT(INOUT) :: edge_counter
  INTEGER(I4B), INTENT(IN) :: layer_number
  

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE  :: this_layer
  INTEGER(I4B) :: layersize, k_count, neighbour_count, all_count, k, maxsize
  INTEGER(I4B) :: stopval, i_sort,i_element,k_element,i_count, tempval, renumindex, ren_pos, sum
  INTEGER(I4B) :: old_edgenum, edge_num, pos, neighboursize, search
  INTEGER(I4B) :: same_layer
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::order, index,neighbour_connlist
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::neighbours, updated_this_layer
  INTEGER(I4B), DIMENSION(num_edges) :: edge_in_new_layer

  maxsize = SIZE(edgeconnectedge,2)
  allocate(neighbours(maxsize))

  CALL count_nonzeros(SIZE(current_layer),current_layer,layersize)
  
  allocate(index(layersize)) 
  allocate(neighbour_connlist(num_edges))
  allocate(this_layer(layersize))
  this_layer  = 0
  allocate(updated_this_layer(layersize))
  updated_this_layer = 0
  this_layer = current_layer
  
  ! disabled for now, this will be uncommented when ready
  
  ! FJCM (The 4 lines below was NOT commented out as indicated in the text a
  !IF (layer_number.EQ.1) THEN 
  !  CALL layer_the_layer(layersize,this_layer,updated_this_layer,layer_number)
  !  this_layer = updated_this_layer
  !END IF 
  
  !*********************************************************
  ! walk through layer and renumber
  !*********************************************************
  WRITE(TESTFILE,*) 'this_layer in exchange_layeredges : ', this_layer
  
  reorder_layer: DO ren_pos = 1,layersize 
    edge_counter = edge_counter +1    
    old_edgenum = this_layer(ren_pos)
    edge_num = edge_counter

    IF (old_edgenum.GE.edge_num) THEN
      !WRITE(TESTFILE,*) 'new: ', edge_num
      !WRITE(TESTFILE,*) 'old: ', old_edgenum
      IF (old_edgenum.NE.edge_num) THEN
        CALL exchange_edges(old_edgenum,edge_num)
      END IF  
      DO search = ren_pos+1,layersize
        IF (this_layer(search).EQ.edge_num) THEN
          this_layer(search) = old_edgenum
        END IF  
      END DO   
      this_layer(ren_pos) = edge_num      
    ELSE
      ! no exchange was made
      edge_counter = edge_counter - 1
    END IF
  END DO reorder_layer
  
  !**********************************************************
  ! identify the elements of the next layer 
  ! and update current_layer
  !**********************************************************
  next_layer = 0
  IF (edge_counter.LT.num_edges) THEN
    pos = 0
    sum = 0
    edge_in_new_layer = 0
    DO i_element = 1, layersize
      neighbours = edgeconnectedge(this_layer(i_element),:) 
      maxsize = size(neighbours)
      CALL count_nonzeros(SIZE(neighbours),neighbours,neighboursize)
      DO k_element = 1,neighboursize
        IF (neighbours(k_element).GT.edge_counter) THEN  
          IF (edge_in_new_layer(neighbours(k_element)).EQ.0) THEN
            pos = pos + 1
            next_layer(pos) = neighbours(k_element)
            edge_in_new_layer(neighbours(k_element)) = 1
          END IF
        END IF
      END DO
    END DO
  END IF 

  deallocate(index)   
  deallocate(neighbour_connlist)
  deallocate(this_layer)
  deallocate(updated_this_layer)

END SUBROUTINE exchange_layeredges



SUBROUTINE COUNT_NONZEROS(listsize,list,nonzeros)
  USE nrtype 
  IMPLICIT NONE
!*************************************************************************
!  Subroutine description
!*************************************************************************
! This routine counts the number of consecutive
! non-zero elements in a list. At the first zero,
! counting is stopped, and any other non_zero 
! elements following are ignored.
!
!*************************************************************************
! AUTHOR
!*************************************************************************
!
!Riana Helena Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
!9 May 99 00:25 RHG
!*************************************************************************
  INTEGER(I4B),INTENT(IN) :: listsize
  INTEGER(I4B), DIMENSION(listsize),INTENT(IN) :: list
  INTEGER(I4B), INTENT(OUT) :: nonzeros

  INTEGER(I4B) nonzero_counter, pos, k 

  pos = 0
 
  DO k = 1,listsize
    IF (list(k).EQ.0) THEN
      pos = k
      EXIT  
    END IF
  END DO
 
  IF (pos.EQ.0) THEN
    nonzeros = listsize
  ELSE  
    nonzeros = k-1
  END IF 
  
END SUBROUTINE COUNT_NONZEROS
!*************************************************************************


SUBROUTINE layer_the_layer(this_layer_size, this_layer,updated_this_layer,layer_number)

!*************************************************************************
!  Subroutine description
!*************************************************************************
! THIS ROUTINE IS NOT COMPLETED YET AND IS CURRENTLY NOT USED...
!
! Each layer is treated as a separate layer in itself. The layers within
! a layer are referred to as rings, to distinguish between the two types 
! of layers.
!
!*************************************************************************
! AUTHOR
!*************************************************************************
!
!Riana Helena Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
! 9 May 99 00:25 RHG
! 17 June - June RHG  further improvement on bandwidth reduction
! 17 Nov 14h45 DBD. USE statement changed to conform to MIPS 7.2.1
! 
!
!*************************************************************************
! Input
!*************************************************************************
!
! The layer, to be layered.
!
!*************************************************************************
! Output
!*************************************************************************
!
! The same layer and the same edges in it, only in a different order.
!
!*************************************************************************


USE geometry
USE renumbering
USE feminterface, ONLY: COUNT_NONZEROS
USE unit_numbers
USE problem_info
USE math_tools, ONLY: FIND_IN_LIST

IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: this_layer_size
INTEGER(I4B), DIMENSION(this_layer_size), INTENT(IN) :: this_layer
INTEGER(I4B), DIMENSION(this_layer_size), INTENT(OUT) :: updated_this_layer
INTEGER(I4B), INTENT(IN) :: layer_number

INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: this_ring_and_more, index
INTEGER(I4B), DIMENSION(:), ALLOCATABLE::new_list,list
REAL(SP), DIMENSION(:), ALLOCATABLE::edge_mincor
INTEGER(I4B), DIMENSION(this_layer_size) :: temp_layer, this_ring, temp_this_ring
INTEGER(I4B), DIMENSION(num_edges) :: in_prev_ring
INTEGER(I4B) first_edge,  edge_search, pos, start, insert_i, listsize,maxconnedges
INTEGER(I4B) iedge
INTEGER(I4B) listpos, i_search, more_listsize, this_ring_size, temp_layer_pos, walk_count
INTEGER(I4B) closest_edge ,k
REAL(SP)  first_nodecor, second_nodecor, maxcor, mincor, value, distance, prev_distance
REAL(SP), DIMENSION(3) :: startpoint, corpos

 
 !************************************************************************
 !
 ! Step 1 is to find the edge that must be the first edge in the current
 ! layer.  When this is the first layer, it is selected as the closest edge
 ! to a specified point, however for the other layers, it is simply the
 ! edge that is the first one in the list.
 !
 !************************************************************************

 first_edge = 0  ! initialisation
 
! IF (layer_number.EQ.1) THEN
  
   WRITE(TESTFILE,*) 'select a first edge for the first layer: '
   startpoint(maxpos) = min_array(maxpos) 
   startpoint(minpos) = min_array(minpos) + &
       &   0.5_SP*(max_array(minpos)-min_array(minpos))
   startpoint(midpos) = min_array(midpos) 
   WRITE(TESTFILE,*) 'Looking for closest edge to this point : ', startpoint
   closest_edge = this_layer(1)
   corpos(1) = MAX(vertices(edges(closest_edge)%nodes(1))%coord(1) , & 
               vertices(edges(closest_edge)%nodes(2))%coord(1))
   corpos(2) = MAX(vertices(edges(closest_edge)%nodes(1))%coord(2) , & 
               vertices(edges(closest_edge)%nodes(2))%coord(2))
   corpos(3) = MAX(vertices(edges(closest_edge)%nodes(1))%coord(3) , & 
               vertices(edges(closest_edge)%nodes(2))%coord(3))
   distance = SQRT((corpos(1) - startpoint(1))*(corpos(1) - startpoint(1)) + & 
     & (corpos(2) - startpoint(2))*(corpos(2) - startpoint(2)) + &
     & (corpos(3) - startpoint(3))*(corpos(3) - startpoint(3)))
   prev_distance = distance
   
   DO edge_search = 2, this_layer_size
     iedge = this_layer(edge_search)
     !corpos(1) = (vertices(iedges(iedge)%nodes(1),1) + vertices(iedges(iedge)%nodes(2),1))/2
     !corpos(2) = (vertices(iedges(iedge)%nodes(1),2) + vertices(iedges(iedge)%nodes(2),2))/2
     !corpos(3) = (vertices(iedges(iedge)%nodes(1),3) + vertices(iedges(iedge)%nodes(2),3))/2 
     corpos(1) = MAX(vertices(edges(iedge)%nodes(1))%coord(1), & 
                 vertices(edges(iedge)%nodes(2))%coord(1))
     corpos(2) = MAX(vertices(edges(iedge)%nodes(1))%coord(2), &
                 vertices(edges(iedge)%nodes(2))%coord(2))
     corpos(3) = MAX(vertices(edges(iedge)%nodes(1))%coord(3), &
                 vertices(edges(iedge)%nodes(2))%coord(3)) 
     distance = SQRT((corpos(1) - startpoint(1))*(corpos(1) - startpoint(1)) + &
       & (corpos(2) - startpoint(2))*(corpos(2) - startpoint(2)) + &
       & (corpos(3) - startpoint(3))*(corpos(3) - startpoint(3)))
       
     IF (distance.LT.prev_distance) THEN
       
       closest_edge = iedge
       prev_distance = distance
     END IF   
   END DO
   first_edge = closest_edge 
   IF (first_edge.EQ.0) THEN
     STOP 'first edge in layer not found.'
   END IF
 !ELSE
 !  WRITE(TESTFILE,*) 'This is not the first layer, and the first edge is known, it is :', this_layer(1)
 !  first_edge = this_layer(1)
 !END IF

 maxconnedges = SIZE(edgeconnectedge,2) ! Corrected DBD
 pos = 0
 allocate(this_ring_and_more(maxconnedges*150))
 this_ring_and_more = 0 ! Array allocation, force all unconnected.
 this_ring = 0
 temp_this_ring = 0
 allocate(new_list(num_edges))
 allocate(list(maxconnedges))
 temp_layer = 0
 new_list = 0
 
 this_ring_and_more(1:maxconnedges) = edgeconnectedge(first_edge,:)
!print *,edgeconnectedge(first_edge,:)
!print *,this_ring_and_more
 
 in_prev_ring = 0
 in_prev_ring(first_edge) = 1
 pos = 0
 temp_layer(1) = first_edge
 temp_layer_pos = 2
 walk_count = 0
 walk_rings: DO
   walk_count = walk_count + 1
   IF (walk_count.GT.80) THEN  ! 80 is an arbitrary choice, but 10 should be max
     STOP 'Too many rings in a single layer.  Sonething wrong here... '
   END IF  
   CALL count_nonzeros(SIZE(this_ring_and_more),this_ring_and_more,more_listsize)
   !*****************************************************************************
   ! sort these edges with respect to the minpos coordinate...
   !*****************************************************************************

   allocate(edge_mincor(more_listsize))
   DO k = 1,more_listsize
     edge_mincor(k)=(vertices(edges(this_ring_and_more(k))%nodes(1))%coord(minpos)+ &
                     vertices(edges(this_ring_and_more(k))%nodes(2))%coord(minpos))/2
   END DO
   allocate(index(more_listsize))
   IF (ascending) THEN
     ascending = .FALSE.
   ELSE 
     ascending = .TRUE.
   END IF
   WRITE(TESTFILE,*) 'edge_mincor : ', edge_mincor
   CALL sort_reallist(more_listsize,edge_mincor,index)
   this_ring_and_more(1:SIZE(index)) = this_ring_and_more(index)
   deallocate(index)
   deallocate(edge_mincor)
   
   pos = 0
   DO i_search = 1, more_listsize
     CALL FIND_IN_LIST(this_layer,this_layer_size,this_ring_and_more(i_search),listpos) 
     IF (listpos.NE.0) THEN
       IF (this_layer(listpos).NE.this_ring_and_more(i_search)) THEN
          STOP 'FIND_IN_LIST error'
       END IF   
       ! the edge is found in the current layer 
       IF (in_prev_ring(this_layer(listpos)).EQ.0) THEN
         !the edge has not occurred in a previous ring of this layer
         pos = pos + 1
         in_prev_ring(this_layer(listpos)) = 1
         this_ring(pos) = this_layer(listpos) 
       END IF
     END IF  
   END DO
   
   IF (pos.EQ.0) THEN
     exit walk_rings
   END IF   
   CALL count_nonzeros(SIZE(this_ring),this_ring,this_ring_size)
   IF ((temp_layer_pos+this_ring_size-1).GT.this_layer_size) THEN
     STOP 'error in layer_the_layer'
   END IF  
   
   IF (temp_layer_pos.GT.this_layer_size) THEN
     
     EXIT walk_rings
   END IF  
   temp_layer(temp_layer_pos : temp_layer_pos + this_ring_size - 1) = this_ring(1:this_ring_size)
   temp_layer_pos = temp_layer_pos + this_ring_size
   
   this_ring_and_more = 0
   start = 1
   DO insert_i = 1,this_ring_size
     ! FJCM
     ! Was
     ! list = edgeconnectedge(this_ring(insert_i),:)          
     ! is
     ! list = 0  ! also added by FJCM   
     list(1:SIZE(edgeconnectedge, DIM=2)) = edgeconnectedge(this_ring(insert_i),:)
     CALL count_nonzeros(SIZE(list),list,listsize) 
     this_ring_and_more(start:start+listsize-1) = list(1:listsize)
     start = start + listsize
   END DO 
   !eliminate multiple entries...
   new_list = 0
   CALL count_nonzeros(SIZE(this_ring_and_more),this_ring_and_more,more_listsize)
   CALL remove_multiple(more_listsize,this_ring_and_more,new_list)
   this_ring_and_more = 0

   ! FJCM
   ! Was
   ! this_ring_and_more = new_list      
   ! is
   this_ring_and_more(1:SIZE(new_list)) = new_list      

   temp_this_ring = this_ring
   this_ring = 0
 END DO walk_rings

 updated_this_layer = temp_layer
 CALL count_nonzeros(SIZE(updated_this_layer),updated_this_layer,listsize)
 IF (listsize.NE.this_layer_size) THEN
   
   updated_this_layer = this_layer
 END IF
 deallocate(this_ring_and_more)
 deallocate(new_list)
 deallocate(list)

 RETURN  

CONTAINS ! Internal subprogram.

SUBROUTINE SORT_REALLIST(n,list,index)
USE nrtype
USE renumbering

IMPLICIT NONE

INTEGER(I4B), INTENT(IN):: n
REAL(SP) , DIMENSION(n), INTENT(INOUT) :: list
INTEGER(I4B), DIMENSION(n), INTENT(OUT) :: index


INTEGER(I4B) i_sort, stopval, k, tempval_int
! Error corrected DBD 1 May 2003 REAL(I4B) tempval_real
REAL(SP) tempval_real
 index = 0
 index = (/ (k, k=1,n) /)
 
 IF (ascending) THEN
  sort_asc : DO stopval = n, 2, -1
   DO i_sort = 1,stopval-1
     IF (list(i_sort+1).GE.list(i_sort)) THEN
       tempval_real = list(i_sort+1)
       list(i_sort+1) = list(i_sort)
       list(i_sort) = tempval_real
       tempval_int = index(i_sort+1)
       index(i_sort+1) = index(i_sort)
       index(i_sort) = tempval_int
     END IF
   END DO
 END DO sort_asc   

 ELSE
 sort_desc : DO stopval = n, 2, -1
   DO i_sort = 1,stopval-1
     IF (list(i_sort+1).LE.list(i_sort)) THEN
       tempval_real = list(i_sort+1)
       list(i_sort+1) = list(i_sort)
       list(i_sort) = tempval_real
       tempval_int = index(i_sort+1)
       index(i_sort+1) = index(i_sort)
       index(i_sort) = tempval_int
     END IF
   END DO
 END DO sort_desc   
 END IF
END SUBROUTINE SORT_REALLIST


END SUBROUTINE layer_the_layer
