! Last changed 24-27 Apr 2003: DBD. Skyline storage scheme added to permit use
!                              of Compaq Extended Math Library. New option 
!                              added to MAVECPROD_REAL
! Last changed 14 Dec 2002: DBD. New mixed potential iterative solvers added.  
! Last changed 05 Mar 2002: MATRIX_SPARSE_ALLOCATE changed to support
!                           QT/QN elements. 
! Last changed 27 Apr 2001: restart_GMRES added in GMRES_SOLVE. DBD. 
! Last changed 25 Apr 2001: USE_PRE_CONDITIONER added in ITER_SOLVE. DBD
! Last changed 16 Apr 2001: QMR possibility added in ITER_SOLVE. MMB
! Last changed 09 Feb 2001 MMB - FJCM updates
! Last changed 31 Jan 2001 MMB - Error handling 
! BiCG_SOLVE changed to ITER_SOLVE
! MATVECPROD changed  


SUBROUTINE CRS_TO_SKYLINE
   USE nrtype
   USE unit_numbers
   USE matrix
   USE problem_info
   IMPLICIT NONE              

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine converts a SYMMETRIC matrix stored (fully) in compressed row storage to one
! stored in skyline storage, for subsequent use in skyline solvers. 
! Note that it is ASSUMED that the matrix is symmetric, since checking this is 
! complex once it is stored in sparse format.
! The matrix entries are also converted to double precision, since the
! XML libraries operate on double-precision data.
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! David B Davidson
!*************************************************************************
!  Revision history:
!*************************************************************************
! First version: DBD 28 April 2003
!  
!*************************************************************************
!  Input
!*************************************************************************
!
! CRS data: val (AR), col_ind (JA) and row_ind (IA) (parenthetic names
! give Compaq Extended Math Library - XML- nomenclature).
!
!*************************************************************************
!  Output
!*************************************************************************
! Profile-in skyline matrix of the the upper-triangular part
! stored by columns (or equivalently, since the matrix is symmetric, 
! the lower-triangular part stored by rows - this routine actually uses 
! this approach). 
! Two arrays store the data:
! AU: This stores the matrix by rows from the first non-zero element to the diagonal.
! IAUDIAG: This stores the index of the diagonal elements in AU.
! Note that this routine assumes that there is AT LEAST a diagonal
! element on each row (if not, the matrix cannot be factored
! in any case). 
!  
!*************************************************************************
! References 
!*************************************************************************
! Compaq Extended Math Library Reference Guide July 1999,
! Chapter 11: "Using the Direct Solvers for Sparse Linear Systems".
!
!*************************************************************************
  INTEGER(I4B) :: counter, nnz, nnz_symmetrical,full_matrix_sym  
                      
  INTEGER(I4B) :: col,row, rowstart, rowstop, ii

!  INTEGER(I4B) :: thiscol,nextcol ! for debugging

  ! Firstly, establisH the overall size of the AU matrix.
  ! This is done on a row-by-row basis
  num_nonzeros_sky = 0
  DO row = 1,dof
    rowstart = row_ind(row)
    num_nonzeros_sky = num_nonzeros_sky+(row-col_ind(rowstart)+1) 	 
  END DO
  ! Allocate and initialize the skyline matrix. 
  ALLOCATE(Asparse_skyline(num_nonzeros_sky))
  ALLOCATE(IAUdiag(dof))
  Asparse_skyline = 0.0_DP ! Note this matrix is stored in Double Precision!
  ! Now fill this matrix from the CRS matrix. Note that 
  ! the structural zeros are taken care of by the above initialization.
  ! Only the non-zero entries are now processed. 

  IF (DEBUG_SKYLINE) THEN 
    WRITE(FILEOUT,*) 
    WRITE(FILEOUT,*) 'Matrix in CRS form'
	WRITE(FILEOUT,*) 'Row index IA',row_ind
	WRITE(FILEOUT,*) 'Col index JA',col_ind
	WRITE(FILEOUT,*) 'CRS matrix:',Asparse
  END IF
  
  counter = 1
  ! First row special treatment - only diagonal
  Asparse_skyline(counter) = DBLE(Asparse(1)) 
  IAUdiag(counter) = 1
  ! Now deal with all other rows
  row_loop: DO row = 2,dof
    rowstart = row_ind(row)
    rowstop  = row_ind(row+1)-1
	counter = counter + 1 
    Asparse_skyline(counter) = DBLE(Asparse(rowstart))
    IF (rowstart.EQ.rowstop.OR.col_ind(rowstart).EQ.row) THEN ! Only diagonal present
	                                                          ! or first element is diagonal
      IAUdiag(row) = counter	
    ELSE 
	  ! Only execute this loop if the first element isn't the diagonal.
      DO ii = rowstart,rowstop-1 
		! nextcol = col_ind(ii+1)
		! thiscol = col_ind(ii)
        counter = counter + col_ind(ii+1) - col_ind(ii) 
		IF (counter.GT.num_nonzeros_sky) THEN 
		  STOP 'Internal error in CRS_TO_SKYLINE'
		END IF
        Asparse_skyline(counter) = DBLE(Asparse(ii+1))
		IF (col_ind(ii+1).EQ.row) THEN  
          IAUdiag(row) = counter
		  EXIT ! Jump out of this loop, all entries up to the diagonal have been processed
	    END IF
	  END DO  
	END IF
  END DO row_loop

  IF (DEBUG_SKYLINE) THEN 
    WRITE(FILEOUT,*) 'Matrix in skyline storage form - only symmetric half'
	WRITE(FILEOUT,*) 'Diagonal index IAUdiag',IAUdiag
	WRITE(FILEOUT,*) 'Upper (lower) matrix:',Asparse_skyline
  END IF

 ! Compute and write some matrix statistics

  nnz             = row_ind(dof+1)-1    ! Number of non-zeros. 
  nnz_symmetrical = (nnz-dof)/2 + dof   ! Size required if symmetry used - lower bound on storage needed.
  full_matrix_sym = dof*(dof+1)/2       ! Storage required for full matrix, if only symmetrical U or L stored.

  WRITE(FILEOUT,'(//,18X,A)') 'PROFILE-IN STORAGE MODE STATISTICS'
  WRITE(FILEOUT,'(1X,A,T100,I12)') 'Number of entries for full symmetrical matrix, stored as either U (or L) : ',&
                              full_matrix_sym
  WRITE(FILEOUT,'(1X,A,T100,I12)') 'Number of non-zero entries (including structural zeros) of U (L) stored: ',&
                              num_nonzeros_sky
  WRITE(FILEOUT,'(1X,A,T100,I12)')  'Num of non-zeros entries as stored in Compressed Row Storage: ',nnz
  WRITE(FILEOUT,'(1X,A,T100,I12)')  'Num of non-zeros entries (assuming symmetry used): ',nnz_symmetrical
  WRITE(FILEOUT,'(1X,A,T100,F12.5,A)')  & 
        'Percentage non-zeros (using symmetry):', 100*REAL(nnz_symmetrical)/REAL(full_matrix_sym),'%'
  WRITE(FILEOUT,'(1X,A,T100,F12.5,A)')  & 
        'Percentage non-zeros using profile storage (ditto):',100*REAL(num_nonzeros_sky)/REAL(full_matrix_sym),'%'
  WRITE(FILEOUT,'(1X,A,T100,F12.5)')  & 
        'Overhead of profiled symmetrical storage cf. non-zero only - but no symmetry - storage: ',&
        REAL(num_nonzeros_sky)/REAL(nnz)
	
END SUBROUTINE CRS_TO_SKYLINE



SUBROUTINE SPARSEIG_RI  
   USE nrtype
   USE unit_numbers
   USE matrix
   USE problem_info
   USE bandwidth
   USE FEMINTERFACE, ONLY: fill_ones, find_greater_than, find_all_in_list, &
                        convertfillcor, convertcor, symbolic_fill
   IMPLICIT NONE              
 
!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine solves the generalised eigenvalue problem in the ARPACK
! regular inverse mode or shift-invert mode
! 
!  S and T are stored as sparse matrices
!
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
!  March 8 2000
!
!*************************************************************************
!  Input
!*************************************************************************
!
!  , Tsparse and indexing matrices
!
!*************************************************************************
!  Output
!*************************************************************************
!
!  
!
!*************************************************************************

 INTEGER(I4B) :: lworkl,info, sig, flag, val, k
 REAl(SP) :: tol , ss
 Logical(LGT),DIMENSION(:), allocatable:: select    
 Real(SP), DIMENSION(:,:), allocatable   ::  v, z
   
 Real(SP), DIMENSION(:), ALLOCATABLE  :: resid, d
 Integer(I4B), DIMENSION(11):: iparam
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: iwork
 Real(SP), DIMENSION(:), allocatable:: workl, workd                 
 character howmny
 REAL(SP), DIMENSION(:), ALLOCATABLE :: mvec, xvec, yvec,dummy
 
 REAL(SP) timerstart, timer1, timer2, timer2b, timer3, timer4, timer5 
 INTEGER(I4B) count, count_rate, count_max
 

! JPS 1 Aug 2005: Commented out calls to subroutines required for
! direct solution of ARPACK shift invert linear systems. These
! systems are now solved using preconditioned conjugate gradients

    CALL SYSTEM_CLOCK(count,count_rate,count_max)
    timerstart = REAL(count, SP)/REAL(count_rate, SP)   ! Initialize timer

!    CALL symbolic_fill

    CALL SYSTEM_CLOCK(count,count_rate)
    timer1 = REAL(count, SP)/REAL(count_rate, SP) ! Wall clock for mesh input and sorting.

!    CALL FILL_CONVERT  

    CALL SYSTEM_CLOCK(count,count_rate)
    timer2 = REAL(count, SP)/REAL(count_rate, SP) ! Wall clock for mesh input and sorting.

    CALL SYSTEM_CLOCK(count,count_rate,count_max)
    timer2b = REAL(count, SP)/REAL(count_rate, SP)  ! Initialize timer

!    CALL DOOLITTLE_LU  

    CALL SYSTEM_CLOCK(count,count_rate)
    timer3 = REAL(count, SP)/REAL(count_rate, SP) ! Wall clock for mesh input and sorting.

!    CALL CCStoCRS

    CALL SYSTEM_CLOCK(count,count_rate)
    timer4 = REAL(count, SP)/REAL(count_rate, SP) ! Wall clock for mesh input and sorting.


  howmny = 'A'
  
  iparam = 0
  iparam(1) = 1      ! exact shift is used
  iparam(3) = maxitr ! maximum number of iterations
  iparam(4) = 1      ! this must be 1 
  iparam(7) = mode    ! this is the mode

  
  tol = residual_norm
  allocate(select(ncv))
  
  allocate(resid(dof))
  resid = 0_SP
  
  allocate(iwork(dof))
  iwork = 0
  allocate(d(nev))
  d = 0
  
  lworkl = ncv*ncv + 8*ncv
  allocate(workl(lworkl))
  allocate(workd(dof*3))
  workd = 0_SP
  workl = 0_SP
  allocate(z(dof,nev))
  z = 0_SP
  allocate(v(dof,ncv))
  v = 0_SP
  info = 0          

  PRINT*, 'howmny ',  howmny 
  !PRINT*, 'select ',  SELECT 
  !PRINT*, 'd      ',  d      
  !PRINT*, 'z      ',  z      
  PRINT*, 'dof    ',  dof    
  PRINT*, 'nev    ',  nev    
  !PRINT*, 'resid  ',  resid  
  PRINT*, 'ncv    ',  ncv    
  !PRINT*, 'v      ',  v      
  PRINT*, 'dof    ',  dof    
  PRINT*, 'iparam ',  iparam 
  !PRINT*, 'workd  ',  workd  
  !PRINT*, 'workl  ',  workl  
  PRINT*, 'lworkl ',  lworkl 
  !PRINT*, 'iwork  ',  iwork  
  PRINT*, 'info   ',  info   

  CALL sparse_ARPACK(howmny, SELECT, d, z, dof, &
     nev, resid, ncv, v, dof, &
     iparam, workd, workl, lworkl, iwork, info)

  PRINT*, "after CALL sparse_ARPACK "
  CALL SYSTEM_CLOCK(count,count_rate)
  timer5 = REAL(count, SP)/REAL(count_rate, SP) ! Wall clock for mesh input and sorting.
  
  ! jps changed count/count_rate to the above to correct timing - integer division was incorrectly truncated
  ! jps changed format in write statements to allow time to be printed as a real
  WRITE (*, '(A)') '*****************Timing results*****************************'
  WRITE (*, '(A,I8)') 'dof                 ', dof
  WRITE (*, '(A,F8.6)') 'fill-in             ',timer1 - timerstart
  WRITE (*, '(A,F8.6)') 'convert fill-in     ',timer2 - timer1
  WRITE (*, '(A,F8.6)') 'Doolittle           ',timer3 - timer2b
  WRITE (*, '(A,F8.6)') 'convert CCS to CRS  ',timer4 - timer3
  WRITE (*, *) 'ARPACK eigenvalues  ',timer5 - timer4
  WRITE (*, *) 'total time (in s) req.', (timer3 - timer2b)+(timer5-timer4)+(timer4 - timer3)+(timer2 -timerstart)


END SUBROUTINE SPARSEIG_RI  

!**********************************************************************
! other routines used mostly by sparseig
!
!  find_all_in_list(listsize,list,element,poslist,dof)
!  fill_ones(listsize,fillsize,fillist,list)
!  find_greater_than(rowsize,rowind,elemk,nextrow)
!  convertfillcor(row,col)
!  
!**********************************************************************

SUBROUTINE find_all_in_list(listsize,list,element,poslist,dof)


USE unit_numbers
USE problem_info
USE nrtype

IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: listsize,dof
INTEGER(I4B),  INTENT(IN), DIMENSION(listsize):: list
INTEGER(I4B), INTENT(IN)  :: element
INTEGER(I4B),  INTENT(OUT), DIMENSION(dof) :: poslist
INTEGER(I4B) :: search 
INTEGER(I4B), DIMENSION(dof)  :: templist

INTEGER(I4B) pos
 

templist = 0  ! This is the default value
               ! and indicates that element is not found in list
 pos = 0                 
 DO search = 1,listsize
   IF (list(search).EQ.element) THEN
     pos = pos + 1
     templist(pos) = search
   END IF
 END DO
 
 poslist = templist

END SUBROUTINE find_all_in_list

!**********************************************************************

SUBROUTINE fill_ones(listsize,fillsize,fillist,list)

USE unit_numbers
USE problem_info
USE nrtype
IMPLICIT NONE
INTEGER(I4B),INTENT(IN) :: listsize, fillsize
INTEGER(I4B),DIMENSION(fillsize), INTENT(IN)  :: fillist
INTEGER(I4B), DIMENSION(listsize), INTENT(INOUT)  :: list
INTEGER(I4B) k



 DO k = 1,fillsize
    list(fillist(k)) = 1
 END DO
 
END SUBROUTINE fill_ones

!**********************************************************************
SUBROUTINE find_greater_than(rowsize,rowind,elemk,nextrow)

USE unit_numbers
USE problem_info
USE nrtype
IMPLICIT NONE

INTEGER(I4B),INTENT(IN) :: elemk, rowsize
INTEGER(I4B), DIMENSION(rowsize), INTENT(IN) :: rowind
INTEGER(I4B),  INTENT(OUT) :: nextrow
INTEGER(I4B) k

 ! this routine assumes the list is in ascending order
 nextrow = 0 
 DO k = 1,rowsize
   IF (rowind(k).GT.elemk) THEN
      nextrow = k
      EXIT
   END IF 
 END DO

END SUBROUTINE find_greater_than


!**********************************************************************
FUNCTION convertfillcor(row,col)
USE matrix
USE unit_numbers
USE problem_info
USE nrtype         
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: row,col
INTEGER(I4B) convertfillcor
INTEGER(I4B) :: rowstart, rowstop,rowsize, pos, search, row_pos
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: row_indices
 
  
  rowstart = urowind(row)
  rowstop  = urowind(row+1) -1
  rowsize = rowstop - rowstart + 1
  
  allocate(row_indices(rowsize))
  row_indices = ucolind(rowstart:rowstop)
  pos = 0  ! this indicates element not found in list...
  DO search = 1,rowsize
    if (row_indices(search).EQ.col) then
      row_pos = search
      pos = row_pos+rowstart-1
      EXIT  
    end if
  END DO
  deallocate(row_indices)
  convertfillcor = pos

END FUNCTION convertfillcor 

!**********************************************************************

FUNCTION convertcor(row,col,flag)
USE CBAA_data
USE matrix
USE unit_numbers
USE problem_info
USE nrtype         
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: row,col
INTEGER(I4B), INTENT(IN), OPTIONAL  :: flag  
INTEGER(I4B) convertcor
INTEGER(I4B) :: rowstart, rowstop,rowsize, pos, search, row_pos
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: row_indices
  
  IF ((.NOT.PRESENT(flag)).OR.(PRESENT(flag).AND.(flag.EQ.1))) THEN ! Use FE indices
    rowstart = row_ind(row)
    rowstop  = row_ind(row+1) -1
    rowsize = rowstop - rowstart + 1
    allocate(row_indices(rowsize))
    row_indices = col_ind(rowstart:rowstop)

  ELSE IF (flag.EQ.2) THEN ! Use BE/FMM indices
    rowstart = CBAA_BE_rowind(row)
    rowstop  = CBAA_BE_rowind(row+1) -1
    rowsize = rowstop - rowstart + 1
    allocate(row_indices(rowsize))
    row_indices = CBAA_BE_colind(rowstart:rowstop)
  END IF  

  pos = 0  ! this indicates element not found in list...
  DO search = 1,rowsize
    if (row_indices(search).EQ.col) then
      row_pos = search
      pos = row_pos+rowstart-1
      EXIT  
    end if
  END DO
  deallocate(row_indices)
  convertcor = pos

END FUNCTION convertcor 

!**********************************************************************

SUBROUTINE MATVECPROD_REAL(flag,x,y,n,sig)
  USE FEMINTERFACE, ONLY: count_nonzeros
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype         
  IMPLICIT NONE

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine computes various products (controlled by flag)
! of the real-valued sparse  matrix(es) S, T, A, B or C matrices with the input vector x, 
! producing output vector y. 
! Optional offset sig can also be provided. 
! This is the specific real implementation of the generic procedure
! MATVECPROD.
! Note that the implementation for flag > 5 differs between the REAL and COMPLEX
! implementations.
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
!  First version May 2000.  RHG.
!  July 2000.    Extended to include A. RHG&MMB.
!  20 July 2000. Specific name changed to MATVECPROD_REAL. sig made optional. DBD. 
!
!*************************************************************************
!  Input
!*************************************************************************
!
!  Ssparse and Tsparse;  or Asparse, Bsparse or Csparse; x. Optionally: sig
!  The sparse matrices are stored in Compressed Row Storage form.
!  The arrays rowind and colind provide respectively the starting points 
!  of each row and the column indices of each row. 
!  These, as well as the sparse matrices, are passed via the "matrix" module.
!
!*************************************************************************
!  Output
!*************************************************************************
!  y
!*************************************************************************
  INTEGER(I4B), INTENT(IN) :: flag, n
! Error corrected DBD 29 April - fortunately this actually worked, although it was wrong....
!  REAL(I4B), DIMENSION(n), INTENT(IN) :: x
!  REAL(I4B), DIMENSION(n), INTENT(OUT) :: y
  REAL(SP), DIMENSION(n), INTENT(IN) :: x
  REAL(SP), DIMENSION(n), INTENT(OUT) :: y

  REAL(SP), INTENT(IN), OPTIONAL  :: sig
  INTEGER(I4B), DIMENSION(n) :: colind
  INTEGER(I4B) :: colsize, row, rowsize, rowstart, rowstop

  IF ((flag.EQ.3.OR.flag.EQ.4).AND..NOT.PRESENT(sig)) THEN
    STOP 'IE: <sig> required in MATVECPROD_REAL.'
  END IF

 !write *,'entering matvecprod '
 !write *,'input(1:20) = ',x(1:20) 

 !remember that S and T share the same sparsity pattern

  !S*x calculated*****************************************

  if (flag.EQ.1) then
     do row = 1,n
    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y(row) = sum(Ssparse(rowstart:rowstop)*x(colind(1:rowsize)))

     end do   

  !T*x calculated*****************************************
   else if (flag.EQ.2) then
     do row = 1,n
    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y(row) = sum(Tsparse(rowstart:rowstop)*x(colind(1:rowsize))) 

     end do   

  !(S- sig*T) calculated********************************
  else if (flag.EQ.3) then
    do row = 1,n

    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y(row) = sum((Ssparse(rowstart:rowstop)- sig*Tsparse(rowstart:rowstop))*x(colind(1:rowsize))) 

     end do   

  !(S - sig*I) calculated*******************************
  else if (flag.EQ.4) then

    do row = 1,n

    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y(row) = sum((Ssparse(rowstart:rowstop)- sig)*x(colind(1:rowsize))) 

     end do   

  !A*x calculated***************************************
  else if (flag.EQ.5) then
    DO row = 1,n
      y(row) = SUM(Asparse(row_ind(row):row_ind(row+1)-1)* &
                  x(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO


! Added DBD 27 April 2003

  ! No equivalent of flag=6 (conjugate operation for complex).

  !B*x calculated***************************************
  else if (flag.EQ.7) then
    DO row = 1,n
      y(row) = SUM(Bsparse(row_ind(row):row_ind(row+1)-1)* &
                  x(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO

  !C*x calculated***************************************
  else if (flag.EQ.8) then
    DO row = 1,n
      y(row) = SUM(Csparse(row_ind(row):row_ind(row+1)-1)* &
                  x(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO

  else 
    STOP 'Internal error in MATVECPROD_REAL: called with invalid flag'
! End added DBD 27 April 2003
  end if
 
END SUBROUTINE MATVECPROD_REAL

!*************************************************************************

SUBROUTINE MATVECPROD_DP(flag,x_DP,y_DP,n,sig)
  USE FEMINTERFACE, ONLY: count_nonzeros
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype         
  IMPLICIT NONE

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine computes various products (controlled by flag)
! of the double-precision real-valued sparse  matrix(es) A, B or C matrices with
! the input vector x,  producing output vector y. 
! Optional offset sig can also be provided. 
! This is the specific real implementation of the generic procedure
! MATVECPROD.
! Only flag=5 implemented at present.
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! David B Davidson
!*************************************************************************
!  Last revised:
!*************************************************************************
!
!  First version 1 May 2003.  DBD.
!
!*************************************************************************
!  Input
!*************************************************************************
!
!  Asparse_DP, Bsparse or Csparse; x. Optionally: sig
!  The sparse matrices are stored in Compressed Row Storage form.
!  The arrays rowind and colind provide respectively the starting points 
!  of each row and the column indices of each row. 
!  These, as well as the sparse matrices, are passed via the "matrix" module.
!
!*************************************************************************
!  Output
!*************************************************************************
!  y
!*************************************************************************
  INTEGER(I4B), INTENT(IN) :: flag, n
  REAL(DP), DIMENSION(n), INTENT(IN) :: x_DP
  REAL(DP), DIMENSION(n), INTENT(OUT) :: y_DP

  REAL(DP), INTENT(IN), OPTIONAL  :: sig ! Not used in this routine, but retained 
                                         ! for compatibility with generic routine
  INTEGER(I4B), DIMENSION(n) :: colind
  INTEGER(I4B) :: colsize, row, rowsize, rowstart, rowstop

  IF ((flag.EQ.3.OR.flag.EQ.4).AND..NOT.PRESENT(sig)) THEN
    STOP 'IE: <sig> required in MATVECPROD_REAL.'
  END IF


  !A*x calculated***************************************
  
  IF (flag.EQ.5) THEN
    DO row = 1,n
      y_DP(row) = SUM(Asparse_DP(row_ind(row):row_ind(row+1)-1)* &
                  x_DP(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO
  ELSE 
    STOP 'Internal error in MATVECPROD_REAL: called with invalid flag'
  END IF
END SUBROUTINE MATVECPROD_DP

!**********************************************************************

SUBROUTINE MATVECPROD_COMPLEX(flag,x_c,y_c,n,sig)
  USE FEMINTERFACE, ONLY: count_nonzeros
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype         
  IMPLICIT NONE

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine computes various products (controlled by flag)
! of the complex-valued sparse  matrix(es) S, T or A with the input vector x, 
! producing output vector y. 
! Optional offset sig can also be provided. 
! This is the specific complex implementation of the generic procedure
! MATVECPROD
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
!  First version May 2000.  RHG.
!  July 2000.    Extended to include A. RHG&MMB.
!  20 July 2000. Specific name changed to MATVECPROD_COMPLEX. sig made optional. DBD. 
!
!*************************************************************************
!  Input
!*************************************************************************
!
!  Ssparse_c and Tsparse_c;  or Asparse_c; x_c. Optionally: sig.
!  Further, see documentation for MATVECPROD_REAL.
!
!*************************************************************************
!  Output
!*************************************************************************
!  y
!*************************************************************************


  INTEGER(I4B), INTENT(IN) :: flag, n
  COMPLEX(SPC), DIMENSION(n), INTENT(IN) :: x_c
  COMPLEX(SPC), DIMENSION(n), INTENT(OUT) :: y_c
  REAL(SP), INTENT(IN), OPTIONAL  :: sig
  INTEGER(I4B), DIMENSION(n) :: colind
  INTEGER(I4B) :: colsize, row, rowsize, rowstart, rowstop

  IF ((flag.EQ.3.OR.flag.EQ.4).AND..NOT.PRESENT(sig)) THEN
    STOP 'IE: <sig> required in MATVECPROD_COMPLEX.'
  END IF

  !remember that S, T and A share the same sparsity pattern

  !S*x calculated*****************************************

  if (flag.EQ.1) then
     do row = 1,n
    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y_c(row) = sum(Ssparse_c(rowstart:rowstop)*x_c(colind(1:rowsize)))

     end do   

  !T*x calculated*****************************************
   else if (flag.EQ.2) then
     do row = 1,n
    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y_c(row) = sum(Tsparse_c(rowstart:rowstop)*x_c(colind(1:rowsize))) 

     end do   

  !(S- sig*T) calculated********************************
  else if (flag.EQ.3) then
    do row = 1,n

    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y_c(row) = sum((Ssparse_c(rowstart:rowstop)- sig*Tsparse_c(rowstart:rowstop))*x_c(colind(1:rowsize))) 

     end do   

  !(S - sig*I) calculated*******************************
  else if (flag.EQ.4) then

    do row = 1,n

    rowstart = row_ind(row)
    rowstop = row_ind(row+1) -1
    rowsize = rowstop - rowstart +1
    colind(1:rowsize) = col_ind(rowstart:rowstop)
    y_c(row) = sum((Ssparse_c(rowstart:rowstop)- sig)*x_c(colind(1:rowsize))) 

     end do   
  !A*x calculated***************************************
  else if (flag.EQ.5) then
    DO row = 1,n
      y_c(row) = SUM(Asparse_c(row_ind(row):row_ind(row+1)-1)* &
                  x_c(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO
  !Conjg(A)*x calculated***************************************
  else if (flag.EQ.6) then
    DO row = 1,n
      y_c(row) = SUM(CONJG(Asparse_c(row_ind(row):row_ind(row+1)-1))* &
                  x_c(col_ind(row_ind(row):row_ind(row+1)-1)))
    END DO
! Added DBD 27 April 2003
  else 
    STOP 'Internal error in MATVECPROD_REAL: called with invalid flag'
! End added DBD 27 April 2003
  end if

END SUBROUTINE MATVECPROD_COMPLEX


SUBROUTINE lowersolve2(n,m,x)
!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This is the version of lowersolve that works with the L factor stored
! separately and in CRS format. 
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
! 15 August 2000 RHG
!
!*************************************************************************

USE FEMINTERFACE, ONLY: convertfillcor
USE matrix
USE unit_numbers
USE problem_info
USE nrtype         
IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: n
REAL(SP), DIMENSION(n), INTENT(IN) :: m
REAL(SP), DIMENSION(n), INTENT(OUT) :: x

INTEGER(I4B) :: rowstart, row, kkpos, rowstop, rowsize
REAL(SP), DIMENSION(n) :: Lks
  
  
  Lks = 0
  x = 0
  x(1) = m(1)/Lowerval(1)
  do row = 2,n
    rowstart = upper_rowind(row)
    rowstop = upper_rowind(row+1)-1 
    rowsize = rowstop - rowstart + 1
    if (rowsize.GT.1) then
      Lks(1:rowsize-1) = Lowerval(rowstart:rowstop-1)
  x(row) = 1/Lowerval(rowstop)*(m(row) - sum(Lks(1:rowsize-1)*x(upper_colind(rowstart:rowstop-1))))
    else
      x(row) = 1/Lowerval(rowstart)*m(row)
    end if
  end do

 
END SUBROUTINE lowersolve2



SUBROUTINE uppersolve2(n,x,y)
!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This is the version of uppersolve that works with the U factor stored
! separately and in CRS format. 
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
! 15 August 2000 RHG
!
!*************************************************************************


USE FEMINTERFACE, ONLY: convertfillcor
USE math_tools, ONLY: FIND_IN_LIST
USE matrix
USE unit_numbers
USE problem_info
USE nrtype         
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: n
REAL(SP), DIMENSION(n), INTENT(IN) :: x
REAL(SP), DIMENSION(n), INTENT(OUT) :: y

INTEGER(I4B) :: lastpos,rowstop, row, rowstart, rowsize
REAL(SP), DIMENSION(n) :: Uks

  
  Uks = 0_SP
  y = 0_SP
  lastpos = CRSupper_row(n+1) -1
  y(n) = x(n)/CRSval(lastpos)

  DO row = n-1,1,-1
    rowstart = CRSupper_row(row)
    rowstop = CRSupper_row(row+1)-1 
    rowsize = rowstop - rowstart + 1
    
    IF (rowsize.GT.1) THEN
      
      Uks(1:rowsize-1) = CRSval(rowstart+1:rowstop)
    
      y(row) = 1/CRSval(rowstart)*(x(row)-sum(Uks(1:rowsize-1)*y(CRSupper_col(rowstart+1:rowstop))))
    ELSE
      y(row) = 1/CRSval(rowstart)*x(row)
    END IF

   
  END DO  


END SUBROUTINE uppersolve2


SUBROUTINE CCStoCRS
USE matrix
USE unit_numbers
USE problem_info
USE nrtype        
IMPLICIT NONE

!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine is very specific, it converts an upper matrix 
! stored in CCS into upper matrix stored in CRS
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
! 15 August 2000 RHG
!
!*************************************************************************


 INTEGER(I4B) :: step, num, row, rowstart, rowstop, rowsize, kstep, d_loop, stat
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: row_entries, rowlist
 REAL(SP), DIMENSION(:), ALLOCATABLE :: rowvallist, dummy


 
 num = SIZE(upper_colind) 
 
 allocate(CRSupper_row(dof+1))
 CRSupper_row = 0
 
 allocate(CRSupper_col(num))
 CRSupper_col = 0
 allocate(row_entries(dof))
 row_entries = 0

 
 !calculate length of each row...
 
 do step = 1,num
    row_entries(upper_colind(step)) = row_entries(upper_colind(step)) + 1
 end do  
 
 CRSupper_row(1) = 1
  
 do kstep = 2,dof
   CRSupper_row(kstep) = sum(row_entries(1:kstep-1)) + 1
 end do
 CRSupper_row(dof+1) = sum(row_entries(1:dof)) + 1
 
 allocate(CRSval(num))
 CRSval = 0.0_SP
 deallocate(row_entries)
 allocate(row_entries(dof+1))
 row_entries = CRSupper_row

  do row = 1,dof
    rowstart = upper_rowind(row)
    rowstop = upper_rowind(row+1) -1
    rowsize = rowstop - rowstart + 1
    allocate(rowlist(1:rowsize))
    allocate(rowvallist(1:rowsize))
    rowlist = 0
    rowvallist = 0.0_SP
    rowlist(1:rowsize) = upper_colind(rowstart:rowstop)
    rowvallist(1:rowsize) = Upperval(rowstart:rowstop)
    CRSupper_col(row_entries(rowlist)) = row
    
    CRSval(row_entries(rowlist)) = rowvallist
    row_entries(rowlist) = row_entries(rowlist) + 1 
    deallocate(rowlist)
    deallocate(rowvallist) 
  end do
  

 deallocate(Upperval)
  
 deallocate(row_entries)


END SUBROUTINE CCStoCRS

SUBROUTINE DOOLITTLE_LU
!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! Calculates the L and U factors of the specified matrix, using the
! Doolittle algoritm.
! 
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  
!*************************************************************************
!  Last revised:
!*************************************************************************
!
! 15 August 2000 RHG
!
!*************************************************************************



USE FEMINTERFACE, ONLY: convertcor
USE matrix
USE unit_numbers
USE problem_info
USE nrtype        
IMPLICIT NONE

INTEGER(I4B) :: crowstart, crowstop, crowsize, rrowstart, rrowsize, rrowstop, apos, lupos
INTEGER(I4B) :: prevdiagpos, col, row, rowcount, d_loop, pos
! Changed DBD 26 March 2003
REAL(SP), DIMENSION(:), ALLOCATABLE :: cind,rind, dummy !, Csparse
! End changed DBD 26 March 2003
REAL(SP) :: sigma_D
! Changed DBD 26 March 2003
!REAL(SP), DIMENSION(:,:), ALLOCATABLE :: c_mat
! End changed DBD 26 March 2003

 WRITE (*, '(A)') 'LU decomposition in progress'
 allocate(Upperval(0:size(upper_colind)))
 allocate(Lowerval(0:size(upper_colind)))
 allocate(cind(dof))
 allocate(rind(dof)) 
 !initialise
 Upperval = 0_SP
 Lowerval = 0_SP
 IF (mode.EQ.3) THEN
   sigma2 = 1_SP
   sigma_D = sigma
 END  IF
 IF (mode.EQ.2) THEN
   sigma2 = 0_SP
   sigma_D = -1_SP
 END IF
 Upperval(1) =  sigma2*Ssparse(1) - sigma_D*Tsparse(1) 
 Lowerval(1) =  1.0_SP
 do col = 2,dof
   crowstart = upper_rowind(col)
   crowstop = upper_rowind(col+1) -1
   crowsize = crowstop - crowstart +1
   do rowcount = 1,crowsize 
     row = upper_colind(crowstart+rowcount-1)
     rrowstart = upper_rowind(row)
     rrowstop =  upper_rowind(row+1)-1
     rrowsize =  rrowstop - rrowstart 
     apos = CONVERTCOR(row,col)
     lupos = crowstart + rowcount -1
     prevdiagpos = upper_rowind(row+1) -1   
     if ((crowsize.GT.1).AND.(row.GT.1)) then
       cind = 0
       rind = 0
       rind(upper_colind(rrowstart:rrowstop)) = Lowerval(rrowstart:rrowstop)
       cind(upper_colind(crowstart:rowcount+crowstart-2)) = Upperval(crowstart:crowstart+rowcount-2)
       Upperval(lupos) = sigma2*Ssparse(apos) - sigma_D*Tsparse(apos)  - sum(rind*cind)
       cind = 0
       rind = 0
       rind(upper_colind(rrowstart:rrowstop)) = Upperval(rrowstart:rrowstop)
       cind(upper_colind(crowstart:rowcount+crowstart-2)) = Lowerval(crowstart:crowstart+rowcount-2) 
       Lowerval(lupos) = (sigma2*Ssparse(apos) - sigma_D*Tsparse(apos)  - sum(rind*cind))/Upperval(prevdiagpos)
     else
       Upperval(lupos) = (sigma2*Ssparse(apos) - sigma_D*Tsparse(apos))
       Lowerval(lupos) = (sigma2*Ssparse(apos) - sigma_D*Tsparse(apos))/Upperval(prevdiagpos)
     end if
    end do
 end do 
 deallocate(cind)
 deallocate(rind)
 WRITE (*, '(A)') ' DOOLITTLE LU factorisation completed '  
END SUBROUTINE DOOLITTLE_LU

SUBROUTINE SYMBOLIC_FILL
USE FEMINTERFACE, ONLY: find_all_in_list, &
                        count_nonzeros, fill_ones
USE matrix
USE math_tools, ONLY: FIND_IN_LIST
USE unit_numbers
USE problem_info
USE nrtype  
USE bandwidth       
IMPLICIT NONE


INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::  filledpos, tempucolind, dummy
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: tempfill, fillist 
INTEGER(I4B), DIMENSION(dof) :: storelist, rowpos
INTEGER(I4B) :: startsize, therow, nextrow, rowsize, diagpos, rowstop, rowpossize, k
INTEGER(I4B) :: therowstart, therowstop, therowsize, ustart, row, insert
INTEGER(I4B) :: elemk, search,pos,  storesize, rowstart, listsize, nnz
INTEGER(I4B) :: filledpossize, oldindsize, newindsize, newstartsize
INTEGER(I4B) :: firstrowstart, firstrowstop, firstrowsize, too_small_counter,d_loop
REAL(SP), DIMENSION(:), ALLOCATABLE :: x, y1, y2 
INTEGER(I4B) :: i, zerocount
 
 !**************************************************************
 ! 
 ! calculate the position of fill-ins and set up indices 
 !
 !**************************************************************
zerocount = 0 
  
  ! assume here that about 20 percent of total possible entries 
  ! will be necessary 

  WRITE (*, '(A)') 'start of fill '
  startsize = ceiling(0.20*dof*dof)
  IF (startsize.LE.dof) THEN
     startsize = dof*2
  END IF 
  too_small_counter = 0  ! this will count the number of times the column index
                         ! was enlarged
  allocate(ucolind(startsize))
  allocate(urowind(dof+1))
  allocate(tempfill(dof))
  allocate(filledpos(dof))
  ucolind = 0
  urowind = 0
  filledpos = 0  
  tempfill = 0
  firstrowstart = row_ind(1)
  firstrowstop = row_ind(2) -1
  firstrowsize = firstrowstop - firstrowstart + 1
  ucolind(1:firstrowsize) = col_ind(firstrowstart:firstrowstop)
  urowind = row_ind(1:2)

  ustart = firstrowstop + 1
  DO row = 2,dof
     IF (ustart.GE.(startsize-dof)) THEN
       
       newstartsize = startsize + dof + startsize*0.05
       !probably we will run out of space soon...
       too_small_counter = too_small_counter + 1
       WRITE (*, '(A,I8)') 'too_small_counter = ', too_small_counter
       allocate(tempucolind(startsize))
       tempucolind = 0
       tempucolind(1:startsize) = ucolind(1:startsize)
       deallocate(ucolind)
       allocate(ucolind(newstartsize))
       ucolind = 0  ! initialisation
       ucolind(1:startsize) = tempucolind(1:startsize)
       deallocate(tempucolind)
       startsize = newstartsize 
     END IF
     tempfill = 0
     rowstart = row_ind(row)
     rowstop  = row_ind(row+1) -1
     rowsize  = rowstop - rowstart + 1
     rowpos = 0 
    
     tempfill(col_ind(rowstart:rowstop)) = 1
    
     CALL find_all_in_list(urowind(row)-1,ucolind,row,rowpos,dof)
     CALL COUNT_NONZEROS(dof,rowpos,rowpossize)
     DO k = 1,rowpossize
       elemk = rowpos(k)
       CALL find_greater_than(ustart,urowind,elemk,nextrow) 
       IF (nextrow.GT.0) then
          therow = nextrow -1
          therowstart = urowind(therow)
          therowstop  = urowind(therow+1) - 1
          therowsize  = therowstop - therowstart
          CALL FIND_IN_LIST(ucolind(therowstart:therowstop),therowsize, therow,diagpos)


!!$             tempfill(ucolind(therowstart+diagpos-1:therowstop)) = & 
!!$             (ucolind(therowstart+diagpos-1:therowstop)-10)/(ucolind(therowstart+diagpos-1:therowstop)-10)
! JPS: removed -10 from ucolind index to prevent divide by zero error. See original statement above

             tempfill(ucolind(therowstart+diagpos-1:therowstop)) = & 
             (ucolind(therowstart+diagpos-1:therowstop))/(ucolind(therowstart+diagpos-1:therowstop))
 100      END IF 
     END DO
     storelist = 0
     pos = 0
     DO search = 1,dof
       IF (tempfill(search).GT.0) then
         pos = pos + 1
         storelist(pos) = search
       END IF
     END DO 
     storesize = pos
     ucolind(ustart:ustart+storesize-1) = storelist(1:storesize)
     ustart = ustart + storesize
     urowind(row+1) = ustart
    
  END DO
  WRITE(*,*) 'zerocount is ', zerocount
  CALL count_nonzeros(size(ucolind),ucolind,nnz)
  urowind(dof+1) = nnz + 1
  WRITE (*, '(A)') ' symbolic fill completed'  !temporary debug data  
   
  deallocate(tempfill)
  deallocate(filledpos)   
  newindsize = ustart
  allocate(tempucolind(nnz))
  tempucolind(1:nnz) = ucolind(1:nnz)
  deallocate(ucolind)
  allocate(ucolind(ustart))
  ucolind = 0
  ucolind(1:nnz) = tempucolind(1:nnz)  
  deallocate(tempucolind) 
  CALL COUNT_NONZEROS(size(col_ind),col_ind,oldindsize) 
  

  
  WRITE (*, '(A)') 'fill-in calculated'
  WRITE (*, '(A,I8)') 'entries in T matrix : ', oldindsize
  WRITE (*, '(A,I8)') 'entries in LU matrix : ', newindsize
  WRITE (*, '(A,I8)') 'percentage full - T matrix  : ', 100*oldindsize/dof/dof
  WRITE (*, '(A,I8)') 'percentage full -LU matrix  : ', 100*newindsize/dof/dof
  

  !write *,'ucolind = ', ucolind
  !write *,'urowind = ', urowind
END SUBROUTINE SYMBOLIC_FILL
!**********************************************************************


SUBROUTINE FILL_CONVERT
!**********************************************************************
!
!Discards the lower part of the fill-in indices.
!
!**********************************************************************
USE FEMINTERFACE, ONLY: find_all_in_list
USE matrix
USE math_tools, ONLY: FIND_IN_LIST
USE unit_numbers
USE problem_info
USE nrtype  
USE bandwidth       
IMPLICIT NONE

INTEGER(I4B) :: row, diagpos, rowsize, rowstart, rowstop, nextstart, newnextstart
INTEGER(I4B) :: ucolindsize, redu_size, d_loop
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: dummy

! This routine converts the fill-in data to fill-in indices for the
! upper matrix  only

  
  CALL COUNT_NONZEROS(size(ucolind),ucolind,ucolindsize)    
  ALLOCATE(upper_colind(((ucolindsize-dof)/2) + dof)) 
  ALLOCATE(upper_rowind(dof+1))
  upper_colind = 0
  upper_rowind = 0 
  nextstart = 1
  newnextstart = 0

  DO row = 1, dof
    rowstart = urowind(row)
    rowstop = urowind(row+1) -1
    rowsize = rowstop - rowstart + 1    
    CALL FIND_IN_LIST(ucolind(rowstart:rowstop),rowsize,row,diagpos)
    redu_size = diagpos
   
    newnextstart = nextstart + redu_size
    upper_colind(nextstart:newnextstart-1) = ucolind(rowstart:diagpos+rowstart-1)
    upper_rowind(row) = nextstart
   
    nextstart = newnextstart 

  END DO   
  upper_rowind(dof+1) = nextstart 
  
  
  deallocate(ucolind)
  deallocate(urowind)
  
  WRITE (*, '(A)') 'end of fill_convert'

END SUBROUTINE fill_convert
!*******************************************************************************


SUBROUTINE symtest(passedtest)
!*******************************************************************************
! The routine checks that the sparse S matrix is symmetrical (in which case 
! the t matrix is also correct)
!*******************************************************************************
  USE nrtype
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE unit_numbers
  USE feminterface, ONLY: convertcor
  IMPLICIT NONE 
  
  LOGICAL(LGT) , INTENT(OUT) :: passedtest
  INTEGER(I4B) testcol, testrow,pos1, pos2, indpos, indpos_sym, errorval
  
  errorval = 0 
  passedtest = .TRUE.  
  DO testrow = 1,dof
    DO testcol = 1,dof
    
    pos1 = convertcor(testrow,testcol)
    pos2 = convertcor(testcol, testrow)
    IF ((pos1.EQ.0).AND.(.NOT.(pos2.EQ.0))) THEN
     errorval = 1
  
    END IF
    IF ((pos2.EQ.0).AND.(.NOT.(pos1.EQ.0))) THEN
     errorval = 1
    END IF

    IF (.NOT.(Ssparse(pos1).EQ.(Ssparse(pos2)))) THEN
     
    END IF
    END DO
  END DO
  
  IF (errorval.EQ.1) THEN
    CALL ERROR_FEMFEKO(1,4104)
  END IF
END SUBROUTINE symtest

!*******************************************************************************
SUBROUTINE sortcolind(collistsize, collist,  sortedlistsize, sortedlist)
!*******************************************************************************
! This is a very specific non general purpose sorting routine
! intended to help set up the column indices for the S and
! T matrix
! INPUT  : a non sorted list interspersed with zeros
! OUTPUT : the same list sorted, with zeros removed
!*******************************************************************************
 
  USE nrtype
  USE matrix
  IMPLICIT NONE  

  INTEGER(I4B), INTENT(IN) :: collistsize
  INTEGER(I4B), DIMENSION(collistsize), INTENT(IN) :: collist
  INTEGER(I4B),  INTENT(OUT) :: sortedlistsize
  INTEGER(I4B), DIMENSION(collistsize), INTENT(OUT) :: sortedlist

  INTEGER(I4B), DIMENSION(collistsize) :: newlist
  INTEGER(I4B) :: search, pos, sortk, sortj, tempval

  ! first remove zeros before sorting...
  
  
  newlist = 0
  pos = 0
  
  do search = 1,collistsize
    if (.NOT.(collist(search).EQ.0)) then
      pos = pos + 1
      newlist(pos) = collist(search)
    end if
  end do

  sortedlistsize = pos

  ! now sort the list...

  do sortk = 1,sortedlistsize
   do sortj = sortk+1,sortedlistsize

      if (newlist(sortk).GT.(newlist(sortj))) then
        tempval = newlist(sortk)
        newlist(sortk) = newlist(sortj)
        newlist(sortj) = tempval
      end if

   end do
 end do
 sortedlist = 0
 sortedlist(1:sortedlistsize)= newlist(1:sortedlistsize)

END SUBROUTINE sortcolind


!*******************************************************************************
SUBROUTINE remove_zeros(listsize, list, newlist)
!*******************************************************************************

USE unit_numbers
USE problem_info
IMPLICIT NONE
 
INTEGER(I4B) :: listsize
INTEGER(I4B), DIMENSION(listsize), INTENT(IN)  :: list 
INTEGER(I4B), DIMENSION(listsize), INTENT(OUT) :: newlist

INTEGER(I4B) :: k_search,insert_i, pos
INTEGER(I4B), DIMENSION(listsize) :: indexlist


indexlist = 1 ! array assignment
DO k_search = 1,listsize
  IF (list(k_search).EQ.0) THEN 
     indexlist(k_search) = 0        
  END IF 
END DO

pos = 1
newlist = 0  ! initialisation

DO insert_i = 1,listsize
  IF (indexlist(insert_i).EQ.1) THEN
     newlist(pos) = list(insert_i)
     pos = pos + 1
  END IF 
END DO

END SUBROUTINE remove_zeros


!*******************************************************************************
SUBROUTINE enlarge_sparse_matrices(estim_nnz,rowstart)
!*******************************************************************************
USE matrix
IMPLICIT NONE

INTEGER(I4B), INTENT(INOUT) :: estim_nnz
INTEGER(I4B), INTENT(IN)  :: rowstart 
INTEGER(I4B) :: prev_estim
INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: tempind 

if (rowstart.GE.estim_nnz-dof) then
  WRITE (*,'(A,I8)') 'increasing size'
  prev_estim = estim_nnz
  estim_nnz = estim_nnz + dof

  allocate(tempind(prev_estim))
  tempind = 0

  tempind(1:prev_estim) = col_ind(1:prev_estim)

  deallocate(col_ind)
  allocate(col_ind(estim_nnz))
  col_ind = 0  
  col_ind(1:prev_estim) = tempind(1:prev_estim)  

  deallocate(tempind) 
end if
 
END SUBROUTINE enlarge_sparse_matrices
!*******************************************************************************


SUBROUTINE MATRIX_SPARSE_ALLOCATE
  USE feminterface, ONLY : COUNT_NONZEROS, sortcolind, & 
                           enlarge_sparse_matrices, &
                           remove_zeros, remove_multiple
  USE geometry
  USE math_tools, ONLY: FIND_IN_LIST
  USE matrix 
  USE problem_info
  USE nrtype         
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
!
! The sparse matrices are allocated for storage.  These are:
! col_ind, row_ind, Ssparse and Tsparse or Asparse_c.
! 
! Additionally, the column and row indices are pre- prepared to indicate in which 
! matrix positions entries will exist.These were constructed in FAST_CONNECT_ELEMENT_TO_FACE 
! (only needed for LT/QN elements and higher). 
!
! Input required: "edges%free", established in routine NUMBER_DOF.
!
! Note: THIS ROUTINE ASSUMES ORDER 2 AND LOWER ELEMENTS.
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! R.H. Geschke.
!
!*******************************************************************************
! Originally written May 2000, as part of the eigenanalysis routine.
! Last revised: 27 June 2000, MMB & DBD. Generalized to CBAA and GW analysis.
! Revised again: 05 March 2002, DBD. Extended to include QT/QN elements. 
! Revised again: 25 March 2003, DBD. Extended to support time-domain analysis. 
! 
!*******************************************************************************
  INTEGER(I4B) rowstart
  INTEGER(I4B) rowpos, sortedsize, kedge, connedges, connfaces
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: sortedcolind, sortedlist, &
                                 conn_dofs, new_conn_dofs, tempcolind
  INTEGER(I4B) ielem, iedge, facepos, &
         edgenodes, pos, iface, &
         estim_nnz, dof_count, & 
         edgenum, facenum, elem_i, num_conn, newsize
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: element_dofs
  INTEGER(I4B) elem_dof

  ! Set up a list of dofs associated with every element <element_dofs>:
  CALL GET_ELEMENT_DOF ! internal subroutine

  elem_dof = ELEM_TET_MATRIX_SIZE ! Temporary working variable, to save space below.
                                  ! Max. number of dof's per element. 
  
  estim_nnz =  3*dof

  ALLOCATE(col_ind(estim_nnz))
  ALLOCATE(row_ind(dof+1))
  
  col_ind = 0
  row_ind = 0
  row_ind(1) = 1
  rowstart = 1 
 
  ! Cycle through all dofs, determine every dof's type and depending
  ! if it is an edge or face dof, compile a list of all the dofs that
  ! interacts with the dof under consideration. (A face dof interacts 
  ! with all the dofs associated with the 2 elements that share the 
  ! face. An edge dof interacts with all the dofs associated with the
  ! elements that share the edge.)
  DOF_LOOP: DO dof_count = 1,dof
      
    SELECT CASE(dof_type(dof_count))
      
    CASE(1)  ! this dof is edge type 1
      CALL FIND_IN_LIST(renumbered_e1,num_edges, dof_count,edgenum) 
      CALL COUNT_NONZEROS(SIZE(edgeconnectelem(edgenum,:)),edgeconnectelem(edgenum,:),num_conn)
      ALLOCATE(conn_dofs(elem_dof*num_conn))
      conn_dofs = 0
      DO elem_i = 1,num_conn
        conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(edgeconnectelem(edgenum,elem_i),1:elem_dof)
      END DO
       
    CASE(2) ! this dof is edge type 2
      CALL FIND_IN_LIST(renumbered_e2,num_edges, dof_count,edgenum) 
      CALL COUNT_NONZEROS(SIZE(edgeconnectelem(edgenum,:)),edgeconnectelem(edgenum,:),num_conn)
      allocate(conn_dofs(elem_dof*num_conn))
      conn_dofs = 0
      DO elem_i = 1,num_conn
        conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(edgeconnectelem(edgenum,elem_i),1:elem_dof)
      END DO
     
    CASE(3) ! this dof is face type 1
      CALL FIND_IN_LIST(renumbered_f1,num_faces, dof_count,facenum) 
      allocate(conn_dofs(2*elem_dof))
      conn_dofs = 0
      DO elem_i = 1,2 ! A face can be shared by at most one other element.
        IF (faceconnectelem(facenum,elem_i).GT.0) THEN  
          conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(faceconnectelem(facenum,elem_i),1:elem_dof)
        END IF 
      END DO
        
    CASE(4) ! this dof is face type 2
      CALL FIND_IN_LIST(renumbered_f2,num_faces, dof_count,facenum) 
      allocate(conn_dofs(2*elem_dof))
      conn_dofs = 0
      DO elem_i = 1,2 ! Ditto
        IF (faceconnectelem(facenum,elem_i).GT.0) THEN  
          conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(faceconnectelem(facenum,elem_i),1:elem_dof)
        END IF 
      END DO

	CASE(5) ! this dof is edge type 3
      CALL FIND_IN_LIST(renumbered_e3,num_edges, dof_count,edgenum) 
      CALL COUNT_NONZEROS(SIZE(edgeconnectelem(edgenum,:)),edgeconnectelem(edgenum,:),num_conn)
      allocate(conn_dofs(elem_dof*num_conn))
      conn_dofs = 0
      DO elem_i = 1,num_conn
        conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(edgeconnectelem(edgenum,elem_i),1:elem_dof)
      END DO

    CASE(6) ! this dof is face type 3
      CALL FIND_IN_LIST(renumbered_f3,num_faces, dof_count,facenum) 
      allocate(conn_dofs(2*elem_dof))
      conn_dofs = 0
      DO elem_i = 1,2 ! Ditto case(4)
        IF (faceconnectelem(facenum,elem_i).GT.0) THEN  
          conn_dofs(elem_i*elem_dof -(elem_dof-1):elem_i*elem_dof) = element_dofs(faceconnectelem(facenum,elem_i),1:elem_dof)
        END IF 
      END DO
        
    CASE DEFAULT
      STOP 'IE: Invalid dof_type in MATRIX_SPARSE_ALLOCATE.'
    END SELECT
      
    allocate(sortedcolind(SIZE(conn_dofs)))
    allocate(new_conn_dofs(SIZE(conn_dofs)))
    sortedcolind = 0
    new_conn_dofs = 0
    CALL remove_zeros(SIZE(conn_dofs), conn_dofs, new_conn_dofs)
      
    CALL sortcolind(SIZE(new_conn_dofs), new_conn_dofs, sortedsize, sortedcolind) 
    deallocate(new_conn_dofs)
    CALL remove_multiple(sortedsize, sortedcolind, sortedcolind)
    CALL COUNT_NONZEROS(sortedsize,sortedcolind,newsize)
      
    col_ind(rowstart:rowstart+newsize-1) = sortedcolind(1:newsize)
 
    deallocate(conn_dofs)
    deallocate(sortedcolind)
    rowstart = rowstart+newsize
    row_ind(dof_count+1) = rowstart
    IF (dof_count.LT.dof) THEN
      CALL enlarge_sparse_matrices(estim_nnz,rowstart)
    END IF

  END DO DOF_LOOP
 
  allocate(tempcolind(estim_nnz))
  tempcolind = 0
  tempcolind = col_ind

  CALL COUNT_NONZEROS(estim_nnz,col_ind,estim_nnz)
  
  deallocate(col_ind)  

  allocate(col_ind(estim_nnz))

  col_ind = 0
  col_ind(1:estim_nnz) = tempcolind(1:estim_nnz)

  deallocate(tempcolind)
  IF (CBAA_ANALYSIS.OR.GW_ANALYSIS.OR.FD_SCAT_ANALYSIS) THEN
    ALLOCATE(Asparse_c(estim_nnz))
    Asparse_c = (0.0,0.0)
  ELSE IF (REAL_EIGEN_ANALYSIS) THEN
    ALLOCATE(Ssparse(0:estim_nnz))
    ALLOCATE(Tsparse(0:estim_nnz))
    Ssparse = 0
    Tsparse = 0
! Added DBD 25 March 2003
  ELSE IF(TD_ANALYSIS) THEN
    ALLOCATE(Asparse(estim_nnz))
    ALLOCATE(Bsparse(estim_nnz))
    ALLOCATE(Csparse(estim_nnz))
    Asparse = (0.0,0.0)
    Bsparse = (0.0,0.0)
    Csparse = (0.0,0.0)
! End added DBD 25 March 2003
  ELSE 
    STOP 'IE: Invalid analysis type in routine MATRIX_SPARSE_ALLOCATE.'
  END IF

  DEALLOCATE(element_dofs)  ! allocated in GET_ELEMENT_DOF

CONTAINS

  SUBROUTINE GET_ELEMENT_DOF
    USE feminterface, ONLY: remove_zeros
    IMPLICIT NONE
    INTEGER(I4B), DIMENSION(6) :: element_edges, edgetype1_dof, & 
	                              edgetype2_dof, edgetype3_dof
    INTEGER(I4B), DIMENSION(4) :: element_faces, facetype1_dof, & 
	                              facetype2_dof, facetype3_dof
    INTEGER(I4B), DIMENSION(ELEM_TET_MATRIX_SIZE) :: temp_dofs
    INTEGER(I4B) :: elem_i

    ALLOCATE(element_dofs(num_elements,ELEM_TET_MATRIX_SIZE))  
	                     ! This is only used for LT/QN and QT/QN case,
                         ! CT/LN is very simple 
    element_dofs = 0
    DO elem_i = 1,num_elements
      element_edges(1:6) = elements(elem_i)%edges(1:6)
      edgetype1_dof(1:6) = renumbered_e1(element_edges)
      edgetype2_dof(1:6) = renumbered_e2(element_edges)
      edgetype3_dof(1:6) = renumbered_e3(element_edges)
      element_faces = elements(elem_i)%faces(1:4)
      facetype1_dof(1:4) = renumbered_f1(element_faces)   
      facetype2_dof(1:4) = renumbered_f2(element_faces)    
      facetype3_dof(1:4) = renumbered_f3(element_faces)    
      temp_dofs(1:6) = edgetype1_dof
      temp_dofs(7:12) =  edgetype2_dof
      temp_dofs(13:16) =  facetype1_dof 
      temp_dofs(17:20) =  facetype2_dof
      temp_dofs(21:26) =  edgetype3_dof
      temp_dofs(27:30) =  facetype3_dof
      CALL remove_zeros(ELEM_TET_MATRIX_SIZE, temp_dofs, & 
	                    element_dofs(elem_i,1:ELEM_TET_MATRIX_SIZE)) 
    END DO
  END SUBROUTINE GET_ELEMENT_DOF

END SUBROUTINE MATRIX_SPARSE_ALLOCATE
!*******************************************************************************

SUBROUTINE ITER_SOLVE(time_taken)
  USE CBAA_data
  USE feminterface, ONLY: MATVECPROD 
  USE cbaa_sys, ONLY: MATVECPROD_FMM
  USE geometry
  USE math_tools, ONLY: TIME_DIFFERENCE
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! Iterative solvers - for both full and sparse complex matrices. 
! Solves [A][x] = [b].  
! Note that the matrix [A} is assumed to be SYMMETRIC. THIS IS NOT CHECKED!
! Presently implemented:
! - Conjugate gradient 
! - Biconjugate gradient
!
! References: CG: Algorithm version 3. J-M Jin. "The Finite Element Method
! in Electromagnetics", (p400ff)
! Bi-CG: algorithm version 1, [ibid.],
! but read paragraph 11.2.2 too!!)  - or equivalently fig 9.18 (p328) Volakis.
! 
! Complex matrix A is stored in A_mat_c; 
! RHS vector b is stroed in b_vec_c;
! LHS vector x is x_vec_c.
! NB! All must have been allocated before this routine is called. This is tested.
!
! Note that the inner product in the Bi-CG and CG cases is NOT the same. 
! The conventional inner product in the CG is computed by dot-producting
! the vector with itself (the F90 DOT_PRODUCT function conjugates the first
! argument - strictly speaking, this yields the conjugate of the innner product 
! as in [Jin,eqn.11.48,p.396], but when the two arguments are the same as is
! always the case here, ! the result is real-valued and identical).
! When the Euclidean norm is required, for (numerical) safety, 
! the ABSolute value of this is computed 
! and the square root of this now real-valued number is then taken.
!
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! M M Botha and D B Davidson
!
!*******************************************************************************
! Initial version: 27 June 2000. 
! 
! Originally written: 
! May 2000, as part of the eigenanalysis routine. MMB.
! Last revised: 
! 27 June 2000. Generalized to CBAA and GW analysis. MMB and DBD.
! 12 July 2000. Error norm changed to usual L2 definition.  Test for converged
!               iteration added. 
! 20 July 2000. CG solver added. 
! 16 April 2001. Quasi-minimal residual (QMR) solver added. MMB.
! 17 April 2001. GMRES solver added. MMB.
!  3 Dec   2002. Experimental mixed potential preconditioner added. DBD.
!*******************************************************************************
! ERROR DIRECTORY
!*******************************************************************************
! 4102: IE: Internal error: unallocated variable referenced in ITER_SOLVE
!
!*******************************************************************************
  REAL(SP), INTENT(OUT) :: time_taken     ! time taken by this routine
! DBD start 11 Dec 02.
! Variables required by mixed potential preconditioner routines. Allocated in 
! internal routine MPPreCondCG if needed.
  INTEGER(I4B) :: num_free_edges
  COMPLEX(SPC),DIMENSION(:),ALLOCATABLE  :: pa,pv,qa,qv 
  INTEGER(I4B),DIMENSION(:),ALLOCATABLE  :: Gmin1,Gplus1
! End DBD 11 Dec 02
  ! Internal variable names (following Jin):
! MMB start 11 Apr, added temp_BE_vector:
  COMPLEX(SPC),DIMENSION(:),ALLOCATABLE  :: rk,pk,zk,Apk,Aark,temp_BE_vector
! MMB end 11 Apr
  COMPLEX(SPC) :: alphak,betak,gammak,rkzk1,rkzk2
  INTEGER(I4B) :: maxiterate,itcount,k


  REAL(SP) :: bnorm,error
  LOGICAL(LGT) success_termination ! Flag to check successful termination.
  INTEGER(I4B), DIMENSION(8) :: time_start,time_finish ! temp variables for timing all events
! MMB start 19 Apr, added temp_BE_vector:
  LOGICAL(LGT) :: preconditioning ! flag to indicate whether this option is enabled.
  COMPLEX(SPC),DIMENSION(:),ALLOCATABLE :: precond_sparse
! MMB end 19 Apr

! DBD start 6 Dec
! Following required for mixed-potential preconditioner. 
! Following line should eventually be moved to GEOMETRY for production code.
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: Node_edge_ptr, Node_edge_ind, Node_edge_ind_sign  
! DBD end 6 Dec

! MMB start 19 Apr, added temp_BE_vector and further change: 
  preconditioning = USE_PRE_CONDITIONER
  IF (preconditioning .AND. (.NOT.(REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS))) CALL SETUP_PRECONDITIONER
! MMB end 19 Apr

  ! Record starting time:
  CALL DATE_AND_TIME (values=time_start)
  ! problem lies here for eigenanalysis
  success_termination  = .FALSE. 
  IF((.NOT.SPARSE.AND..NOT.ALLOCATED(A_mat_c)).OR.& 
    (SPARSE.AND..NOT.ALLOCATED(Asparse_c)).OR.    &
    (.NOT.ALLOCATED(x_vec_c)).OR.                 & 
    (.NOT.ALLOCATED(b_vec_c))) THEN 
    CALL ERROR_FEMFEKO(1,4105)
  END IF 
  ALLOCATE(rk(dof))
  ALLOCATE(zk(dof))
  ALLOCATE(pk(dof))
  ALLOCATE(Apk(dof))
  ALLOCATE(Aark(dof)) ! Only for CG
! MMB start 11 Apr, added temp_BE_vector:
  ALLOCATE(temp_BE_vector(ap_dof))
! MMB end 11 Apr

  maxiterate = NINT(max_iter_factor*dof)
  ! Precompute the norm of b:
!  babs = ABS(SUM(b_vec_c*CONJG(b_vec_c)))
  bnorm = SQRT(ABS(DOT_PRODUCT(b_vec_c,b_vec_c)) )! Note that DOT_PRODUCT conjugates
                                                  ! first argument.
  ! Initialise the iterated variables:
  x_vec_c = (0.0,0.0) ! Arbitrary
  
  

  SOLVER_CHOICE: IF (SOLVER_TYPE.EQ.1) THEN ! Bi-CG
    ! This is [Jin,p.400] CG Algorithm version 1, with B and h replaced 
    ! by A and b respectively, and the inner product defined 
    ! as the SYMMETRIC inner product, eqn. (11.82), [Jin p.404], to yield
    ! the Bi-CG algorithm. 
    ! It assumes a complex symmetric matrix. 

    ! Initialize rk:  
    IF (SPARSE) THEN
      ! Sparse FEM matmult.:
      CALL MATVECPROD(5,x_vec_c,rk,dof)
      IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
        IF (CBAA_FMM_storage) THEN
          ! FMM MoM matmult.:
          CALL MATVECPROD_FMM(x_vec_c(1:ap_dof),temp_BE_vector)
          rk(1:ap_dof) = rk(1:ap_dof) + temp_BE_vector
        ELSE
        ! Full MoM matmult.:
          rk(1:ap_dof) = rk(1:ap_dof) + MATMUL(CBAA_BE_mat,x_vec_c(1:ap_dof))
        END IF
! MMB 11 Apr 2001 end
      END IF
      rk = b_vec_c - rk
    ELSE
      rk = b_vec_c - MATMUL(A_mat_c,x_vec_c)
    END IF
    IF (preconditioning) THEN 
      zk = precond_sparse*rk
    ELSE
      zk = rk
    END IF
    pk = zk ! NOT arbitrary!
    rkzk2 = SUM(rk*zk)

    ! Itererate
    BCG_iterate: DO itcount = 1,maxiterate
      rkzk1 = rkzk2
      ! Calculation of [A]{pk}:
      IF (SPARSE) THEN
        ! Sparse FEM matmult.:
        CALL MATVECPROD(5,pk,Apk,dof)
        IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
          IF (CBAA_FMM_storage) THEN
            ! FMM MoM matmult.:
            CALL MATVECPROD_FMM(pk(1:ap_dof),temp_BE_vector)
            Apk(1:ap_dof) = Apk(1:ap_dof) + temp_BE_vector
          ELSE
            ! Full MoM matmult.:
            Apk(1:ap_dof) = Apk(1:ap_dof) + MATMUL(CBAA_BE_mat,pk(1:ap_dof))
          END IF
! MMB 11 Apr 2001 end
        END IF
      ELSE
        Apk = MATMUL(A_mat_c,pk)
      END IF

      alphak = rkzk1/SUM(Apk*pk)
      x_vec_c = x_vec_c + alphak*pk
      rk = rk - alphak*Apk
      IF (preconditioning) THEN 
        zk = precond_sparse*rk
      ELSE
        zk = rk
      END IF
      rkzk2 = SUM(rk*zk)
      gammak = rkzk2/rkzk1
      pk = zk + gammak*pk
      error = SQRT(ABS(SUM(rk*CONJG(rk))))/bnorm
      IF (ON_SCREEN_REPORTING) THEN
      print *, itcount,'/',maxiterate, error
      END IF
      IF (error.LT.residual_norm) THEN
        success_termination = .TRUE.
        EXIT BCG_iterate
      END IF
    END DO BCG_iterate

   ELSE IF(SOLVER_TYPE.EQ.2) THEN ! CG
    ! This is [Jin,p.402] CG Algorithm version 3.
    ! It assumes a complex symmetric matrix. 
    ! The inner product is the conventional Euclidean norm.
    ! CAUTION - NOT TESTED FOR CBAA CASE!!

    ! Initialize rk:  
    IF (SPARSE) THEN
      ! Sparse FEM matmult.:
      CALL MATVECPROD(5,x_vec_c,rk,dof)
      IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
        IF (CBAA_FMM_storage) THEN
          ! FMM MoM matmult.:
          CALL MATVECPROD_FMM(x_vec_c(1:ap_dof),temp_BE_vector)
          rk(1:ap_dof) = rk(1:ap_dof) + temp_BE_vector
        ELSE
          ! Full MoM matmult.:
          rk(1:ap_dof) = rk(1:ap_dof) + MATMUL(CBAA_BE_mat,x_vec_c(1:ap_dof))
        END IF
! MMB 11 Apr 2001 end
      END IF
      rk = b_vec_c - rk
    ELSE
      rk = b_vec_c - MATMUL(A_mat_c,x_vec_c)
    END IF
    !Initialize pk:
    IF (SPARSE) THEN
      ! Sparse FEM matmult.:
      CALL MATVECPROD(6,rk,Aark,dof) ! Conjg(A)
      IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
        IF (CBAA_FMM_storage) THEN
          ! FMM MoM matmult.:
          CALL MATVECPROD_FMM(CONJG(rk(1:ap_dof)),temp_BE_vector)
          Aark(1:ap_dof) = Aark(1:ap_dof) + CONJG(temp_BE_vector)
        ELSE
          ! Full MoM matmult.:
          Aark(1:ap_dof) = Aark(1:ap_dof) + MATMUL(CONJG(CBAA_BE_mat),rk(1:ap_dof))
        END IF
! MMB 11 Apr 2001 end
      END IF
    ELSE
      Aark  = MATMUL(CONJG(A_mat_c),rk)
    END IF
    pk = Aark/DOT_PRODUCT(Aark,Aark) ! See note above re DOT_PRODUCT

    ! Iterate
    CG_iterate: DO itcount = 1,maxiterate
      ! Calculation of [A]{pk}:
      IF (SPARSE) THEN
    ! Sparse FEM matmult.:
    CALL MATVECPROD(5,pk,Apk,dof)
    IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
      IF (CBAA_FMM_storage) THEN
        ! FMM MoM matmult.:
        CALL MATVECPROD_FMM(pk(1:ap_dof),temp_BE_vector)
        Apk(1:ap_dof) = Apk(1:ap_dof) + temp_BE_vector
      ELSE
        ! Full MoM matmult.:
        Apk(1:ap_dof) = Apk(1:ap_dof) + MATMUL(CBAA_BE_mat,pk(1:ap_dof))
      END IF
! MMB 11 Apr 2001 end
    END IF
      ELSE
    Apk = MATMUL(A_mat_c,pk)
      END IF
      alphak = 1/DOT_PRODUCT(Apk,Apk)
      x_vec_c = x_vec_c + alphak*pk ! xk now overwritten with xk+1
      rk = rk - alphak*Apk          ! rk now overwritten with rk+1

      ! Calculation of [A^a][rk+1]:
      IF (SPARSE) THEN
    ! Sparse FEM matmult.:
    CALL MATVECPROD(6,rk,Aark,dof)
    IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
! MMB 11 Apr 2001 start, changed the action of this IF statement:
      IF (CBAA_FMM_storage) THEN
        ! FMM MoM matmult.:
        CALL MATVECPROD_FMM(CONJG(rk(1:ap_dof)),temp_BE_vector)
        Aark(1:ap_dof) = Aark(1:ap_dof) + CONJG(temp_BE_vector)
      ELSE
        ! Full MoM matmult.:
        Aark(1:ap_dof) = Aark(1:ap_dof) + MATMUL(CONJG(CBAA_BE_mat),rk(1:ap_dof))
      END IF
! MMB 11 Apr 2001 end
    END IF
      ELSE
    Aark = MATMUL(CONJG(A_mat_c),rk)
      END IF
      betak = 1/DOT_PRODUCT(Aark,Aark)
      pk = pk + betak*Aark          ! pk now overwritten with pk+1
      error = SQRT(ABS(DOT_product(rk,rk)))/bnorm
      IF (ON_SCREEN_REPORTING) THEN
       print *, itcount,'/',maxiterate, error 
      END IF
      IF (error.LT.residual_norm) THEN
    success_termination = .TRUE.
    EXIT CG_iterate
      END IF
    END DO CG_iterate

! begin MMB added 16 April 2001
  ELSE IF(SOLVER_TYPE.EQ.3) THEN ! QMR
    CALL QMR_SOLVE

  ELSE IF(SOLVER_TYPE.EQ.4) THEN ! GMRES
    CALL GMRES_SOLVE
! end MMB added 16 April 2001

   ELSE IF(SOLVER_TYPE.EQ.5) THEN ! CG
    CALL MPPreCondCG
  END IF SOLVER_CHOICE 

  ! Record finishing time:
  CALL DATE_AND_TIME (values=time_finish)
  time_taken = TIME_DIFFERENCE(time_start,time_finish)

  ! Write solution information to the output file:
  WRITE (FILEOUT,'(//,20X,A)')   'ITERATIVE SOLUTION OF THE SYSTEM MATRIX EQUATION'
  SELECT CASE (SOLVER_TYPE)
  CASE (1)
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Bi-conjugate gradient'
  CASE (2)
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Conjugate gradient'
! Added by MMB, 16 Apr 2001
  CASE (3)
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Quasi-Minimal Residual'
  CASE (4)
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Generalized Minimum Residual'
    WRITE (FILEOUT,'(1X,A,I4)')     'Iterations before restart:',restart_GMRES
  CASE (5)
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Mixed Potential Conjugate gradient'
  CASE DEFAULT
    STOP 'IE: Invalid iterative solver requested.'
! end MMB addition
  END SELECT
! Start DBD addition 25 Apr 2001
  IF(USE_PRE_CONDITIONER) THEN
    WRITE (FILEOUT,'(1X,A)')     'Diagonal preconditioning used.'
  ELSE
    WRITE (FILEOUT,'(1X,A)')     'No preconditioning used.'
  END IF     
! End DBD addition 25 Apr 2001
  WRITE (FILEOUT,'(1X,A,E9.3)')    'Specified residual norm:         ', residual_norm
  WRITE (FILEOUT,'(1X,A,I12)')     'Specified max. iterations:       ', maxiterate
  WRITE (FILEOUT,'(1X,A,E9.3)')    'Final residual norm:             ', error
  WRITE (FILEOUT,'(1X,A,I12)')     'Final num. iterations:           ', itcount
  WRITE (FILEOUT,'(1X,A,F12.3)')   'Time taken by iter. solve (sec): ', time_taken

  ! Deallocate local, allocatable variables:
  DEALLOCATE(rk)
  DEALLOCATE(zk)
  DEALLOCATE(pk)
  DEALLOCATE(Apk)
  DEALLOCATE(Aark)
! MMB start 11 Apr, added temp_BE_vector:
  DEALLOCATE(temp_BE_vector)
! MMB end 11 Apr
! MMB start 19 Apr, added temp_BE_vector:
  IF (preconditioning) DEALLOCATE(precond_sparse)
! MMB end 19 Apr

CONTAINS
  
! begin MMB 16 Apr, added
!*******************************************************************************
  SUBROUTINE QMR_SOLVE
!*******************************************************************************
! This is the Quasi-Minimal Residual Algorithm.
! References: See Barrett et al. 'Templates...', SIAM 1994.
! Chapter 9, Fig 9.6 in Volakis et al. 'FEM for EM:...'. The algorithm in
! Barrett can be rewritten as in Volakis for the symmetric case.
! A complex symmetric matrix is assumed. 
! Volakis et al. has a typo on line 7 of the iteration loop - must be theta_n = rho_{n+1}...
! MMB 16 April 2001
!*******************************************************************************
    ! QMR specific variables:
    COMPLEX(SPC),DIMENSION(:),ALLOCATABLE :: dk,sk,vk
    COMPLEX(SPC),DIMENSION(2) :: rho_k,theta_k,gamma_k
    COMPLEX(SPC) :: beta_k,delta_k,epsilon_k,eta_k

    ALLOCATE(dk(dof))
    ALLOCATE(sk(dof))
    ALLOCATE(vk(dof))

    ! Initialization:
    x_vec_c = (0.0,0.0) ! arbitrary
    ! Initialize rk:  
    IF (SPARSE) THEN
      ! Sparse FEM matmult.:
      CALL MATVECPROD(5,x_vec_c,rk,dof)
      IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
        IF (CBAA_FMM_storage) THEN
          ! FMM MoM matmult.:
          CALL MATVECPROD_FMM(x_vec_c(1:ap_dof),temp_BE_vector)
          rk(1:ap_dof) = rk(1:ap_dof) + temp_BE_vector
        ELSE
        ! Full MoM matmult.:
          rk(1:ap_dof) = rk(1:ap_dof) + MATMUL(CBAA_BE_mat,x_vec_c(1:ap_dof))
        END IF
      END IF
      rk = b_vec_c - rk
    ELSE
      rk = b_vec_c - MATMUL(A_mat_c,x_vec_c)
    END IF
    rho_k(2) = SQRT(SUM(rk*CONJG(rk)))
    vk = (1.0/rho_k(2))*rk
    pk = (0.0,0.0)
    gamma_k(2) = 1.0
    eta_k = -1.0
    epsilon_k = 1.0
    dk = (0.0,0.0)
    sk = (0.0,0.0)
    theta_k(2) = (0.0,0.0)

    ! Itererate:
    QMR_iterate: DO itcount = 1,maxiterate
      gamma_k(1) = gamma_k(2)
      rho_k(1)   = rho_k(2)
      theta_k(1) = theta_k(2)
      IF (rho_k(1).EQ.0.0) STOP 'IE: QMR algorithm crashed in ITER_SOLVE.'
      delta_k = SUM(vk*vk)
      IF (delta_k.EQ.0.0) STOP 'IE: QMR algorithm crashed in ITER_SOLVE.'
      pk = vk - (rho_k(1)*delta_k/epsilon_k)*pk
      ! Calculation of [A]{pk}:
      IF (SPARSE) THEN
        ! Sparse FEM matmult.:
        CALL MATVECPROD(5,pk,Apk,dof)
        IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
          IF (CBAA_FMM_storage) THEN
            ! FMM MoM matmult.:
            CALL MATVECPROD_FMM(pk(1:ap_dof),temp_BE_vector)
            Apk(1:ap_dof) = Apk(1:ap_dof) + temp_BE_vector
          ELSE
            ! Full MoM matmult.:
            Apk(1:ap_dof) = Apk(1:ap_dof) + MATMUL(CBAA_BE_mat,pk(1:ap_dof))
          END IF
        END IF
      ELSE
        Apk = MATMUL(A_mat_c,pk)
      END IF
      epsilon_k = SUM(pk*Apk)
      IF (epsilon_k.EQ.0.0) STOP 'IE: QMR algorithm crashed in ITER_SOLVE.'
      beta_k = epsilon_k/delta_k
      vk = Apk - beta_k*vk
      rho_k(2) = SQRT(SUM(vk*CONJG(vk)))
      vk = (1/rho_k(2))*vk
      theta_k(2) = rho_k(2)/(gamma_k(1)*ABS(beta_k))
      gamma_k(2) = 1/(SQRT(1+theta_k(2)**2))
      IF (gamma_k(2).EQ.0.0) STOP 'IE: QMR algorithm crashed in ITER_SOLVE.'
      eta_k = - eta_k*rho_k(1)*(gamma_k(2)**2)/(beta_k*(gamma_k(1)**2))
      dk = eta_k*pk + ((theta_k(1)*gamma_k(2))**2)*dk
      sk = eta_k*Apk + ((theta_k(1)*gamma_k(2))**2)*sk
      x_vec_c = x_vec_c + dk
      rk = rk - sk
      error = SQRT(ABS(SUM(rk*CONJG(rk))))/bnorm
      IF (ON_SCREEN_REPORTING) THEN
        print *, itcount,'/',maxiterate, error
      END IF
      IF (error.LT.residual_norm) THEN
        success_termination = .TRUE.
        EXIT QMR_iterate
      END IF
    END DO QMR_iterate

    DEALLOCATE(dk)
    DEALLOCATE(sk)
    DEALLOCATE(vk)

  END SUBROUTINE QMR_SOLVE
!*******************************************************************************

  SUBROUTINE SETUP_PRECONDITIONER
    USE feminterface, ONLY: CONVERTCOR
!*******************************************************************************
! Creates a sparsely stored preconditioner (actually its inverse). 
! Presently only diagonal. MMB. 19 Apr 2001
!*******************************************************************************
    INTEGER(I4B) :: precount,row

    ALLOCATE(precond_sparse(dof))
    DO precount = 1,dof
      IF (SPARSE) THEN
        row = CONVERTCOR(precount,precount)
        precond_sparse(precount) = Asparse_c(row)
        IF (CBAA_ANALYSIS.AND.(precount.LE.ap_dof)) THEN
          IF (CBAA_FMM_storage) THEN
            row = CONVERTCOR(precount,precount,2)
            precond_sparse(precount) = precond_sparse(precount) + CBAA_BE_val(row)
          ELSE
            precond_sparse(precount) = precond_sparse(precount) + CBAA_BE_mat(precount,precount)
          END IF
        END IF
      ELSE
        precond_sparse(precount) = A_mat_c(precount,precount)
      END IF
      precond_sparse(precount) = 1.0/precond_sparse(precount)
    END DO

  END SUBROUTINE SETUP_PRECONDITIONER
!*******************************************************************************

  SUBROUTINE GMRES_SOLVE
!*******************************************************************************
! This is the Generalised Minimal Residual Algorithm.
! References: See Barrett et al. 'Templates...', SIAM 1994.
! Chapter 9, Fig 9.7 in Volakis et al. 'FEM for EM:...'. 
! See also Saad, 'Iterative methods for sparse linear systems', PWS, 1996.
! In Saad the least squares solution is fully discussed, as it is implemented here.
! MMB 17 April 2001
!*******************************************************************************
    COMPLEX(SPC),DIMENSION(:),ALLOCATABLE :: wvec,svec,yvec,diag_precond
    COMPLEX(SPC),DIMENSION(:,:),ALLOCATABLE :: v_vectors,h_mat,rotations
    REAL(SP) :: rnorm
    COMPLEX(SPC) :: ctemp1,ctemp2
    INTEGER(I4B) :: max_vec,mcount1,mcount2

! DBD changed 27 Apr 2001
    max_vec = restart_GMRES ! Number of cycles before a restart
! End DBD changed 27 Apr 2001
    ALLOCATE(v_vectors(dof,max_vec))
    ALLOCATE(h_mat(max_vec+1,max_vec))
    ALLOCATE(wvec(dof))
    ALLOCATE(svec(max_vec+1))
    ALLOCATE(rotations(max_vec,2))
    ALLOCATE(yvec(max_vec))

    ! Initialization:
    x_vec_c = (0.0,0.0) ! arbitrary initial guess
    h_mat = (0.0,0.0)

    ! Iterate:
    GMRES_iterate: DO itcount = 1,maxiterate

      ! Initialize rk:  
      IF (SPARSE) THEN
        ! Sparse FEM matmult.:
        CALL MATVECPROD(5,x_vec_c,rk,dof)
        IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
          IF (CBAA_FMM_storage) THEN
            ! FMM MoM matmult.:
            CALL MATVECPROD_FMM(x_vec_c(1:ap_dof),temp_BE_vector)
            rk(1:ap_dof) = rk(1:ap_dof) + temp_BE_vector
          ELSE
          ! Full MoM matmult.:
            rk(1:ap_dof) = rk(1:ap_dof) + MATMUL(CBAA_BE_mat,x_vec_c(1:ap_dof))
          END IF
        END IF
        rk = b_vec_c - rk
      ELSE
        rk = b_vec_c - MATMUL(A_mat_c,x_vec_c)
      END IF
      IF (preconditioning) rk = precond_sparse*rk

      rnorm = SQRT(ABS(SUM(rk*CONJG(rk))))
      v_vectors(1:dof,1) = (1/rnorm)*rk
      svec = (0.0,0.0)
      svec(1) = rnorm

      DO mcount1 = 1,max_vec
        ! Calculation of wvec = [A]{v_mcount1}:
        IF (SPARSE) THEN
          ! Sparse FEM matmult.:
          CALL MATVECPROD(5,v_vectors(1:dof,mcount1),wvec,dof)
          IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
            IF (CBAA_FMM_storage) THEN
              ! FMM MoM matmult.:
              CALL MATVECPROD_FMM(v_vectors(1:ap_dof,mcount1),temp_BE_vector)
              wvec(1:ap_dof) = wvec(1:ap_dof) + temp_BE_vector
            ELSE
              ! Full MoM matmult.:
              wvec(1:ap_dof) = wvec(1:ap_dof) + MATMUL(CBAA_BE_mat,v_vectors(1:ap_dof,mcount1))
            END IF
          END IF
        ELSE
          wvec = MATMUL(A_mat_c,v_vectors(1:dof,mcount1))
        END IF
        IF (preconditioning) wvec = precond_sparse*wvec

        DO mcount2 = 1,mcount1
          h_mat(mcount2,mcount1) = SUM(CONJG(wvec)*v_vectors(1:dof,mcount2))
          wvec = wvec - h_mat(mcount2,mcount1)*v_vectors(1:dof,mcount2)
        END DO

        h_mat(mcount1+1,mcount1)   = SQRT(ABS(SUM(wvec*CONJG(wvec))))
        IF (mcount1.LT.max_vec) THEN
          v_vectors(1:dof,mcount1+1) = (1/h_mat(mcount1+1,mcount1))*wvec
        END IF
 
        ! Apply all existing rotations to the current column of <h_mat>,
        ! calculate the next pair of rotation constants and apply
        ! them as well (to h_mat and svec).
        ! -
        ! apply the old rotations to this column:
        DO mcount2 = 1,mcount1-1
          ctemp1 =   rotations(mcount2,1)*h_mat(mcount2,mcount1)   &
                   + rotations(mcount2,2)*h_mat(mcount2+1,mcount1)
          ctemp2 = - rotations(mcount2,2)*h_mat(mcount2,mcount1)   &
                   + rotations(mcount2,1)*h_mat(mcount2+1,mcount1)
          h_mat(mcount2,mcount1)   = ctemp1
          h_mat(mcount2+1,mcount1) = ctemp2
        END DO
        ! calculate the new rotations:
        rotations(mcount1,1) = &
          h_mat(mcount1,mcount1)/SQRT(h_mat(mcount1,mcount1)**2+h_mat(mcount1+1,mcount1)**2)
        rotations(mcount1,2) = &
          h_mat(mcount1+1,mcount1)/SQRT(h_mat(mcount1,mcount1)**2+h_mat(mcount1+1,mcount1)**2)
        ! apply the new rotations to this column:
        ctemp1 =   rotations(mcount1,1)*h_mat(mcount1,mcount1)   &
                 + rotations(mcount1,2)*h_mat(mcount1+1,mcount1)
        ctemp2 = - rotations(mcount1,2)*h_mat(mcount1,mcount1)   &
                 + rotations(mcount1,1)*h_mat(mcount1+1,mcount1)

        h_mat(mcount1,mcount1)   = ctemp1
        h_mat(mcount1+1,mcount1) = ctemp2
        ! apply the new rotations to <svec>:
        ctemp1 =   rotations(mcount1,1)*svec(mcount1) + rotations(mcount1,2)*svec(mcount1+1)
        ctemp2 = - rotations(mcount1,2)*svec(mcount1) + rotations(mcount1,1)*svec(mcount1+1)
        svec(mcount1)   = ctemp1
        svec(mcount1+1) = ctemp2

      END DO

      ! Calculate the V_vectors' coefficients (solve [H_m*m]{yvec} = {svec}):
      ! (H is upper triagular, so use back substitution.)
      DO mcount1 = max_vec,1,-1
        ctemp1 = (0.0,0.0)
        DO mcount2 = mcount1+1,max_vec
          ctemp1 = ctemp1 + h_mat(mcount1,mcount2)*yvec(mcount2)
        END DO
        yvec(mcount1) = (1/h_mat(mcount1,mcount1))*(svec(mcount1)-ctemp1)
      END DO

      ! Calculate the new, approximate answer:
      x_vec_c = x_vec_c + MATMUL(v_vectors,yvec)

!      error = ABS(svec(max_vec+1))/bnorm ! see page 163, eq.6.36 of Saad
      ! Calculate error the expensive way (but including the preconditioner!):
      IF (SPARSE) THEN
        ! Sparse FEM matmult.:
        CALL MATVECPROD(5,x_vec_c,rk,dof)
        IF (CBAA_ANALYSIS.AND.(ap_dof.GT.0)) THEN
          IF (CBAA_FMM_storage) THEN
            ! FMM MoM matmult.:
            CALL MATVECPROD_FMM(x_vec_c(1:ap_dof),temp_BE_vector)
            rk(1:ap_dof) = rk(1:ap_dof) + temp_BE_vector
          ELSE
          ! Full MoM matmult.:
            rk(1:ap_dof) = rk(1:ap_dof) + MATMUL(CBAA_BE_mat,x_vec_c(1:ap_dof))
          END IF
        END IF
        rk = b_vec_c - rk
      ELSE
        rk = b_vec_c - MATMUL(A_mat_c,x_vec_c)
      END IF
      error = SQRT(ABS(SUM(rk*CONJG(rk))))/bnorm


      IF (ON_SCREEN_REPORTING) THEN
        print *, itcount,'/',maxiterate, error
      END IF
      IF (error.LT.residual_norm) THEN
        success_termination = .TRUE.
        EXIT GMRES_iterate
      END IF
    END DO GMRES_iterate

    DEALLOCATE(v_vectors)
    DEALLOCATE(h_mat)
    DEALLOCATE(wvec)
    DEALLOCATE(svec)
    DEALLOCATE(rotations)
    DEALLOCATE(yvec)

  END SUBROUTINE GMRES_SOLVE
!*******************************************************************************
! end MMB 16 Apr, added

! Added DBD 4 Dec 02
  SUBROUTINE MPPreCondCG
    USE geometry
!*******************************************************************************
! This is the conjugate gradient algorithm using a mixed-potential preconditioner,
! based on theory in:
! Dyczij-Edlinger, Peng and Lee, "A fast vector-potential method using tangentially
! continuous vector finite elements", IEEE T-AP, Jun 1998, p.863-868.
! Note that the prescribed BC's (zero) on the nodal-based potential are handled explicitly by zero-ing
! these, rather than trying to remove them a priori. The reason is that it is possible
! for a prescribed node to be associated with a free edge (however, a prescribed edge 
! implies a prescribed node). Note also that ONLY an homogeneous (zero) Dirichlet 
! BC has been implemented at present.
!*******************************************************************************
    INTEGER(I4B) :: num_free_vertices
    INTEGER(I4B) :: iedge,inode,kk,free_edge 
    COMPLEX(SPC),DIMENSION(:),ALLOCATABLE  :: pkV,pkA,pkMP,rkV,rkA,temp_vec

    CALL COUNT_FREE_EDGES(num_free_edges)
    CALL FIND_AND_COUNT_FREE_VERTICES(num_free_vertices)

    IF(num_free_edges.NE.dof) THEN
      STOP 'Error in experimental code. Only implemented for CT/LN elements.'
    END IF
    
    ! Additional variables specific to this method
    ALLOCATE(pA(num_free_edges)) ! Input vector for matrix-vector product
    ALLOCATE(pV(num_nodes))      ! ditto
    ALLOCATE(qA(num_free_edges)) ! Output vector from matrix-vector product
    ALLOCATE(qV(num_nodes))      ! ditto
    ALLOCATE(Gmin1(num_free_edges))
    ALLOCATE(Gplus1(num_free_edges))

    ALLOCATE(pkA(num_free_edges))
    ALLOCATE(pkV(num_nodes))
    ALLOCATE(pkMP(num_free_edges+num_nodes))
    ALLOCATE(rkA(num_free_edges))
    ALLOCATE(rkV(num_nodes))
    ALLOCATE(temp_vec(num_free_edges+num_nodes))

    ! Initialize all these working arrays to zero.
    pA       = ZERO_C ! All arrays.
    pV       = ZERO_C
    qA       = ZERO_C
    qV       = ZERO_C
    pkA      = ZERO_C
    pkV      = ZERO_C
    pkMP     = ZERO_C
    rkA      = ZERO_C
    rkV      = ZERO_C
    temp_vec = ZERO_C
    Gmin1    = 0
    Gplus1    = 0

    ! Set up G indexes. Note that this datastructure couples FREE edges
	! and nodes. 
    free_edge = 0
    DO iedge = 1,num_edges
      IF (edges(iedge)%free) THEN
		free_edge = free_edge +1
        Gmin1(free_edge) = edges(iedge)%nodes(1)
        Gplus1(free_edge) = edges(iedge)%nodes(2)
	  END IF
    END DO
    IF (free_edge.NE.num_free_edges) THEN
	  STOP 'Internal error in MPPreCondCG. Inconsistent edge count.' ! Consistency check
	END IF

    ! This sparse G matrix cannot be simply transposed.
	! Instead, a sparse matrix with the necessary data is created.
    CALL NODEEDGE_INDEXLIST_MAKE

    IF(.NOT.GW_ANALYSIS) THEN
	  STOP 'Experimental code not tested for other applications.'
    END IF

    x_vec_c = (0.0,0.0) ! This is x_E = x_A in DEPL reference.

    ! Initialize rk: (r_E in )
    IF (SPARSE) THEN
      ! Sparse FEM matmult.:
      CALL MATVECPROD(5,x_vec_c,rkA,dof) ! Note rkA used as temporary storage here...
      rkA = b_vec_c - rkA                ! but now has its intended meaning. 
                                         !Note: rk (= rkE) = rkA 
    ELSE
	  STOP 'Not implemented'
    END IF
    pA(1:num_free_edges) = rkA(1:num_free_edges)
    pV = (0.0, 0.0) ! Not needed, but zeroed anyway.
    CALL MIXED_POTENTIAL_MATRIX_PRODUCT(num_free_edges,num_nodes,pA,pV,qA,qV,&
      adjoint_flag=.FALSE.,v_only_flag=.TRUE.)
    rkV(1:num_nodes) = qV(1:num_nodes)

    ! Initialize pk. Requires product of adjoint of A and r_k, k=1.
    pA(1:num_free_edges) = rkA(1:num_free_edges)
    pV(1:num_nodes) =rkV(1:num_nodes) 
    CALL MIXED_POTENTIAL_MATRIX_PRODUCT(num_free_edges,num_nodes,pA,pV,qA,qV,&
      adjoint_flag=.TRUE.,v_only_flag=.FALSE.)
    ! Build temporary vector (Adjoint of A times rk, k=1):
    temp_vec(1:num_free_edges) = qA(1:num_free_edges)
    temp_vec(num_free_edges+1:num_free_edges+num_nodes) =qV(1:num_nodes) 
    ! A^a r_1/<A^a r_1,A^a r_1>
    pkMP = temp_vec/DOT_PRODUCT(temp_vec,temp_vec) ! See note above re DOT_PRODUCT


!print *,'r1',rkA,rkV
!print *,'p1',pkMP
   
    ! Iterate
    MixedPot_CG_iterate: DO itcount = 1,maxiterate
      ! Find alpha_k = 1/<A p_k,A p_k>
      pkA = pkMP(1:num_free_edges)
      pkV = pkMP(num_free_edges+1:num_free_edges+num_nodes)
      ! Calculation of [A]{p_k}:
      pA(1:num_free_edges) = pkA(1:num_free_edges)
      pV(1:num_nodes)      = pkV(1:num_nodes) 
      CALL MIXED_POTENTIAL_MATRIX_PRODUCT(num_free_edges,num_nodes,pA,pV,qA,qV,&
        adjoint_flag=.FALSE.,v_only_flag=.FALSE.)
      temp_vec(1:num_free_edges) = qA(1:num_free_edges)
      temp_vec(num_free_edges+1:num_free_edges+num_nodes) =qV(1:num_nodes) 

!print *,'k=',itcount,'Apk',temp_vec

      alphak = 1/DOT_PRODUCT(temp_vec,temp_vec)

  	  ! Update solution vector and residual. Note that these are both for the E solution!
      x_vec_c = x_vec_c + alphak*pkA ! xk now overwritten with xk+1
      rkA = rkA - alphak*temp_vec(1:num_free_edges)  ! rk now overwritten with rk+1
      rkV = rkV - alphak*temp_vec(num_free_edges+1:num_free_edges+num_nodes)  

      ! Find beta_k = 1/<A^a r_k+1,A^a r_k+1>. Note r_k already overwrirten with r_k+1.
      ! Calculation of [A^a][rk+1]:
      pA(1:num_free_edges) = rkA(1:num_free_edges)
      pv(1:num_nodes)      = rkV(1:num_nodes) 
      CALL MIXED_POTENTIAL_MATRIX_PRODUCT(num_free_edges,num_nodes,pA,pV,qA,qV,&
        adjoint_flag=.TRUE.,v_only_flag=.FALSE.)
      temp_vec(1:num_free_edges) = qa(1:num_free_edges)
	  temp_vec(num_free_edges+1:num_free_edges+num_nodes) =qv(1:num_nodes) 
      betak = 1/DOT_PRODUCT(temp_vec,temp_vec)
      pkMP = pkMP + betak*temp_vec          ! pk now overwritten with pk+1

      error = SQRT(ABS(DOT_product(rkA,rkA)))/bnorm
      IF (ON_SCREEN_REPORTING) THEN
        print *, itcount,'/',maxiterate, error
!WRITE (FILEOUT,*) itcount,'/',maxiterate, error
      END IF

!WRITE (FILEOUT,*) 'Residual vector',rkA,rkV
!WRITE (FILEOUT,*) 'Direction vector',pkMP
       IF (error.LT.residual_norm) THEN
        success_termination = .TRUE.
        EXIT MixedPot_CG_iterate
      END IF
    END DO MixedPot_CG_iterate
    CALL NODEEDGE_INDEXLIST_CLEAN
    DEALLOCATE(pA)
    DEALLOCATE(pV)
    DEALLOCATE(qA)
    DEALLOCATE(qV)
    DEALLOCATE(Gmin1)
    DEALLOCATE(Gplus1)
    DEALLOCATE(pkA)
    DEALLOCATE(pkMP)
    DEALLOCATE(pkV)
    DEALLOCATE(rkA)
    DEALLOCATE(rkV)
    DEALLOCATE(temp_vec)

  END SUBROUTINE MPPreCondCG

  SUBROUTINE MIXED_POTENTIAL_MATRIX_PRODUCT(dim_vec,dim_scalar,pA,pV,qA,qV,adjoint_flag,v_only_flag)
    USE geometry
    IMPLICIT NONE
!*******************************************************************************
!*** Implement the matrix-vector product for the mixed-potential matrix-vector
!*** product.                                                               ****
!*******************************************************************************
! Uses the algorithm in Dyczij-Edlinger, Peng and Lee, 
! "A fast vector-potential method using tangentially
! continuous vector finite elements", IEEE T-AP, Jun 1998, p.863-868.
! Specifically, it implements equation 39 via eqns. 40,41, and 42.
! Note that the transpose assumes that the matrix is symmetric, and hence 
! only the complex conjugate is required. THIS IS NOT TESTED!
! Calling the routine without the optional flag implements the 
! matrix-vector product.
! Input:
! vector pA,pV (must be defined in calling routine)
! Output:
! vector qA (=pA if flag v_only is set, see below),qV.
! This is the product as in (39) with adjoint_flag false 
! or the product as in (39) but with the adjoint implemented (assuming a 
! symmetric matrix M) as the complex conjugate. 
! 
! Another option provided by this routine (with optional flag v_only_flag set to true)
! is to evaluate ONLY the product of G^T with the input A vector. 
! 
! Note that the scalar potential is assumed to be ZERO at prescribed nodes as appropriate for 
! PEC boundaries. Again, this is assumed, and NOT tested for!
!
! Author: DB Davidson, 11 Dec 2002.
!*******************************************************************************
    INTEGER (I4B), INTENT(IN) :: dim_vec,dim_scalar
    LOGICAL (LGT), OPTIONAL  :: adjoint_flag,v_only_flag
    INTEGER (I4B) inode,kk,free_edge ! counters.
    COMPLEX(SPC), DIMENSION(dim_vec), INTENT (IN) :: pA
    COMPLEX(SPC), DIMENSION(dim_scalar), INTENT (IN) :: pV
    COMPLEX(SPC), DIMENSION(dim_vec), INTENT (OUT) :: qA
    COMPLEX(SPC), DIMENSION(dim_scalar), INTENT (OUT) :: qV
    COMPLEX(SPC),DIMENSION(:),ALLOCATABLE  :: pe


    IF(.NOT.PRESENT(adjoint_flag)) adjoint_flag = .false.
    IF(.NOT.PRESENT(v_only_flag)) v_only_flag = .false.
    ALLOCATE(pe(dim_vec))

    IF (.NOT.v_only_flag) THEN
      DO free_edge =1,dim_vec
        pe(free_edge) = pA(free_edge) - pV(Gmin1(free_edge)) + pV(Gplus1(free_edge))  
      END DO
      IF (SPARSE) THEN
        ! Sparse FEM matmult.:
        IF (.NOT.adjoint_flag) THEN
          CALL MATVECPROD(5,pe,qA,dof) 
        ELSE
          CALL MATVECPROD(6,pe,qA,dof) ! Conjg. A.
        END IF
      END IF
    ELSE 
      ! Product of G with input vector only 
! Correction DBD 14 Dec.
      qA(1:dim_vec) = pA(1:dim_vec) 
    END IF

!print *,'Gmin1',Gmin1
!print *,'Gplus1',Gplus1
!print *,'pE',pE


    ! Complete matrix product for scalar component: 
    qV = (0.,0.)
    DO inode = 1,num_nodes
      DO kk = 1+Node_edge_ptr(inode), Node_edge_ptr(inode+1)
        IF(renumbered_e1(Node_edge_ind(kk)).NE.0) THEN ! Free edge
          qV(inode) = qV(inode) + qA(renumbered_e1(Node_edge_ind(kk))) *  Node_edge_ind_sign(kk)
!print *,'qv(',inode,')= ',qv(inode)
        END IF
      END DO
      IF(.NOT.vertices(inode)%free) THEN
        qV(inode) = (0.,0.) ! Zero value of potential and all related values at prescribed nodes.
!print *,'qv(',inode,') set to zero'
      END IF
    END DO
    ! Cleanup temporary variables

    DEALLOCATE(pe)

    END SUBROUTINE MIXED_POTENTIAL_MATRIX_PRODUCT


!*******************************************************************************
! Added DBD 4-11 Dec 02
  SUBROUTINE FIND_AND_COUNT_FREE_VERTICES(num_free_vertices)
    USE geometry
    IMPLICIT NONE
!*******************************************************************************
!*** Count free nodes.                                                       
!*******************************************************************************
! If the mixed-potential preconditioner is used, the free nodes 
! must be flagged. 
!*******************************************************************************
    INTEGER(I4B), INTENT(OUT) :: num_free_vertices
    INTEGER(I4B) :: iedge, inode ! counters

    ! Flag all vertices as free initially:
    vertices(1:num_nodes)%free = .TRUE.

    ! Each prescribed edge implies that the associated vertices are also
    ! precribed:

    ! Cycle edges:
    DO iedge = 1,num_edges 
      IF (.NOT.edges(iedge)%free) THEN 
        vertices(edges(iedge)%nodes(1))%free = .FALSE.
        vertices(edges(iedge)%nodes(2))%free = .FALSE.
      END IF
    END DO

!print *,'Vertices free (T/F):'
!print *,vertices(1:num_nodes)%free

    ! Now, count the free vertices:
    num_free_vertices = 0
    DO inode = 1,num_nodes 
      IF (vertices(inode)%free) THEN 
        num_free_vertices = num_free_vertices + 1
      END IF
    END DO

!print *,'Num free vertices:'
!print *,num_free_vertices

  END SUBROUTINE FIND_AND_COUNT_FREE_VERTICES

  SUBROUTINE COUNT_FREE_EDGES(num_free_edges)
    USE geometry
	IMPLICIT NONE
!*******************************************************************************
!*** Count free edges.                                                       ****
!*******************************************************************************
! Required for mixed-potential preconditioner.
!*******************************************************************************
    INTEGER(I4B), INTENT(OUT) :: num_free_edges
    INTEGER(I4B) iedge

    num_free_edges = 0

    ! Cycle edges:
    DO iedge = 1,num_edges 
      IF (edges(iedge)%free) THEN 
        num_free_edges = num_free_edges + 1
      END IF
    END DO

  END SUBROUTINE COUNT_FREE_EDGES
! End DBD added 4 Dec 02

! Start DBD added 6 Dec 02


   SUBROUTINE NODEEDGE_INDEXLIST_MAKE
     USE nrtype
     USE geometry
     IMPLICIT NONE   
!*******************************************************************************
! This subroutine builds a node-edge index list, indicating which edges
! are connected to each node, and which are the starting node and which 
! the end node (the edge is always directed from the lower to higher global 
! node number in FEMFEKO).
! Output is:
! * An array Node_edge_ptr(:) pointing to the start of each node's data in the node-edge list
! * An array Node_edge_nind(:) containing the edges associated with each node starting
!   at node one and moving upwards. This is the node-edge list.  
! * An array Node_edge_nind_sign(:) which flags the node as the start (-1)
!   or end (+1) of the edge. Otherwise as for Node_edge_nind
!
! Original version 6 Dec 2002 - DB Davidson. Routine adapted from Meyer's routine
!                               NODEELEMENT_INDEXLIST_MAKE
! For production code, this and next routines should be moved to GEOMETRY once tested.
!*******************************************************************************
   INTEGER(I4B) ii,kk ! more counters    
! Following line should eventually be moved to GEOMETRY
!   INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: Node_edge_ptr, Node_edge_ind, Node_edge_ind_sign  
   
   ! First check to see if the data structures are already created:
   IF (ALLOCATED(Node_edge_ind).AND.ALLOCATED(Node_edge_ptr)) RETURN

   ! Create a node-starting edge list
   ALLOCATE(Node_edge_ptr(0:num_nodes+1))     
   Node_edge_ptr(:) = 0
   DO ii=1,num_edges
     DO kk = 1,2
       Node_edge_ptr(edges(ii)%nodes(kk)) = Node_edge_ptr(edges(ii)%nodes(kk)) + 1
	 END DO
   END DO

   ! CSR format from long array
   DO ii=2,num_nodes
     Node_edge_ptr(ii) = Node_edge_ptr(ii) + Node_edge_ptr(ii-1)
   END DO
   DO ii=num_nodes+1, 2, -1   
     Node_edge_ptr(ii) = Node_edge_ptr(ii-1)
   END DO
   Node_edge_ptr(1) = 0 
     
   ! Create a node indexing array
   ALLOCATE(Node_edge_ind(Node_edge_ptr(num_nodes+1)))     
   ALLOCATE(Node_edge_ind_sign(Node_edge_ptr(num_nodes+1)))     
   DO ii=1, num_edges
     DO kk = 1,2
       ! Add one for we have to move one up in the index each time the same
       ! node is worked with
       Node_edge_ptr(edges(ii)%nodes(kk)) = Node_edge_ptr(edges(ii)%nodes(kk)) + 1     
       Node_edge_ind(Node_edge_ptr(edges(ii)%nodes(kk))) = ii
       SELECT CASE (kk)
	     CASE(1) ! Start node
           Node_edge_ind_sign(Node_edge_ptr(edges(ii)%nodes(kk))) = -1
	     CASE(2) ! End node
           Node_edge_ind_sign(Node_edge_ptr(edges(ii)%nodes(kk))) = +1
		 CASE DEFAULT
		   STOP 'Internal error in SUBROUTINE NODEEDGE_INDEXLIST_MAKE'
	   END SELECT 
	 END DO
   END DO
   
   ! Interesting. This restores the original nptr values.  
   DO ii=num_nodes+1, 2, -1   
     Node_edge_ptr(ii) = Node_edge_ptr(ii-1)
   END DO
   Node_edge_ptr(1) = 0


! write(fileout,*)'Node_edge_ptr',Node_edge_ptr
! write(fileout,*)'Node_edge_ind',Node_edge_ind
! write(fileout,*)'Node_edge_ind_sign',Node_edge_ind_sign



END SUBROUTINE NODEEDGE_INDEXLIST_MAKE
!******************************************************************************

SUBROUTINE NODEEDGE_INDEXLIST_CLEAN
  IMPLICIT NONE
!******************************************************************************
! Purpose: Clean up arrays of node-edge index list (allocated in subroutine
! MAKE_NODEEDGE_INDEXLIST).
! Author:  DBD 6 Dec 2002
!******************************************************************************

  IF (ALLOCATED(Node_edge_ind)) DEALLOCATE(Node_edge_ind)    
  IF (ALLOCATED(Node_edge_ind_sign)) DEALLOCATE(Node_edge_ind_sign)    
  IF (ALLOCATED(Node_edge_ptr)) DEALLOCATE(Node_edge_ptr)

END SUBROUTINE NODEEDGE_INDEXLIST_CLEAN
!******************************************************************************



! End DBD added 4-11 Dec 02


END SUBROUTINE ITER_SOLVE
!*******************************************************************************



SUBROUTINE ITER_SOLVE_DP(time_taken,smallest_its,largest_its,average_its)
  USE feminterface, ONLY: MATVECPROD
  USE geometry
  USE matrix
  USE math_tools, ONLY: TIME_DIFFERENCE
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! Iterative solvers - for sparse real double-precision matrices. 
! Solves [A][x] = [b].  
! Note that the matrix [A} is assumed to be SYMMETRIC. THIS IS NOT CHECKED!
! Presently implemented:
! - Conjugate gradient 
!
! References: CG: Algorithm version 3. J-M Jin. "The Finite Element Method
! in Electromagnetics", (p400ff)
! Bi-CG: algorithm version 1, [ibid.],
! but read paragraph 11.2.2 too!!)  - or equivalently fig 9.18 (p328) Volakis.
! 
! Real, double precision matrix A is stored in Asparse_DP; 
! RHS vector b is stroed in b_vec_DP;
! LHS vector x is x_vec_DP.
! NB! All must have been allocated before this routine is called, and A and b stored. 
! The parameters are passed via the module "matrix".
!
! The routine is implemented in double precision to permit the various precondioners
! supported by CXML to be used (these are all DP implementations). 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson
!
!*******************************************************************************
! Initial version: 1 May 2003
!
!*******************************************************************************
  REAL(SP), INTENT(OUT) :: time_taken     ! time taken by this routine. Note: SP retained here.
  INTEGER(I4B), INTENT(OUT), OPTIONAL :: smallest_its,largest_its,average_its
  INTEGER(I4B),SAVE :: cumulative_its
  INTEGER(I4B),SAVE :: num_calls
  
  LOGICAL(LGT) :: first_call = .TRUE.        ! Initialization confers SAVE attribute.
  ! Internal variable names (following Jin):
  REAL(DP), DIMENSION(:),ALLOCATABLE  :: rk,pk,zk,Apk,Aark
  REAL(DP) :: alphak,betak,gammak,rkzk1,rkzk2
  INTEGER(I4B) :: maxiterate,itcount,k

  REAL(DP) :: bnorm,error
  LOGICAL(LGT) success_termination ! Flag to check successful termination.
!  INTEGER(I4B), DIMENSION(8) :: time_start,time_finish ! temp variables for timing all events
  LOGICAL(LGT) :: preconditioning ! flag to indicate whether this option is enabled.
  REAL(DP),DIMENSION(:),ALLOCATABLE :: precond_sparse
  INTEGER(I4B), DIMENSION(8) :: time_start,time_finish ! temp variables for timing all events

  ! Record starting time:
  CALL DATE_AND_TIME (values=time_start)

  preconditioning = USE_PRE_CONDITIONER
  IF (preconditioning) STOP ! CALL SETUP_PRECONDITIONER

  success_termination  = .FALSE. 
  IF((SPARSE.AND..NOT.ALLOCATED(Asparse_DP)).OR.    &
     (.NOT.ALLOCATED(x_vec_DP)).OR.                 & 
     (.NOT.ALLOCATED(b_vec_DP))) THEN 
     CALL ERROR_FEMFEKO(1,4107)
  END IF 
  ALLOCATE(rk(dof))
  ALLOCATE(zk(dof))
  ALLOCATE(pk(dof))
  ALLOCATE(Apk(dof))
  ALLOCATE(Aark(dof)) ! Only for CG

  maxiterate = NINT(max_iter_factor*dof)
  ! Precompute the norm of b:
  bnorm = SQRT(ABS(DOT_PRODUCT(b_vec_DP,b_vec_DP)) )
  ! Initialise the iterated variables:
  x_vec_DP = 0.0_DP ! Arbitrary initialization, zero chosen as convenient


  SOLVER_CHOICE: IF (SOLVER_TYPE.EQ.1) THEN 
    STOP 'Unimplemented iterative solver type in USER_DITSOL_PCG'
   ELSE IF(SOLVER_TYPE.EQ.2) THEN ! CG
    ! This is [Jin, 2nd edn,p.610] CG Algorithm version 3.
    ! It assumes a real symmetric matrix. 
    ! The inner product is the conventional Euclidean norm.

    ! Initialize rk:  
    CALL MATVECPROD(5,x_vec_DP,rk,dof)
    rk = b_vec_DP - rk
    !Initialize pk:
    CALL MATVECPROD(5,rk,Aark,dof) ! For a real symmetric matrix, A_adjoint = A.
    pk = Aark/DOT_PRODUCT(Aark,Aark) 

    ! Iterate
    CG_iterate: DO itcount = 1,maxiterate
      ! Calculation of [A]{pk}:
      CALL MATVECPROD(5,pk,Apk,dof)
      alphak = 1/DOT_PRODUCT(Apk,Apk)
      x_vec_DP = x_vec_DP + alphak*pk ! xk now overwritten with xk+1
      rk = rk - alphak*Apk            ! rk now overwritten with rk+1

      ! Calculation of [A^a][rk+1]:
      CALL MATVECPROD(5,rk,Aark,dof)
      betak = 1/DOT_PRODUCT(Aark,Aark)
      pk = pk + betak*Aark          ! pk now overwritten with pk+1
      error = SQRT(ABS(DOT_product(rk,rk)))/bnorm
      IF (ON_SCREEN_REPORTING) THEN
        print *, itcount,'/',maxiterate, error
      END IF
      IF (error.LT.residual_norm) THEN
        success_termination = .TRUE.
        EXIT CG_iterate
      END IF
    END DO CG_iterate

  ELSE IF(SOLVER_TYPE.EQ.3) THEN ! QMR
    STOP 'Unimplemented iterative solver type in USER_DITSOL_PCG'
    ! CALL QMR_SOLVE

  ELSE IF(SOLVER_TYPE.EQ.4) THEN ! GMRES
    STOP 'Unimplemented iterative solver type in USER_DITSOL_PCG'
    ! CALL GMRES_SOLVE

   ELSE IF(SOLVER_TYPE.EQ.5) THEN ! CG
    STOP 'Unimplemented iterative solver type in USER_DITSOL_PCG'
  END IF SOLVER_CHOICE 

  ! Record finishing time:
  CALL DATE_AND_TIME (values=time_finish)
  time_taken = TIME_DIFFERENCE(time_start,time_finish)

  IF(FIRST_CALL) THEN 
  ! Write solution information to the output file:
    WRITE (FILEOUT,'(//,20X,A)')   'ITERATIVE SOLUTION OF THE SYSTEM MATRIX EQUATION'
    SELECT CASE (SOLVER_TYPE)
    CASE (1)
      WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Bi-conjugate gradient'
    CASE (2)
      WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Conjugate gradient'
    CASE (3)
      WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Quasi-Minimal Residual'
    CASE (4)
      WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Generalized Minimum Residual'
      WRITE (FILEOUT,'(1X,A,I4)')     'Iterations before restart:',restart_GMRES
    CASE (5)
      WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  Mixed Potential Conjugate gradient'
    CASE DEFAULT
      STOP 'IE: Invalid iterative solver requested.'
    END SELECT
    IF(USE_PRE_CONDITIONER) THEN
      WRITE (FILEOUT,'(1X,A)')     'Diagonal preconditioning used.'
    ELSE
      WRITE (FILEOUT,'(1X,A)')     'No preconditioning used.'
    END IF     
    WRITE (FILEOUT,'(1X,A,E9.3)')    'Specified residual norm:         ', residual_norm
    WRITE (FILEOUT,'(1X,A,I12)')     'Specified max. iterations:       ', maxiterate
    ! WRITE (FILEOUT,'(1X,A,E9.3)')    'Final residual norm:             ', error
    ! WRITE (FILEOUT,'(1X,A,I12)')     'Final num. iterations:           ', itcount
    smallest_its = itcount
    largest_its  = itcount
    cumulative_its  = 0
	num_calls = 0
	FIRST_CALL = .FALSE. ! So that the above is only printed and/or intialized once.
  END IF

  ! Record some statistics about the iterative solver.
  IF(itcount.LT.smallest_its) THEN 
    smallest_its = itcount
  END IF
  IF(itcount.GT.largest_its) THEN 
    largest_its = itcount
  END IF
  cumulative_its = cumulative_its + itcount
  num_calls      = num_calls + 1
  average_its    = cumulative_its / num_calls ! Rounded down
  IF (itcount.GE.maxiterate) THEN
    CALL ERROR_FEMFEKO(0,4108,num_calls) ! Print warning of unconverged iteration,
  END IF  
  ! Deallocate local, allocatable variables:
  DEALLOCATE(rk)
  DEALLOCATE(zk)
  DEALLOCATE(pk)
  DEALLOCATE(Apk)
  DEALLOCATE(Aark)
  IF (preconditioning) DEALLOCATE(precond_sparse)
END SUBROUTINE ITER_SOLVE_DP

