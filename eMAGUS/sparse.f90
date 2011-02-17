! Last changed 09 Feb 2001 MMB - FJCM updates


SUBROUTINE sparse_ARPACK(howmny, SELECT, d, z, n, numev, resid, numcv, &
                          v, ldv, iparam, workd, workl, lworkl, iwork, info)

USE output_error
!*************************************************************************
!  Subroutine description
!*************************************************************************
! 
! This routine supplies the routine ssaupd.f with the matrix-vector products
! as requested by the reverse communication flag ido. 
! Currently it only supports the shift-invert mode of ARPACK.     
!
!*************************************************************************
! AUTHOR
!*************************************************************************
!
! Modified ssband.f supplied with ARPACK
! Riana Helena Geschke
!
!*************************************************************************
! Last revised:
!*************************************************************************
!
!  8 Mar 2000
! 20 Jul 2000: MATVECPROD arguments 4 & 5 switched. DBD.
! 3 August 2005: Added iterative solver capability to shift invert linear system. JPS.
! 10 August 2005: Added sparse direct linear solver using umfpack with c wrapper functions.
!
!*************************************************************************


USE feminterface, ONLY : MATVECPROD, LOWERSOLVE2, UPPERSOLVE2, ITER_SOLVE
USE nrtype
USE matrix
USE unit_numbers
USE problem_info
IMPLICIT NONE


   Character,                 INTENT(IN)      :: howmny      
   ! FJCM
   ! Had to swop the following two command lines for DF to accept
   ! Was
   ! Integer(I4B),              INTENT(IN)      :: n, numev, &
   !                              numcv, ldv, lworkl
   ! Logical(LGT),DIMENSION(numcv),INTENT(INOUT):: select
   ! is 
   Integer(I4B),              INTENT(IN)      :: n, numev, &
                                 numcv, ldv, lworkl
   Logical(LGT),DIMENSION(numcv),INTENT(INOUT):: select
   ! 
   INTEGER(I4B) , INTENT(INOUT)               :: info  
   Integer(I4B), DIMENSION(11), INTENT(INOUT) :: iparam
   Integer(I4B), DIMENSION(n), INTENT(INOUT)  :: iwork
   Real(SP), DIMENSION(lworkl),INTENT(INOUT)  :: workl                 
   Real(SP), DIMENSION(numev),  INTENT(INOUT) :: d
   Real(SP), DIMENSION(n) ,    INTENT(INOUT)  :: resid
   Real(SP), DIMENSION(ldv,numcv), INTENT(INOUT) :: v  
   Real(SP), DIMENSION(3*n), INTENT(INOUT)    :: workd   
   Real(SP), DIMENSION(n,numev), INTENT(INOUT) :: z
   REAL(SP) :: time_taken  ! added by jps to test subroutine ITER_SOLVE
   
   !local variables
   Real(SP), DIMENSION(n)                     :: tempvecA, tempvecB
   Integer(I4B), DIMENSION(14)                :: ipntr
   Integer(I4B)                               :: ido, i, type, ierr, nconv , printk, count  
   REAL(SP)                                   :: tol
   !jps added the following 2 lines for timing
   REAL(SP) :: timer1, timer2, time_ido1, time_ido2
   INTEGER(I4B):: counter, counter_rate, ndim, nnzdim, iter_count
 

   ! added for luinc and pcg_solve
   REAL(SP), ALLOCATABLE :: jluvalues(:)
   INTEGER(I4B), DIMENSION(n) :: diagpntr
   INTEGER(I4B) :: jinfo

  
   ! added for umfpack
   INTEGER(PTR) :: umf_numeric, umf_symbolic, umf_sys
   REAL(DP), DIMENSION(20) :: umf_control
   REAL(DP), DIMENSION(90) :: umf_info
   INTEGER(PTR) :: umf_n
   INTEGER(PTR) :: umf_nnz
   INTEGER(PTR), DIMENSION(n+1) :: umf_col_ptr
   INTEGER(PTR), ALLOCATABLE :: umf_row_index(:)
   REAL(DP), ALLOCATABLE :: umf_values(:)
   REAL(DP), DIMENSION(n) :: umf_b, umf_x
 
   print *,'ENTERING SPARSE_D'

   time_ido1 = 0.0
   time_ido2 = 0.0
   iter_count = 0.0

   count = 0
   print *,'n= ', n
   !initialisation
   tempvecA = 0_SP 
   tempvecB = 0_SP 
   
   resid = 0_SP
   workd = 0_SP
   ipntr = 0  
 

   ! determine which mode will be used.
 
   if ( iparam(7) .eq. 1 ) then
       type = 1
    else if ( iparam(7) .eq. 3 .and. bmat .eq. 'I') then
       type = 2
    else if ( iparam(7) .eq. 2 ) then
       type = 3
    else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
       type = 4     
    else if ( iparam(7) .eq. 4 ) then
       type = 5
    else if ( iparam(7) .eq. 5 ) then 
       type = 6
    else  ! still update this to use error reporting
       print*, ' '
       print*, 'bmat is inconsistent with IPARAM(7).'
       print*, ' ' 
       go to 9000
    end if
 
    print *,'type = ', type 
    print *,'sigma = ', sigma

    ! added by JPS
    ! if we use shift invert mode, allocate and populate sparse matrix S_min_sigT
    ! to be used in the linear system solves at each time step
    ndim = dof
    nnzdim = SIZE(col_ind)

    IF (TYPE .EQ. 4) THEN
       ALLOCATE(S_min_sigT(nnzdim))
       S_min_sigT(1:nnzdim) = Ssparse(1:nnzdim) - sigma*Tsparse(1:nnzdim)
       ! allocate luinc values array
       !ALLOCATE(jluvalues(nnzdim))
       ! perform incomplete LU decomposition
       !CALL luinc(ndim, nnzdim, S_min_sigT, col_ind, row_ind, jluvalues, diagpntr, jinfo)
       !WRITE(*,*) "ILU(0) routine completed with errorcode = ", jinfo
    
       tol = 1.0e-8   ! tolerance explicitly passed to ARPACK
       WRITE(*,*) "tolerance passed to arpack = ", tol


       !############################
       ! umfpack variable initialization and array allocation
       umf_n = ndim
       umf_nnz = nnzdim
       ALLOCATE(umf_values(umf_nnz))
       umf_values = S_min_sigT
       ALLOCATE(umf_row_index(umf_nnz))
       umf_row_index = col_ind-1 !swap row and col for CCS and adjust to 0 based indices
       umf_col_ptr = row_ind-1   !this is a quick and dirty hack that exploits symmetry of matrix
       umf_sys = 0


       !----------------------------------------------------------------
       ! Perform symbolic factorization of the matrix using umfpack
       !----------------------------------------------------------------

       ! set default parameters
       CALL umf4def (umf_control)


       ! pre-order and symbolic analysis
       CALL umf4sym (umf_n, umf_n, umf_col_ptr, umf_row_index, umf_values, umf_symbolic, umf_control, umf_info)

       ! check umf4sym error condition
       IF (umf_info (1) < 0) THEN
          WRITE(*,*) 'Error occurred in umf4sym: ', umf_info (1)
          STOP
       ELSE
          WRITE(*,*) 'Symbolic factorization completed successfully.'
       ENDIF

       !-------------------------------------------------------------
       ! perform numeric LU factorization using umfpack
       !------------------------------------------------------------

       CALL umf4num (umf_col_ptr, umf_row_index, umf_values, umf_symbolic, umf_numeric, umf_control, umf_info)

       ! check umf4num error condition
       IF (umf_info (1) < 0) THEN
          WRITE(*,*) 'Error occurred in umf4num: ', umf_info (1)
          STOP
       ELSE
          WRITE(*,*) 'Numeric factorization completed successfully.'
       ENDIF


       ! free the memory allocated for symbolic analysis
       CALL umf4fsym (umf_symbolic)

       ! end intitial umfpack initialisation and computation
       !####################################################
    END IF
    
    
   
    !initialise the reverse communication flag
     ido   = 0 
    
    !--------------------------------------------%
    !  M A I N   L O O P (reverse communication) |
    !--------------------------------------------%

 90  continue 
      call ssaupd (ido, bmat, dof, which, numev, tol, resid, numcv, &
                   v, ldv, iparam, ipntr, workd, workl, lworkl, &
                   info)

      
      if (ido .eq. -1) then
         if ( type .eq. 3 ) then
         
           !-----------------------------------------%
           ! Perform  y <--- OP*x = inv[T]*S*x       |
           ! to force the starting vector into       | 
           ! the range of OP.                        |
           !-----------------------------------------%
            call MATVECPROD(1,workd(ipntr(1):ipntr(1)+n-1),tempvecA,dof)
            call scopy(n, tempvecA(1), 1, workd(ipntr(1)), 1)
            call lowersolve2(n,workd(ipntr(1):ipntr(1)+n-1),tempvecB)
            call uppersolve2(n,tempvecB,workd(ipntr(2):ipntr(2)+n-1))
          
          else if ( type .EQ. 4) then
           !---------------------------------------------%
           ! Perform  y <--- OP*x = inv[S - sigma*T]*T*x |
           !---------------------------------------------%
           !call MATVECPROD(2,workd(ipntr(1):ipntr(1)+n-1),tempvecA,dof,sigma)
           CALL Mv(ndim, nnzdim, Tsparse(1:nnzdim), col_ind, row_ind, &
                   workd(ipntr(1):ipntr(1)+n-1),tempvecA)
           !CALL scopy(n, tempvecA(1), 1, workd(ipntr(1)), 1)
           !call lowersolve2(n,workd(ipntr(1):ipntr(1)+n-1),tempvecB)
           !call uppersolve2(n,tempvecB,workd(ipntr(2):ipntr(2)+n-1)) 
           !WRITE(*,*) "calling pcg_solve"
           !CALL pcg_solve(ndim, nnzdim, jluvalues, col_ind, row_ind, diagpntr, S_min_sigT(1:nnzdim), col_ind, row_ind, &
           !              workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1))


           !###########################
           ! start umfpack solution
           !############################
           
           ! convert to rhs vector to double precision for use in umfpack solver
           umf_b = REAL(workd(ipntr(1):ipntr(1)+n-1), DP)

           ! solve Ax=b, without iterative refinement
           CALL umf4sol (umf_sys, umf_x, umf_b, umf_numeric, umf_control, umf_info)

           ! check umf4sol error condition
           IF (umf_info (1) < 0) THEN
              WRITE(*,*) 'Error occurred in umf4sol: ', umf_info (1)
              STOP
           ELSE
              !WRITE(*,*) 'Linear system solution completed successfully.'
           ENDIF

           !convert type and store in proper workd space.
           workd(ipntr(2):ipntr(2)+n-1) = REAL(umf_x, SP)

           !######################
           ! end umfpack solution
           !######################


        end if

      else if (ido .eq. 1) then
         if ( type .eq. 3 ) then
           !-----------------------------------------%
           ! Perform  y <--- OP*x = inv[T]*S*x       |
           !-----------------------------------------%
           call MATVECPROD(1,workd(ipntr(1):ipntr(1)+n-1),tempvecA,dof,sigma)
           call scopy(n, tempvecA(1), 1, workd(ipntr(1)), 1)
           call lowersolve2(n,workd(ipntr(1):ipntr(1)+n-1),tempvecB)
           call uppersolve2(n,tempvecB,workd(ipntr(2):ipntr(2)+n-1))
           
         else if ( type .EQ. 4) then
           !---------------------------------------------%
           ! Perform  y <--- OP*x = inv[S - sigma*T]*T*x |
           !---------------------------------------------%
           CALL SYSTEM_CLOCK(counter,counter_rate)                    ! added by jps to time the linear system solution
           timer1 = REAL(counter, SP)/REAL(counter_rate, SP)  
           !call lowersolve2(n,workd(ipntr(3):ipntr(3)+n-1),tempvecB)
           !call uppersolve2(n,tempvecB,workd(ipntr(2):ipntr(2)+n-1)) 
           !CALL pcg_solve(ndim, nnzdim, jluvalues, col_ind, row_ind, diagpntr, S_min_sigT(1:nnzdim), col_ind, row_ind, &
           !              workd(ipntr(3):ipntr(3)+n-1), workd(ipntr(2):ipntr(2)+n-1))
           

           !###########################
           ! start umfpack solution
           !############################
           
           ! convert to double precision for use in umfpack solver
           umf_b = REAL(workd(ipntr(3):ipntr(3)+n-1), DP)

           ! solve Ax=b, without iterative refinement
           CALL umf4sol (umf_sys, umf_x, umf_b, umf_numeric, umf_control, umf_info)

           ! check umf4sol error condition
           IF (umf_info (1) < 0) THEN
              WRITE(*,*) 'Error occurred in umf4sol: ', umf_info (1)
              STOP
           ELSE
              !WRITE(*,*) 'Linear system solution completed successfully.'
           ENDIF

           !convert type and store in proper workd space
           workd(ipntr(2):ipntr(2)+n-1) = REAL(umf_x, SP)

           !######################
           ! end umfpack solution
           !######################


           CALL SYSTEM_CLOCK(counter,counter_rate)
           timer2 = REAL(counter, SP)/REAL(counter_rate, SP)
           time_ido1 = time_ido1 + (timer2 - timer1)

         end if

      else if (ido .eq. 2) then
        !----------------------------------%
        !        Perform y <-- T*x         | 
        !                                  | 
        !----------------------------------%
         CALL SYSTEM_CLOCK(counter,counter_rate)                    ! added by jps to time the linear system solution
         timer1 = REAL(counter, SP)/REAL(counter_rate, SP)  
        
         !CALL MATVECPROD(2,workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),dof,sigma)
         CALL Mv(ndim, nnzdim, Tsparse(1:nnzdim-1), col_ind, row_ind, &
                 workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1))
         iter_count = iter_count + 1
         !WRITE(*,*) "Called subroutine Mv and iter_count = ", iter_count
         CALL SYSTEM_CLOCK(counter,counter_rate)
         timer2 = REAL(counter, SP)/REAL(counter_rate, SP)
         time_ido2 = time_ido2 + (timer2 - timer1)
      else 
        !-----------------------------------------%
        ! Either we have convergence, or there is | 
        ! error.                                  |
        !-----------------------------------------%
         
         if ( info .lt. 0) then
           !--------------------------%
           ! Error message, check the |
           ! documentation in SSAUPD  |
           !--------------------------%

            print *, ' '
            print *, ' Error with _saupd info = ',info
            print *, ' Check the documentation of _saupd '
            print *, ' '
            go to 9000
         else 
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit', &
                       ' Arnoldi update, try increasing numcv.'
               print *, ' '
            end if
           
             
            if (iparam(5) .gt. 0) then
              
     
               call sseupd ( rvec, 'A', select, d, eigenvectors, n, sigma, &
                       bmat, n, which, numev, tol, resid, numcv, v, ldv,& 
                       iparam, ipntr, workd, workl, lworkl, info )            
              
 
               if ( info .ne. 0) then
 
                 !------------------------------------%
                 ! Check the documentation of sneupd. |
                 !------------------------------------%

                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
 
               end if

            end if

         end if

         go to 9000

      end if

!     %----------------------------------------%
!     | L O O P  B A C K to call SSAUPD again. |
!     %----------------------------------------%

      go to 90 

 9000 continue

 IF ( info .eq. 0) THEN
    
    !######################
    ! final umfpack cleanup
    ! free the memory allocated for numeric factorization
    CALL umf4fnum (umf_numeric)
    !#######################

  !****************************************************************************
  !Print out convergence information 
  !****************************************************************************
   
   nconv = iparam(5)
   print *, ' Convergence information***************************'
   print *, ' Number of eigenvalue requested is ', numev
   print *, ' The number of Lanczos vectors generated',' (numcv) is ', numcv
   print *, ' The number of converged Ritz values is ',nconv
   print *, ' What portion of the spectrum ', which
   print *, ' The number of Implicit Arnoldi', ' update taken is ', iparam(3)
   print *, ' The number of OP*x is ', iparam(9)
   print *, ' The convergence tolerance is ', tol
   PRINT *, ' Total time taken for solving the LU factored systems is ', time_ido1
   print *, ' '
   PRINT *, 'ndim = ', ndim
   PRINT *, 'nnzdim = ', nnzdim
   PRINT *, 'SIZE(Tsparse) = ', SIZE(Tsparse)
   PRINT *, 'SIZE(Ssparse) = ', SIZE(Ssparse)
   PRINT *, 'SIZE(col_ind) = ', SIZE(col_ind)
   PRINT *, 'SIZE(row_ind) = ', SIZE(row_ind)
   PRINT *, 'col_ind(1) = ', col_ind(1)
   PRINT *, 'col_ind(nnzdim) = ', col_ind(nnzdim)
   PRINT *, 'Tsparse(nnzdim+1) = ', Tsparse(nnzdim+1)
   PRINT *, 'Ssparse(nnzdim+1) = ', Ssparse(nnzdim+1)
   PRINT *, 'MINVAL(Tsparse) = ', MINVAL(Tsparse(1:nnzdim))
   PRINT *, 'MINVAL(Ssparse) = ', MINVAL(Ssparse(1:nnzdim))
   
   do printk = 1,nconv
     print *,'eigenvalue(',printk,')= ', sqrt(d(printk))
   end do
   nconv_eigenvalues = nconv 
   eigenvalues = 0
   eigenvalues(nev-nconv+1:nev) = d(1:nconv)
   !eigenvalues(1:nconv) = d(1:nconv)


ELSE
   print *,'info = ', info
END IF


END SUBROUTINE SPARSE_ARPACK

!*********************************************************************
! jps added subroutines Mv, cg_solve, pcg_solve, luinc, lusol, tofull 
! to source file sparse.f90. These subroutines are currently only used 
! by subroutine sparse_ARPACK
!*********************************************************************

SUBROUTINE cg_solve(n, nnz, Avalues, colindex, rowindex, b, x)
! subroutine to solve a linear system A*x=b using
! the conjugate gradient algorithm
! note that A is a sparse matrix stored in CRS format
! written by Julian P. Swartz
! last modified on 26 July 2005

  USE nrtype

  IMPLICIT NONE
  
  ! input and output parameter declarations
  INTEGER(I4B), INTENT(IN) :: n, nnz
  REAL(SP), INTENT(IN), DIMENSION(nnz) :: Avalues
  INTEGER(I4B), INTENT(IN), DIMENSION(nnz) :: colindex
  INTEGER(I4B), INTENT(IN), DIMENSION(n+1) :: rowindex
  REAL(SP), INTENT(IN), DIMENSION(n) :: b
  REAL(SP), INTENT(OUT), DIMENSION(n) :: x

  ! internal variable declarations
  INTEGER(I4B) :: index, maxit    ! loop index and max number of iterations
  REAL(SP) :: alpha, beta, tol
  REAL(SP), DIMENSION(n) :: xold, xnew, pold, pnew, rold, rnew, tempv
  
  ! set tol and maxit internally for now
  tol = 1.0E-6
  maxit = 100000
  
  ! initialise variables
  xold(1:n) = 0.0
  CALL Mv(n, nnz, Avalues, colindex, rowindex, xold, tempv)
  rold = b - tempv
  pold = rold
  index = 1

  ! Loop until norm of residual is within tolerance or
  ! maxit is reached
  CG_LOOP: DO
     CALL Mv(n, nnz, Avalues, colindex, rowindex, pold, tempv)
     alpha = DOT_PRODUCT(rold, rold) / DOT_PRODUCT(tempv, pold)
     xnew = xold + alpha*pold
     rnew = rold - alpha * tempv

     ! if the norm (squared) of the residual is less than tolerance, stop
     IF (DOT_PRODUCT(rnew, rnew) < tol) THEN
        WRITE(*,*) "cg_solve: iterations to convergence ", index 
        EXIT CG_LOOP        
     END IF

     ! if maxit reached, terminate and print error message
     IF (index == maxit) THEN
        WRITE(*, *) "cg_solve did not converge within maxit iterations."
        STOP
     END IF
     
     !update xold
     xold = xnew
            
     beta = DOT_PRODUCT(rnew, rnew) / DOT_PRODUCT(rold, rold)

     !update rold
     rold = rnew
            
     pnew = rnew + beta*pold

     !update pold
     pold = pnew
     index = index+1

  END DO CG_LOOP

  x = xnew


END SUBROUTINE cg_solve


SUBROUTINE pcg_solve(n, nnz, Mvalues, mcolind, mrowind, mdiagpntr, Avalues, colindex, rowindex, b, x)
! subroutine to solve a linear system A*x=b using
! the conjugate gradient algorithm
! note that A is a sparse matrix stored in CRS format
! A ILU preconditioner is passed in as a sparse matrix
! with associated values, column and row rowpointer arrays
! along with a pointer array to the diagonal elements
! the L and U factors are stored in the same sparse matrix structure
! written by Julian P. Swartz
! last modified on 1 August 2005

  USE nrtype

  IMPLICIT NONE
  
  ! input and output parameter declarations
  INTEGER(I4B), INTENT(IN) :: n, nnz
  REAL(SP), INTENT(IN), DIMENSION(nnz) :: Avalues, Mvalues
  INTEGER(I4B), INTENT(IN), DIMENSION(nnz) :: colindex, mcolind
  INTEGER(I4B), INTENT(IN), DIMENSION(n+1) :: rowindex, mrowind
  REAL(SP), INTENT(IN), DIMENSION(n) :: b, mdiagpntr
  REAL(SP), INTENT(OUT), DIMENSION(n) :: x

  ! internal variable declarations
  INTEGER(I4B) :: index, maxit    ! loop index and max number of iterations
  REAL(SP) :: alpha, beta, tol
  REAL(SP), DIMENSION(n) :: xold, xnew, pold, pnew, rold, rnew, tempv, zold, znew
  
  ! set tol and maxit internally for now
  tol = 1.0E-6
  maxit = 100000
  
  ! initialise variables
  xold = 0.0
  CALL Mv(n, nnz, Avalues, colindex, rowindex, xold, tempv)
  rold = b - tempv
 
  CALL lusol(n, nnz, rold, zold, Mvalues, mcolind, mrowind, mdiagpntr)
  pold = zold
  index = 1
 
  ! Loop until norm of residual is within tolerance or
  ! maxit is reached
  CG_LOOP: DO
     CALL Mv(n, nnz, Avalues, colindex, rowindex, pold, tempv)
     alpha = DOT_PRODUCT(rold, zold) / DOT_PRODUCT(tempv, pold)
     xnew = xold + alpha*pold
     rnew = rold - alpha * tempv

     ! update z
     CALL lusol(n, nnz, rnew, znew, Mvalues, mcolind, mrowind, mdiagpntr)

     ! if the norm (squared) of the residual is less than tolerance, stop
     IF (DOT_PRODUCT(rnew, rnew) < tol) THEN
        WRITE(*,*) "pcg_solve: iterations to convergence ", index 
        EXIT CG_LOOP        
     END IF

     ! if maxit reached, terminate and print error message
     IF (index == maxit) THEN
        WRITE(*, *) "pcg_solve did not converge within maxit iterations."
        WRITE(*,*) "dot(rnew, rnew) = ", DOT_PRODUCT(rnew, rnew)
        STOP
     END IF
     
     !update xold
     xold = xnew
            
     beta = DOT_PRODUCT(rnew, znew) / DOT_PRODUCT(rold, zold)

     !update rold and zold
     rold = rnew
     zold = znew
            
     pnew = znew + beta*pold
     !update pold
     pold = pnew
     index = index+1

  END DO CG_LOOP
 
  x = xnew
  

END SUBROUTINE pcg_solve


SUBROUTINE Mv(n, nnz, Mvalues, colindex, rowindex, v, w)
! Subroutine to compute sparse matrix-vector product
! w = M*v for sparse square matrices of dimension n
! written by Julian P. Swartz
! last modified on 26 July 2005
  
  USE nrtype

  IMPLICIT NONE

  ! input and output parameter declarations
  INTEGER(I4B), INTENT(IN) :: n, nnz
  REAL(SP), INTENT(IN), DIMENSION(nnz) :: Mvalues
  INTEGER(I4B), INTENT(IN), DIMENSION(nnz) :: colindex
  INTEGER(I4B), INTENT(IN), DIMENSION(n+1) :: rowindex
  REAL(SP), INTENT(IN), DIMENSION(n) :: v
  REAL(SP), INTENT(OUT), DIMENSION(n) :: w

  ! internal variable declarations
  INTEGER(I4B) :: i, startindex, endindex

  DO i = 1, n
     startindex = rowindex(i)
     endindex = rowindex(i+1) - 1
     w(i) = DOT_PRODUCT(Mvalues(startindex:endindex), v(colindex(startindex:endindex)))
  END DO

END SUBROUTINE Mv

SUBROUTINE tofull(n, nnz, values, colindex, rowindex, M)
! Subroutine to convert a sparse matrix stored in CRS
! format to a full representation
! written by Julian P. Swartz
! last modified on 26 July 2005

  USE nrtype

  IMPLICIT NONE
  ! input and output parameter declarations
  INTEGER(I4B), INTENT(IN) :: n, nnz
  REAL(SP), INTENT(IN), DIMENSION(nnz) :: values
  INTEGER(I4B), INTENT(IN), DIMENSION(nnz) :: colindex
  INTEGER(I4B), INTENT(IN), DIMENSION(n+1) :: rowindex
  REAL(SP), INTENT(OUT), DIMENSION(n, n) :: M
  
  ! internal variable declarations
  INTEGER(I4B) :: i, row, col
  
  M = 0.0
  row = 1

  DO i=1, nnz
     IF (i == rowindex(row+1)) THEN
        row = row+1
     END IF
     col = colindex(i)
     M(row, col) = values(i)    
  END DO
  
END SUBROUTINE tofull

SUBROUTINE luinc(n, nnz, Avalues, colindex, rowindex, luval, uptr, errcode)
! function to perform the level zero fill in, ILU(0) factorization
! of the sparse matrix stored in CSR format
! follows algorithm of chapter 10 of "Iterative Methods
! for Sparse Linear Systems" by Y. Saad

  USE nrtype

  IMPLICIT NONE

  ! input and output parameter declarations
  INTEGER(I4B), INTENT(IN) :: n, nnz
  REAL(SP), INTENT(IN), DIMENSION(nnz) :: Avalues
  INTEGER(I4B), INTENT(IN), DIMENSION(nnz) :: colindex
  INTEGER(I4B), INTENT(IN), DIMENSION(n+1) :: rowindex
  REAL(SP), INTENT(OUT), DIMENSION(nnz) :: luval
  INTEGER(I4B), INTENT(OUT) :: errcode
  INTEGER(I4B), INTENT(OUT), DIMENSION(n) :: uptr


  ! internal variable declarations
  INTEGER(I4B) :: i, startindex, endindex, ji, k, flag, jrow, jj, jw
  INTEGER(I4B), DIMENSION(n) :: iw
  REAL(SP) :: tl



  ! initialize luval to nonzero values array of matrix        
  luval(1:nnz) = Avalues(1:nnz)
  ! initialize work array to zero
  iw = 0.0
  ! initialize array to hold diagonal entry indices
  uptr = 0
  ! initialise errcode
  errcode = 0

  !# main loop
  OUTER: DO k = 1,n
     !# find row start and end index in self.values
     startindex = rowindex(k)
     endindex = rowindex(k+1)-1
     DO ji = startindex, endindex
        iw(colindex(ji)) = ji
     END DO
     ji = startindex        
     flag = 0
            
     INNER: DO
        IF (flag /= 0) EXIT INNER
        
        jrow = colindex(ji)
                
        ! exit if diagonal element is reached
        IF (jrow >= k) EXIT INNER
           ! goto 200 section
                    
        !# compute the multiplier for jrow
        tl = luval(ji) * luval(uptr(jrow))   !# note, t-ell not t-one
        luval(ji) = tl

        !# perform linear combination
        DO jj = uptr(jrow)+1, rowindex(jrow+1)-1
           jw = INT(iw(colindex(jj)))
           IF (jw /= 0) THEN
              luval(jw) = luval(jw) - tl*luval(jj)
           END IF
        END DO
        ji = ji+1
        IF (ji <= endindex) THEN
           CYCLE INNER
        ELSE
           !# end of row reached
           flag = 1
           EXIT INNER
        END IF
     END DO INNER
                        
     !# store pointer to diagonal element
     uptr(k) = ji
     IF ((jrow /= k) .OR. (luval(ji) == 0.0)) THEN
        errcode = 1
        WRITE(*,*) "Error in computing ILU(0)"
        EXIT OUTER
     END IF
     luval(ji) = 1.0/luval(ji)
     !# refresh entries of iw to zero
     DO i = startindex, endindex
        iw(colindex(i)) = 0.0
     END DO

     !# exit on error
     IF (errcode == 1) THEN
        WRITE(*,*) "Error in computing ILU(0)"
        EXIT OUTER
     END IF

  END DO OUTER

  !# end function and return result
  IF (errcode == 0) THEN
     !# implement ugly hack to fix bug in ILU(0) code
     !# that causes the diagonal elements to appear inverted
     DO i = 1, SIZE(uptr)
        luval(uptr(i)) = 1.0/luval(uptr(i))
     END DO
     WRITE(*,*) "ILU(0) successful"
  END IF

END SUBROUTINE luinc

SUBROUTINE lusol(n, nnz, rhs, sol, luval, colindex, rowindex, diagptr)
! Subroutine to solve a linear system A*x=b where A = L*U using
! forward and backward substitution.
! L is unit lower triangular and U is upper triangular. Both
! matrices are sparse and are stored in the same sparse matrix
! data structure (omitting the ones on the diagonal of L).
! diagptr is a pointer to the position of the  diagonal elements
! in luval.
! Written by Julian P. Swartz
! last modified on 1 August 2005

  USE nrtype

  IMPLICIT NONE

  ! input and output parameter declarations
  INTEGER, INTENT(IN) :: n, nnz
  REAL, INTENT(IN), DIMENSION(nnz) :: luval
  INTEGER, INTENT(IN), DIMENSION(nnz) :: colindex
  INTEGER, INTENT(IN), DIMENSION(n+1) :: rowindex
  INTEGER, INTENT(IN), DIMENSION(n) :: diagptr
  REAL, INTENT(IN), DIMENSION(n) :: rhs 
  REAL, INTENT(OUT), DIMENSION(n) :: sol

  ! internal variable declarations
  INTEGER :: i, k

  ! perform forward substitution solution
  ! sol(i) = rhs(i) - sum (L(i,j) * sol(j))
  DO i = 1, n
     sol(i) = rhs(i)
     DO k = rowindex(i), diagptr(i)-1
        sol(i) = sol(i) - luval(k)*sol(colindex(k))
     END DO
  END DO
  
  ! perform backward substitution solution
  ! sol(i) = sol(i) - sum(U(i,j)*sol(j))
  DO i = n, 1, -1
     DO k = diagptr(i)+1, rowindex(i+1)-1
        sol(i) = sol(i) - luval(k)*sol(colindex(k))
     END DO
     ! sol(i) = sol(i) / U(i,i)
     sol(i) = sol(i) / luval(diagptr(i))
  END DO

END SUBROUTINE lusol

