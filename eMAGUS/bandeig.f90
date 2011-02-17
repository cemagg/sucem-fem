! Last revised 14 March 00 DBD.


SUBROUTINE BAND_EIGENVALUES
   USE nrtype
   USE output_error
   USE geometry
   USE bandwidth
   USE unit_numbers
   USE matrix
   USE problem_info
   IMPLICIT NONE              
 
!*************************************************************************
!  Subroutine description
!*************************************************************************
!
! This routine solves the generalised eigenvalue problem in the ARPACK
! regular inverse mode or shift-invert mode
! 
! The band storage LAPACK scheme was used to store the S and T matrices
!
!
!*************************************************************************
!  AUTHOR
!*************************************************************************
!
! Riana Helena Geschke  (this is based on a band symmetric 
! Fortran 77 example included in the ARPACK/EXAMPLES/BAND directory)
!
!*************************************************************************
!  Last revised:
!*************************************************************************
! 14 March 00: j -> jj DBD.
!
!
!*************************************************************************
!  Input
!*************************************************************************
!
!  S and T banded matrices
!
!*************************************************************************
!  Output
!*************************************************************************
!
!  The converged eigenvalues, residuals and eigenvectors if desired
!
!*************************************************************************



!     ... Construct the matrix A in LAPACK-style band form.
!     ... Call SSBAND with regular mode to find eigenvalues LAMBDA 
!         such that
!                          A*x = LAMBDA*M*x.
!
!
!
!\Routines called:
!     ssband  ARPACK banded eigenproblem solver.
!     slaset  LAPACK routine to initialize a matrix to zero.
!     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sgbmv   Level 2 BLAS that computes the band matrix vector product

      !********************************************************************
      !  Define leading dimensions for all   
      !  arrays.                             
      !  MAXN   - Maximum size of the matrix 
      !  MAXNEV - Maximum number of          
      !               eigenvalues to be computed 
      !  MAXNCV - Maximum number of Arnoldi  
      !               vectors stored              
      !  MAXBDW - Maximum bandwidth          
      ! 
      !********************************************************************
     
      INTEGER(I4B)             ::  lworkl,k,n
      INTEGER(I4B)    :: maxn, maxnev, maxncv,maxbdw, ldv  , lda

      INTEGER(I4B), DIMENSION(11)                ::   iparam
      INTEGER(I4B), DIMENSION(:), ALLOCATABLE    :: s, diagval
      INTEGER(I4B), DIMENSION(:), ALLOCATABLE    ::  iwork, resid, ax, mx
      LOGICAL(LGT), DIMENSION(:), ALLOCATABLE    ::  select
      REAL(SP), DIMENSION(:,:), ALLOCATABLE       ::  rfac
      REAL(SP), DIMENSION(:), ALLOCATABLE       :: workl,val
      REAL(SP), DIMENSION(:), ALLOCATABLE       ::  workd      
      REAL(SP), DIMENSION(:,:), ALLOCATABLE     :: d  
      REAL(SP), EXTERNAL :: snrm2
      
      INTEGER(I4B)          info, jj, ido, isub2, isup2
      INTEGER(I4B)          isub, isup, idiag, nconv
      REAL(SP)             tol, h, r1, r2
      

      REAL(SP), PARAMETER ::   one = 1.0E+0, zero = 0.0E+0, two = 2.0E+0,four = 4.0E+0, six = 6.0E+0
      
      
      n = dof
      maxn = n 
      maxnev = nev
      maxncv= ncv 
      maxbdw= kl+1
      ldv = n
      lda = 2*kl+ku +1
      allocate(val(10))
      allocate(s(maxn))
      allocate(iwork(maxn))
      allocate(resid(maxn))
      allocate(ax(maxn))
      allocate(mx(maxn))
      allocate(select(maxncv))
      allocate(rfac(lda,maxn))
      allocate(workl(maxncv*maxncv+8*maxncv))
      allocate(workd(3*maxn))
      !allocate(eigenvectors(ldv, maxncv))
      allocate(d(maxncv,2))

      
     !*****************************************************************************
     ! The number N is the dimension of the matrix.  A
     ! generalized eigenvalue problem is solved        
     ! (BMAT = 'G').  NEV is the number of eigenvalues 
     ! to be approximated. The user can modify N, NEV, 
     ! NCV and WHICH to solve problems of different    
     ! sizes, and to get different parts the spectrum. 
     ! However, the following conditions must be       
     ! satisfied:                                      
     !                   N <= MAXN                     
     !                 NEV <= MAXNEV                   
     !           NEV + 1 <= NCV <= MAXNCV               
     !*****************************************************************************

      tol  = residual_norm  ! this should be very small. Zero value corresponds to default.
      IF ( n .GT. maxn ) THEN
         CALL ERROR_FEMFEKO(1,4204)
      ELSE IF ( nev .GT. maxnev ) THEN
         CALL ERROR_FEMFEKO(1,4205)
      ELSE IF ( ncv .GT. maxncv ) THEN
         CALL ERROR_FEMFEKO(1,4206)
      END IF
     
     !*****************************************************************************
     ! The work array WORKL is used in SSAUPD as           
     ! workspace.  Its dimension LWORKL is set as          
     ! illustrated below.  The parameter TOL determines    
     ! the stopping criterion. If TOL<=0, machine          
     ! precision is used.  The variable IDO is used for    
     ! reverse communication, and is initially set to 0.   
     ! Setting INFO=0 indicates that a random vector is    
     ! generated in SSAUPD to start the Arnoldi iteration. 
     !*****************************************************************************

      lworkl  = ncv**2+8*ncv
      
      ido  = 0
      info = 0
      
     !*****************************************************************************
     ! IPARAM(3) specifies the maximum number of Arnoldi 
     ! iterations allowed.  Mode 2 of SSAUPD is used     
     ! (IPARAM(7) = 2). All these options can be changed 
     ! by the user. For details see the documentation in 
     ! SSBAND.                                            
     !****************************************************************************

      iparam(3) = maxitr
      iparam(7) = mode
      
      CALL slaset('A', 2*kl+ku+1, n, zero, zero, rfac, lda)

    
     !*********************************************************************************
     ! Call SSBAND to find eigenvalues and 
     ! eigenvectors.  Eigenvalues are      
     ! returned in the first column of D.  
     ! Eigenvectors are returned in the    
     ! first NCONV (=IPARAM(5)) columns of 
     ! eigenvectors.                                  
     !********************************************************************************
     

      CALL ssband( rvec, 'A', select, d, eigenvectors, ldv, sigma, n, sb_mat, tb_mat, lda, &
     &             rfac, kl, ku, which, bmat, nev, tol, &
     &             resid, ncv, eigenvectors, ldv, iparam, workd, workl, lworkl, &
     &             iwork, info)

      IF ( info .eq. 0) THEN

         nconv = iparam(5)

        !****************************************************************************
        !Print out convergence information 
        !****************************************************************************

        
         print *, ' Convergence information***************************'
         print *, ' Number of eigenvalue requested is ', nev
         print *, ' The number of Lanczos vectors generated',' (NCV) is ', ncv
         print *, ' The number of converged Ritz values is ',nconv
         print *, ' What portion of the spectrum ', which
         print *, ' The number of Implicit Arnoldi', ' update taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence tolerance is ', tol
         print *, ' '

         
        !***************************************************************
        ! Compute the residual norm. 
        !      A*x - lambda*x    
        !***************************************************************

         DO  jj = 1, nconv
            CALL sgbmv('Notranspose', n, n, kl, ku, one, sb_mat(kl+1,1), lda, eigenvectors(1,jj), 1, zero, ax, 1)
            CALL sgbmv('Notranspose', n, n, kl, ku, one, tb_mat(kl+1,1), lda, eigenvectors(1,jj), 1, zero, mx, 1)
            CALL saxpy(n, -d(jj,1), mx, 1, ax, 1)
            d(jj,2) = snrm2(n, ax, 1)
            d(jj,2) = d(jj,2) / abs(d(jj,1))

        END DO 

         CALL smout(FILEOUT, nconv, 2, d, maxncv, -6,&
           &      'Ritz values and relative residuals')
      
         CALL smout(6, nconv, 2, d, maxncv, -6,&
           &      'Ritz values and relative residuals')
       
         nconv_eigenvalues = nconv
         eigenvalues(1:nconv) = d(1:nconv,1)
         
       
      ELSE 

        !*********************************************************************************
        ! Either convergence failed, or there 
        ! is error.  Check the documentation  
        ! for SSBAND.                         
        !*********************************************************************************
        CALL ERROR_FEMFEKO(1,4207,int1=info)

      END IF
      deallocate(val)
      deallocate(s)
      
      deallocate(iwork)
      deallocate(resid)
      deallocate(ax)
      deallocate(mx)
      deallocate(select)
      deallocate(rfac)
      deallocate(workl)
      deallocate(workd)
      deallocate(d)
END SUBROUTINE BAND_EIGENVALUES



    
