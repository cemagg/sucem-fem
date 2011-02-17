! Last changed 23 Feb 2001 DBD - see line 530 for corrected comment
!                                and line 588 for correction.
! Last changed 09 Feb 2001 MMB - FJCM updates


SUBROUTINE EIGEN_SYSMAT(timer2,timer3)
   USE bandwidth
   USE eigen_analysis_data
   USE feminterface, ONLY: EIG_MAKE_FE_STMATRICES
   USE geometry
   USE matrix
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE problem_info
   USE unit_numbers
   IMPLICIT NONE   
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs the system matrix.
!
! Note that ALL unconnected faces are treated as on the OUTSIDE
! of a PEC box; associated edges are treated as an homogenous Dirichlet BC
! (Etan=0). 
!
! The routine implements three different S and T matrix formulations:
! Lee and Mittra's original one; Savage and Peterson's first order hierarchal
! elements (actually equivalent to Lee and Mittra's) and
! Savage and Peterson's second order hierarchal elements. 
! For theoretical details of Savage and Peterson's approach see 
! S_AND_T_MAKE_HIERARCHAL. 
!
! For "2nd" order elements, there are four types of d.o.f., 
! two associated with edges (e1 & e2) and two with
! faces (f1 & f2). Following the notation of Savage and Peterson, these are flagged 
! as follows in the code, uses the variable dof_type to keep track of these:
!
! Dof  Flag
! e1   1 
! e2   2
! f1   3
! f2   4
!
! The number of d.o.f.'s is establised by a search that determines free
! edges and faces (i.e. not constrained by an homogenous Dirichlet BC in this
! case).
! 
! The d.o.f are numbered in the single loop DOF_ASSIGN. This is done 
! by global edge number, to try to minimize the matrix bandwidth.
! Two ancilliary datastructures, viz. face_edgelink and edge_facelink, are
! used to associate faces with edges, to permit this numbering procedure
! to be performed.
!
! Various types of d.o.f.'s are intemingled during this 
! process (for the LT/QN elements). A separate datastructure, dof_type,
! keeps track of the d.o.f. type for the LT/QN elements. 
! 
! Once the d.o.f. numbering has been established, 
! the system S and T matrices are assembled
! by element, calling the routine S_AND_T_MAKE_HIERARCHAL that generates the 
! S and T elemental matrices. 
!
! Then the eigenvalue problem is solved using LAPACK or APRACK routines 
! for general real-valued matrices. Note that at present this routine
! does only LOSSLESS cavity analysis, but there is no underlying theoretical 
! restriction in the FE analysis.
! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson and R.H. Geschke
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 8 April 1997 by DBD. Output data improved slightly and
! error reporting from SSYGV monitored.
!
! Minor changes 16 May 97 by DBD to debug control names.
!
! 21 May 97: Unused variable and function in USE removed by DBD
! 18 June 97: Timing data added.
!
! Major changes starting 8 Dec 1998 by DBD, finished 18 Feb 1999. Higher order
! hierarchal elements added and routine almost entirely re-written. 
!
! Minimising of bandwidth reordering implemented. Asssignment of d.o.f. changed, 
! finished 13 May 1999 
!
! Minor extension to compute eigenvectors 2 June DBD
! 
! Mid-June 1999: extensions to banded matrices RHG 
!
! June 21: minor extension to eigenvector datastructures RHG, corrections DBD.
!
! Nov 17 99: minor changes to USE statement for MIPS 7.2.1 compiler. 
! Nov 26 99: incorrect dellocation of un-allocated x_vec removed.DBD.
! 
! Jan 23 00: documentation updated. Error handling streamlined. edge_facelink
!            storage allocation "bullet-proofed".DBD.
! Feb 4    : trivial error message correction. DBD.
! Feb 22   : Subroutines E1E1,E1E2 etc. extended for complex eigenanalysis
!            Further restructuring via use of internal subprograms.
!            Linear algebra still in progress. DBD.
! 7 Mar    : Updated to use new data structures. DBD.
! 17 Apr   : Updated to include CT/LN sparse solver RHG.
! 2 May    : Error checking for ARPACK added.  ARPACK options removed. 
! 13 May   : DEBUG_DOF moved to NUMBER_DOF. DBD.
! 30 Jan 2001 : MMB. NUMBER_DOF & MESH_INFO_WRITE moved to 'femfeko.f90'.
!               Error handling reworked.
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is a list of eigenvalues of the generalized symmetric definite 
! eigenvalue problem:
!
!    [S] [x] = k_0^2 [T] [x]
!
! where the symbols have their usual FE or EM meaning and x is the vector
! of unknowns (the tangential, edge-base, fields for 1st order elements
! or more abstract edge and face based variables for the 2nd order elements).
!
! Also the matrix bandwidth is calculated.
! 
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! LAPACK routine SSYGV (plus dependencies) for solution of generalized 
! eigenvalue problem.
!
! If not pre-installed on your system, can be downloaded from
!   http://www.netlib.org/lapack/single
! Note that the BLAS routines are also required.
!
! See detailed comments in main program and 
! NOTE THE WARNINGS RE. DEFAULT PRECISION!
!
!*******************************************************************************
   INTEGER(I4B) lda,ldb,lwork                      ! Dimension info for SSYGV
   INTEGER(I4B) info                               ! Flag from SSYGV
   REAL(SP), INTENT(OUT) :: timer2, timer3         ! Timers
   INTEGER(I4B) count, count_rate , k              ! Timers

   ! data structures associated with sparse solution

   INTEGER(I4B) rowpos, browind, sortedsize,estim_nnz, rowstart, prev_estim
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: sortedcol_ind,translated_conn_edges, sortedlist
   LOGICAL(LGT) passedtest
!  Error corrected DBD 1 May 2003 REAL(I4B) ss 
   REAL(SP) ss 
!*******************************************************************************
!*** Print individual S and T matrix entries - for debugging and testing only.
!*******************************************************************************
  
   IF(DEBUG_SYSTEM_ELEMENTS) THEN
     CALL PRINT_ELEMENT_S_AND_T
   END IF

   IF (PREDICT_SPURIOUS_EIGENMODES) THEN
     CALL COUNT_SPURIOUS_EIGENMODES
   END IF
   ! DEBUG_DOF moved to within NUMBER_DOF
  
   IF ((BANDRENUM_STORE).AND.(.NOT. (SPARSE)) ) THEN  !inserted RHG 17 APR
     CALL EIG_MAKE_FE_STMATRICES(.TRUE.) ! Calculate system bandwidth, 
                                         ! without setting up the matrices 
   END IF
         
   EXECUTE_TEST: IF (EXECUTE) THEN 

     CALL EIG_MATRIX_ALLOCATE ! Allocate storage

     CALL EIG_MAKE_FE_STMATRICES(.FALSE.)
  
     CALL SYSTEM_CLOCK(count,count_rate)
     timer2 = count/count_rate ! Wall clock time for system build.

     CALL EIG_MATRIX_SOLUTION ! Solve system matrix
     
   ELSE
     WRITE(FILEOUT,'(//A,A)') '****** EIGENANALYSIS NOT PERFORMED', & 
                             '**************'                            
     WRITE(FILEOUT,'(A)') 'EXECUTE FLAG OFF'
   END IF EXECUTE_TEST
  
   RETURN
      
CONTAINS

!*******************************************************************************
SUBROUTINE COUNT_SPURIOUS_EIGENMODES
  IMPLICIT NONE
!*******************************************************************************
!*** Estimate number of spurious eigenmodes.                                ****
!*******************************************************************************
! If an estimation of the number of spurious (trivial) eigenvalues & modes
! has been requested, this is done by counting the number of free vertices
! for CT/LN elements (these are essentially the degrees of freedom of a
! Lagrangian 1st order nodal-based scheme), and by counting both free vertices
! and mid-point vertices (essentially free edges) for the LT/QN elements 
! (corresponding to a Lagrangian 2nd order nodal-based scheme). 
! See Section 7.3.2, "Iterative & Self-adaptive Finite Elements in 
! Electrogmagnetic Modelling", M. Salazar-Palma et al, Artech House, 1998.
!*******************************************************************************
  INTEGER(I4B) iedge,inode
  num_spurious_eigenmodes = 0

  ! Cycle vertices:
  DO inode = 1,num_nodes ! CT/LN and up
    IF (vertices(inode)%free) THEN 
      num_spurious_eigenmodes =  num_spurious_eigenmodes +1
    END IF
  END DO

  ! Cycle edges:
  DO iedge = 1,num_edges
    IF (edges(iedge)%free) THEN 
      SELECT CASE (edges(iedge)%order)
      CASE (1)    
        CONTINUE
      CASE (2)
        num_spurious_eigenmodes =  num_spurious_eigenmodes +1
      CASE DEFAULT
        STOP 'IE: Invalid hierarchal order in COUNT_SPURIOUS_EIGENMODES.'
      END SELECT
    END IF
  END DO

  first_eigenmode = num_spurious_eigenmodes+1

print *,'num_spurious_eigenmodes =',num_spurious_eigenmodes

END SUBROUTINE COUNT_SPURIOUS_EIGENMODES


SUBROUTINE EIG_MATRIX_ALLOCATE
  USE geometry
  USE problem_info
  USE feminterface, ONLY: MATRIX_SPARSE_ALLOCATE
  IMPLICIT NONE
!*******************************************************************************
! This routine allocates storage as required, depending on storage 
! format and real or complex eigenanalysis.
!*******************************************************************************

  IF (BANDRENUM_STORE) THEN
    IF(REAL_EIGEN_ANALYSIS) THEN
      ALLOCATE(sb_mat(2*kl+ku+1,dof))
      ALLOCATE(tb_mat(2*kl+ku+1,dof))
      ALLOCATE(eigenvalues(dof))
      sb_mat = 0.0_SP ! Array initializations.
      tb_mat = 0.0_SP
    ELSE
      CALL ERROR_FEMFEKO(1,4200)
    END IF

  ELSE IF (SPARSE) THEN
PRINT *, "before MATRIX_SPARSE_ALLOCATE"
    CALL MATRIX_SPARSE_ALLOCATE
PRINT *, "after MATRIX_SPARSE_ALLOCATE"
  ELSE
    IF(REAL_EIGEN_ANALYSIS) THEN
      ALLOCATE(s_mat(dof,dof))
      ALLOCATE(t_mat(dof,dof))
      ALLOCATE(eigenvalues(dof))
      ALLOCATE(work(3*dof-1))
      lda = dof ! Set leading dimensions for LAPACK routine for 'A' (S)
      ldb = dof ! .. and 'B' (T)
      lwork = 3*dof-1 
      ! Build (connected) system S and T matrices. 
      ! Note that only free edges and faces enter into the system.
      s_mat = 0.0_SP ! Array initializations.
      t_mat = 0.0_SP !
    ELSE 
!      ALLOCATE(s_mat_c(dof,dof))
!      ALLOCATE(t_mat_c(dof,dof))
!      ALLOCATE(eigenvalues_c(dof))
!      ALLOCATE(work_c(3*dof-1)) ! May not be needed.
      CALL ERROR_FEMFEKO(1,4201)
    END IF
  END IF


END SUBROUTINE EIG_MATRIX_ALLOCATE
!*******************************************************************************

SUBROUTINE EIG_MATRIX_SOLUTION
  use matrix
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine calls the appropriate linear algebra routines, depending
! on storage type, and real or complex analysis, and writes out the eigenvalues.
!*******************************************************************************
  INTEGER(I4B) col

  ! FJCM
  ! Was
  ! rvec = .T.
  ! is
  rvec = .TRUE.
  bmat = 'G'       ! Generalised eigenvalue problem

  LINEAR_ALG: IF ((BANDRENUM_STORE).AND.(.NOT. SPARSE)) THEN  ! Only real eigenanalysis implemented 
    !******************************************************************
    ! An iterative banded solver is used to obtain the eigenvalues, 
    ! and possibly the eigenvectors.
    ! The eigenvectors are the Ritz vectors and the eigenvalues the 
    ! Ritz values.
    ! LAPACK band storage was implemented. 
    !******************************************************************      

    
    IF (COMPUTE_EIGENMODES) THEN 
      rvec = .TRUE.  ! Ritz vectors are requested

    ELSE
      rvec = .FALSE. ! only the Ritz values are requested 
    END IF    
    !nev eigenvectors stored in columns... 
    allocate(eigenvectors(dof,ncv))  ! eigenvectors also used for temp.storage,
                                     ! even when these are not requested 
    eigenvectors = 0 

!    CALL BAND_EIGENVALUES    
    CALL SYSTEM_CLOCK(count,count_rate)
    timer3 = count/count_rate ! Wall clock time for solve.

    ! FJCM 
    ! Was   ! 
    ! WRITE(FILEOUT,'(G14.6)' ) sqrt(eigenvalues(1:nconv_eigenvalues))
    ! is (next two lines)
    WHERE (eigenvalues .LE. 0) eigenvalues=0 ! Set eignevalues equal to 0 if NaN       

    WRITE(FILEOUT,'(//A,A)') '****** CAVITY EIGENANALYSIS RESULTS', & 
                              '**************'
    WRITE(FILEOUT,'(A,I8/)') 'Root of eigenvalues (=k0) in ascending order: '
    WRITE(FILEOUT,'(G14.6)') sqrt(eigenvalues(1:nconv_eigenvalues))

    DEALLOCATE(sb_mat)
    DEALLOCATE(tb_mat)        
  ELSE IF (SPARSE) THEN
    allocate(eigenvalues(nev))
    eigenvalues = 0   
    ALLOCATE(eigenvectors(dof,ncv))  ! necessary even if eigenvectors not
                               ! requested, also used as temporary space
    eigenvectors = 0
    CALL SPARSEIG_RI
    timer3 = count/count_rate ! Wall clock time for solve.
  ELSE ! BANDED_LA ! DBD - this comment is incorrect, this is full matrix stuff.
    !********************************************************************    

    ! The LAPACK eigenvalue solver for full matrices is used
    !********************************************************************
    REAL_CMPLX: IF (REAL_EIGEN_ANALYSIS) THEN
      IF (COMPUTE_EIGENMODES) THEN 
        ! on exit of SSYGV, s_mat contains the eigenvectors
        CALL SSYGV( (1),'V','L',(dof),s_mat,(lda),t_mat,(ldb),eigenvalues, &
              work,lwork,info)
      ELSE 
        ! eigenvectors not requested
        CALL SSYGV( (1),'N','L',(dof),s_mat,(lda),t_mat,(ldb),eigenvalues, &
              work,lwork,info)
      END IF
      DEALLOCATE(t_mat) ! Deallocate BEFORE copying s_mat to eigenvectors
                      ! to save memory.
      DEALLOCATE(work)
      IF (COMPUTE_EIGENMODES) THEN 
        ALLOCATE(eigenvectors(dof,dof))
        eigenvectors(1:dof,1:dof) = s_mat(1:dof,1:dof)
      END IF
      DEALLOCATE(s_mat)    
      CALL LAPACK_DATA  ! interprets the output flag info
      IF (DEBUG_EIGENVECTORS.AND.COMPUTE_EIGENMODES) THEN
    IF (BANDRENUM_STORE) THEN
          WRITE(FILEOUT,'(A/)') '(Real) eigenvector matrix column by column' 
          DO col=1,nconv_eigenvalues
        WRITE(FILEOUT,'(20(F6.3,1X))') eigenvectors(:,col)
          END DO
    ELSE
          WRITE(FILEOUT,'(A/)') '(Real) eigenvector matrix column by column'
          DO col=1,dof
        WRITE(FILEOUT,'(20(F6.3,1X))') eigenvectors(:,col)
          END DO
    END IF
      END IF
    ELSE ! complex  
      CALL ERROR_FEMFEKO(1,4201)
      DEALLOCATE(t_mat_c) 
      DEALLOCATE(s_mat_c) 
      DEALLOCATE(work_c) 
    END IF REAL_CMPLX
  !CALL SYSTEM_CLOCK(count,count_rate)
  !timer3 = count/count_rate ! Wall clock time for solve.
  
  END IF LINEAR_ALG
  CALL SYSTEM_CLOCK(count,count_rate) !changed RHG 24 April
  timer3 = count/count_rate ! Wall clock time for solve. !changed RHG 24 April

  ! FJCM
  ! Was
  ! WRITE(FILEOUT,'(G14.6)' ) sqrt(eigenvalues) !changed RHG 24 April
  ! is (next two lines)
  WHERE (eigenvalues .LE. 0) eigenvalues=0 ! Set eignevalues equal to 0 if NaN         
  ! DBD - correction, following was 
  ! WRITE(FILEOUT,'(G14.6)' ) sqrt(eigenvalues(1:nconv_eigenvalues))
  ! and should be 
  WRITE(FILEOUT,'(G14.6)' ) sqrt(eigenvalues)

END SUBROUTINE EIG_MATRIX_SOLUTION


!*******************************************************************************
SUBROUTINE PRINT_ELEMENT_S_AND_T
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! This utility routine prints out the  S and T entries for each individual 
! element. For testing and debugging only. 
!*******************************************************************************
  INTEGER(I4B) ielem,jj
  COMPLEX(SPC), DIMENSION(6,6) ::  Se1,Te1        ! Elemental 1st order
                                                  ! FE matrices
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te        
                                                  ! Elemental FE matrices - 
                                                  ! hierarchal.
 
  WRITE(FILEOUT,'(//A)') 'Savage and Peterson formulation used for S and T.'
  WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER

  WRITE(FILEOUT,'(//A,A)') '****************** S matrix elements ',&
                         '******************'                                                     

  DO ielem=1,num_elements

    CALL S_AND_T_MAKE_HIERARCHAL(ielem,elements(ielem)%order,&
	     elements(ielem)%mixed,Se,Te)

    WRITE(FILEOUT,'(A,I8)') 'Element number: ',ielem
    SELECT CASE (elements(ielem)%order) 
    CASE(1)                    
      DO jj = 1,6
        WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
        WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
        Se(jj,1:jj)
      END DO           
    CASE(2)          
      DO jj = 1,20
        WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
    WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
        REAL(Se(jj,1:jj)) ! HACK
      END DO
    CASE DEFAULT
      STOP 'IE: Error in routine EIGEN_SYSMAT - unimplemented hierarchal order'
    END SELECT  

  END DO 

  WRITE(FILEOUT,'(//A,A)') '****************** T matrix elements ',&
                         '******************'

  DO ielem=1,num_elements

    CALL S_AND_T_MAKE_HIERARCHAL(ielem,elements(ielem)%order,&
	     elements(ielem)%mixed,Se,Te)

    WRITE(FILEOUT,'(A,I8)') 'Element number: ',ielem
    SELECT CASE (elements(ielem)%order) 
    CASE(1)                    
      DO jj = 1,6
        WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
        WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
        Te(jj,1:jj)
      END DO           
    CASE(2)          
      DO jj = 1,20
        WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
    WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
        REAL(Te(jj,1:jj)) ! HACK
      END DO
    CASE DEFAULT
      STOP 'IE: Error in routine EIGEN_SYSMAT - unimplemented hierarchal order'
    END SELECT 
  END DO

END SUBROUTINE PRINT_ELEMENT_S_AND_T   

!*******************************************************************************

SUBROUTINE LAPACK_DATA
  IMPLICIT NONE
!*******************************************************************************
! This routine interprets the output flag info , of the LAPACK routine SSYGV
!*******************************************************************************
  ! SSYGV is a FORTRAN 77 LAPACK routine for symmetric real single
  ! precision matrices.
  ! It computes the solution of the problem
  !   [S] [x] = lambda [T] [x] (it uses A and B for S and T respectively)
  ! with the first flag set to 1.
  !
  ! The 2nd flag selects computation of  eigenvalues only 'N'
  ! or 'V' to compute eigenvectors as well. In the latter case, 
  ! the eigenvectors - sorted as per the 9th parameter "eigenvalues"
  ! are returned in the 5th parameter "s_mat" by column as well.
  !   
  ! The 3rd flag indicates that it is stored in lower triangular form
  ! )('U' indicates upper).
  !
  ! On an SGi, a pre-compiled version is availalbe
  ! in the SGI Maths library (not pre-installed, must be installed).
  ! To link, use "-lcomplib.sgimath" when runnin the f90 linker.
  !
  ! Otherwise, it must be downloaded from 
  !   http://www.netlib.org/lapack
  ! and compiled using an f77 compiler. 
  ! If manually compiled, be VERY careful that SP and I4B (both 4 bytes) 
  ! as used in this code corresponds to the f77 compiler defaults! 
  ! Many - espcially  PC compilers - use 2 bytes for integers as the default.


  SELECT CASE(info)
    CASE(:-1) ! < 0
     WRITE(FILEOUT,'(A,I2,A/)') 'Argument number ', ABS(info), &
                              ' of call to SSGYV was invalid'
    CASE(0)
      IF (COMPUTE_EIGENMODES) THEN 
        WRITE(FILEOUT,'(A,/)') & 
        'Routine SSGYV exited successfully, computing eigenvalues and eigenvectors.'
      ELSE
        WRITE(FILEOUT,'(A,/)') & 
        'Routine SSGYV exited successfully, computing eigenvalues only.'
      END IF
    CASE(1:)  ! > 0
      IF (info.LE.dof) THEN 
        CALL ERROR_FEMFEKO(1,4203,int1=info)
      ELSE
        CALL ERROR_FEMFEKO(1,4208,int1=(info-dof))
     END IF
   END SELECT  
END SUBROUTINE LAPACK_DATA

! End of internal subprograms.
 
END SUBROUTINE EIGEN_SYSMAT
!*******************************************************************************


SUBROUTINE EIG_MAKE_FE_STMATRICES(bandwidth_determine)
  USE bandwidth
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE feminterface, ONLY: CONVERTCOR, LOCAL_TO_GLOBAL_INDEX_TET
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Old documentation copied:
! Note - matrix cannot be naively symmetrized, since entries are in fairly
! randomly entered in U and L parts. Entries in the corresponding
! symmetrical submatrices are entered, term by term,  as computed for  
! terms such as Se1f1 etc. (i.e. Sf1e1 in this case). 
! This is not required for sub-matrices such as Se1e1 etc (note that the 
! whole sub-matrix is presently built - for simplicity - 
! for these entries at the cost of a little computational efficiency). 
! ****
! Calculate system bandwidth, without actually setting up the matrices
!*******************************************************************************
! This routine fills the S and T matrices with the elemental, Finite Element
! contributions. It can also be called to only calculate the matrix bandwidth
! properties instead, depending on the value of the argument <bandwidth_determine>.
! Adapted from old routines. MMB 2001-10-05
! NOTE - NOT YET FULLY EXTENDED TO QT/QN ELEMENTS!
!*******************************************************************************
  LOGICAL(LGT), INTENT(IN) :: bandwidth_determine  
                              ! indicates whether the bandwidth is being determined,
                              ! or whether the matrix entries must entered into then matrix

  INTEGER(I4B) ielem,row,col,mm,nn,indpos,i_dof
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te      ! Elemental FE matrices
  
  ! Data structures to determine the matrix bandwidth:
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: row_max, row_min, row_BW
  INTEGER(I4B) :: HBW

  ! Initialize:
  IF (bandwidth_determine) THEN
    ALLOCATE(row_max(dof))
    ALLOCATE(row_min(dof))
    ALLOCATE(row_BW(dof))
    row_max = (/ (i_dof, i_dof=1,dof) /)
    row_min = row_max
  END IF

  ! Cycle through the elements to consider every elemental contribution
  ! individually:  
  ELEMENT_LOOP: DO ielem = 1,num_elements

    ! Compute elemental S and T matrices:
    IF (.NOT.bandwidth_determine) THEN
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),&
	       Se,Te)
    END IF

    ROW_LOOP: DO mm = 1,20
      row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm) ! global row
      
      COLUMN_LOOP: DO nn = 1,20
        col = LOCAL_TO_GLOBAL_INDEX_TET(ielem,nn) ! global column

        ! Check that this member of the elemental matric does contribute
        ! to the system matrices:
        CONTRIBUTE: IF ((row.GT.0).AND.(col.GT.0)) THEN 

          BW_DET: IF (bandwidth_determine) THEN
            IF (col.GT.row_max(row)) THEN
              row_max(row) = col
            ELSE IF (col.LT.row_min(row)) THEN
              row_min(row) = col
            END IF
          ELSE
            IF (REAL_EIGEN_ANALYSIS) THEN
              IF (BANDRENUM_STORE) THEN
                sb_mat(kl+ku+1+row-col,col) = sb_mat(kl+ku+1+row-col,col) + & 
                                              REAL( Se(mm,nn) )
                tb_mat(kl+ku+1+row-col,col) = tb_mat(kl+ku+1+row-col,col) + & 
                                              REAL( Te(mm,nn) )
              ELSE IF (SPARSE) THEN
                indpos = CONVERTCOR(row,col)
                SSparse(indpos) =  SSparse(indpos) + REAL( Se(mm,nn) )
                TSparse(indpos) =  TSparse(indpos) + REAL( Te(mm,nn) )   
              ELSE
                s_mat(row,col) = s_mat(row,col) + REAL( Se(mm,nn) )
                t_mat(row,col) = t_mat(row,col) + REAL( Te(mm,nn) )
              END IF
            ELSE ! Complex
              ! Band renumbered storage not implemented for complex
              ! eigenanlaysis, checked in MESHIN.
              s_mat_c(row,col) = s_mat_c(row,col) + Se(mm,nn)
              t_mat_c(row,col) = t_mat_c(row,col) + Te(mm,nn)
            END IF
          END IF BW_DET

        END IF CONTRIBUTE
      END DO COLUMN_LOOP
    END DO ROW_LOOP
  END DO ELEMENT_LOOP

  
  IF (bandwidth_determine) THEN

    IF (SPARSE) THEN
      IF (nev.GT.ncv-1) THEN
        CALL ERROR_FEMFEKO(1,4202)
      END IF
    END IF
    row_BW = 0
    DO i_dof = 1,dof
      row_BW(i_dof) = MAX(ABS(i_dof-row_max(i_dof)),ABS(i_dof-row_min(i_dof)))
    END DO
    HBW = MAXVAL(row_BW) + 1
    DEALLOCATE(row_BW)
    DEALLOCATE(row_max)
    DEALLOCATE(row_min)
    IF (BANDRENUM) THEN
      WRITE(FILEOUT,'(A)') 'Mesh renumbered'
      WRITE(FILEOUT,'(A,I8)') 'Matrix half bandwidth : ', HBW  
    ELSE   
      WRITE(FILEOUT,'(A)') 'Mesh not renumbered'
      WRITE(FILEOUT,'(A,I8)') 'Matrix half bandwidth : ', HBW  
    END IF
    ! determine the number of super- and sub diagonals
    kl = HBW - 1
    ku = kl ! these are the same in the FEM case (S and T symmetrical)

  ELSE

    DEBUG_MATRICES: IF (DEBUG_SYSTEM_MATRIX) THEN
      WRITE(FILEOUT,'(A)') 'Savage and Peterson formulation used for S and T.'
      WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER
      WRITE(FILEOUT,'(//A,A)') '****** S  matrix ', & 
                               '**************'
      DO row=1,dof
        WRITE(FILEOUT,'(A,I8)') 'row = ',row
        IF(REAL_EIGEN_ANALYSIS) THEN
          WRITE(FILEOUT,'(10(E12.5,1X))') s_mat(row,:)
        ELSE
          WRITE(FILEOUT,'(10(E12.5,1X))') s_mat_c(row,:)
        END IF
      END DO 
      WRITE(FILEOUT,'(//A,A)') '****** T  matrix ', & 
                               '**************'
      DO row=1,dof
        WRITE(FILEOUT,'(A,I8)') 'row = ',row
        IF(REAL_EIGEN_ANALYSIS) THEN
          WRITE(FILEOUT,'(10(E12.5,1X))') t_mat(row,:)       
        ELSE
          WRITE(FILEOUT,'(10(E12.5,1X))') t_mat_c(row,:)       
        END IF
      END DO
    END IF DEBUG_MATRICES
 
  END IF

END SUBROUTINE EIG_MAKE_FE_STMATRICES
!*******************************************************************************





