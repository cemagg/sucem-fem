! Calls to following routines edited out: hack to port to SGI, DBD 8 FEb 05

! DSSKYF; DSSKYS, DMATVEC_GENR




SUBROUTINE TD_SYSMAT
   USE bandwidth
   USE feminterface, ONLY: CRS_TO_SKYLINE
   USE geometry
   USE matrix
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE problem_info
   USE td_source
   USE unit_numbers
   IMPLICIT NONE   
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs the system matrix for the time domain solver.
!
! Boundary conditions should be specified on the outside boundary via the BC and AB
! cards. At present, it is assumed that an external incident plane wave is the source.
! It is imposed via the approach in Section 12.6, J Jin, "The FE in EM", 2nd edn, Wiley 2002. 
!
! This routine uses many of the standard system matrix routines, as well as the conventions 
! of the rest of the code. 
! 
! The sparse full matrix solver uses the Compaq Extended Maths Library, requiring 
! the matrix to be stored using skyline storage. The sparse iterative solver uses 
! CRS and no conversion is needed.
! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 1 April 2003 by DBD. 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the time domain field data at points as specified in the NE
! cards. 

!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! LAPACK routine SGETF2 (plus dependencies) for factorizing matrix.
!
! If not pre-installed on your system, can be downloaded from
!   http://www.netlib.org/lapack/single
! Note that the BLAS routines are also required.
!
! See detailed comments in main program and 
! NOTE THE WARNINGS RE. DEFAULT PRECISION!
!
!*******************************************************************************
  ! Full matrix LAPACK routines:
  INTEGER(I4B) lda                                ! Dimensioning info for SGETF2
  INTEGER(I4B) info                               ! Flag from ditto
  ! Sparse matrix Compaq XML routines:
  INTEGER(I4B) ierror
  INTEGER(I4B), PARAMETER :: IPARAM_LENGTH=100,RPARAM_LENGTH=100
  INTEGER(I4B), DIMENSION(IPARAM_LENGTH) ::  iparam                                
  REAL(DP), DIMENSION(RPARAM_LENGTH) ::  rparam                                

  iparam = 0      ! Array initialization
  rparam = 0.0_DP

  IF (PML%full) THEN
    CALL TD_MATRIX_FILL_EXACT_PML ! Fill FEM system matrix - FE contributions.
  ELSE 
    CALL TD_MATRIX_FILL_APPROX_PML
  END IF 
  CALL TD_ABC_ASSEMBLY ! Augment FEM system matrix with boundary integral contributions. 

  STORAGE_TYPE: IF (SPARSE) THEN
    IF (SOLVER_TYPE.EQ.0) THEN
      ! Convert the A matrix from compressed row storage to (double precision) skyline storage
      CALL CRS_TO_SKYLINE 
      ! Note that Asparse matrix is not used again after this, and can be de-allocated in production code
	  ! (retained at present for testing). 

      ! Now factor matrix, after allocated storage and initializign parameters
      ALLOCATE(rwrk(1)) ! Dummy variable, unused.
      ALLOCATE(iwrk(2*dof)) ! Work array - required UNCHANGED for subsequent substitution stages!
  	  iparam(1) = IPARAM_LENGTH ! Length of IPARAM
      iparam(2) = RPARAM_LENGTH ! Length of RPARAM
      iparam(3) = 2*dof         ! Size of int. work array
      iparam(4) = 0             ! Size of real array (unused)
      iparam(5) = FILEOUT       ! Unit number for errors etc.
      iparam(6) = 2             ! Detailed information reporting.
	  iparam(7) = 0             ! Use defaults otherwise
	  iparam(8) = 1             ! Profile-in storage (default)
	  ! Other integer and real parameters not required or defaults used.
      STOP 'Supply DSSKYF or equivalent'
! Hack DBD 8 FEb 05
!      CALL DSSKYF((dof),Asparse_skyline,IAUdiag,(num_nonzeros_sky),iparam,rparam,iwrk,rwrk,ierror) 
      ! Check succesful factorization:
      SELECT CASE(ierror) ! Check error status on exit.
      CASE(0)
	    WRITE(FILEOUT,'(/,1X,(A),/)') 'Routine DSSKYF (TD_SYSMAT) exited successfully.'
      CASE DEFAULT 
        CALL ERROR_FEMFEKO(1,4705,int1=info)
      END SELECT  
	ELSE ! Iterative solver of some type.
      ALLOCATE(Asparse_DP(row_ind(dof+1)-1)) 
	  Asparse_DP = DBLE(Asparse)
	  ! Subsequently, set up preconditioners etc. here.
	END IF
  ELSE 
    IF (SOLVER_TYPE.EQ.0) THEN
    ! Factor the LHS matrix in [eqn.(12.33),Jin2]
    !********************************************************************    
    ! The LAPACK factorizer for full, real-valued matrices is used
    !********************************************************************
      lda = dof ! Set leading dimensions for LAPACK routine for [A] 

      CALL SGETRF((dof),(dof),A_mat,lda,ipiv,info)
      ! This is the LAPACK routine for real single precision matrices.
      ! It factors the matrix A_mat. 
      ! On exit, A_mat is overwritten by 
      ! P*L*U. (Diag of L not stored). ipiv stored pivot indices. 
      SELECT CASE(info) ! Check error status on exit.
      CASE(:-1) ! < 0
        CALL ERROR_FEMFEKO(1,4701,int1=ABS(info))
      CASE(0)
	    WRITE(FILEOUT,'(/,(A),/)') 'Routine SGETRF (TD_SYSMAT) exited successfully.'
      CASE(1:)  ! > 0
        CALL ERROR_FEMFEKO(1,4702,int1=info)
      END SELECT  
    ELSE
	  STOP 'Iterative full matrix solvers not implemented for TD analysis'
	END IF
  END IF STORAGE_TYPE

CONTAINS

!*******************************************************************************

SUBROUTINE TD_MATRIX_FILL_EXACT_PML
  USE feminterface, ONLY: LOCAL_TO_GLOBAl_INDEX_TET, CONVERTCOR
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE basis_function
  USE frequency_data
  USE problem_info
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! This routine fills the FEM system matrices, by calling a number of 
! routines that first compute, and then assemble, the elemental matrices.
!*******************************************************************************
! Note that this routine calls the general purpose FREQUENCY DOMAIN 
! routines. These return potentially complex valued elemental matrices.
! For time-domain analyis, the real part alone is retained (note that
! the matrices cannot be complex-valued in the time-domain).
!
! Note that for time-domain analysis, the "imaginary" part of eps_r stores
! the electrical conductivity. (This is for coding convenience.)
!
! See subroutine PML_ASSIGN_PROPERTIES for definitions of "zones". 
!
! This routine uses the theory in:
! [Jiao03] D.Jiao, J-M Jin, E Michielssen, D J Riley, 
! "Time-domain finite-element simulation of three-dimensional scattering and radiation
!  problems using perfectly matched layers", IEEE T-AP, Feb 2003. 
!
! Written   DBD 12 July 2003. 
!*******************************************************************************
  INTEGER(I4B) ielem,jface,row,column,col,line,linelength,mm,nn
  INTEGER(I4B) m,n,itemp
  LOGICAL(LGT) matrix_symmetry
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te ! Elemental FE matrices
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se_r,Te_r ! Elemental FE matrices 
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Te_p,Te_q,Te_x,Te_y,Te_z ! Matrices representing PML absorber. See [Jiao03]
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se_x,Se_y,Se_z ! Matrices representing PML absorber. See [Jiao03]
  REAL(SP) :: rtemp1,rtemp2,rtemp3,C1
  REAL(SP) :: epsilon,sigma,mu
  REAL(SP), DIMENSION(3,3) :: L_p, L_q, L_x, L_y, L_z, L_sx, L_sy, L_sz
  REAL(SP) sigma_x, sigma_y, sigma_z
  INTEGER(I4B) zone_label
  REAL(SP), DIMENSION(3) :: elem_centre
  REAL(SP) d, dist_x, dist_y, dist_z


  ! (Re-)initialize arrays and matrix. Unnecessary at present.
  IF (SPARSE) THEN
    Asparse = 0.0_SP ! Array initializations.
    Bsparse = 0.0_SP 
    Csparse = 0.0_SP 
  ELSE 
    A_mat = 0.0_SP
    B_mat = 0.0_SP
    C_mat = 0.0_SP
  END IF

  ! Add the FE contribution to the system matrices:

  FEM_contrib: DO ielem = 1,num_elements

    ! Compute elemental S and T matrices
    CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.) ! unscaled, but complex
	Se_r = REAL(Se)
    Te_r = REAL(Te)
	
!    if (ielem.eq.33) then
!	continue; ! for debugging
!	end if

    PML_TEST: IF(material_type(elements(ielem)%material).EQ.1) THEN ! Target or homogeneous embedding.
      epsilon = eps_0*REAL(eps_r(elements(ielem)%material))
      sigma   = AIMAG(eps_r(elements(ielem)%material))
      mu      = mu_0*REAL(mu_r(elements(ielem)%material))
    ELSE IF(material_type(elements(ielem)%material).EQ.3) THEN ! PML
      zone_label = elements(ielem)%material - MAX_MATERIALS
      epsilon = eps_r(HOMOG_MEDIUM)*eps_0
	  mu = mu_r(HOMOG_MEDIUM)*mu_0
	  elem_centre = ELEMENT_CENTRE(ielem)
	  d = PML%thickness
	  dist_x = MIN(ABS(elem_centre(1)-PML%x_min ),ABS(elem_centre(1)-PML%x_plus))
      dist_y = MIN(ABS(elem_centre(2)-PML%y_min ),ABS(elem_centre(2)-PML%y_plus))
      dist_z = MIN(ABS(elem_centre(3)-PML%z_min ),ABS(elem_centre(3)-PML%z_plus))
	  ! Note: values for sigma_x,y and z are computed irrespective of whether the points lie inside or
	  ! outside the PML layer. However, they are only USED if the element's centre lie in the appropriate PML zones, 
	  sigma_x = 0
	  sigma_y = 0
	  sigma_z = 0
	  IF (PML%absorb_x) sigma_x = PML%sigma_max * (dist_x/d)**PML%m
 	  IF (PML%absorb_y) sigma_y = PML%sigma_max * (dist_y/d)**PML%m
	  IF (PML%absorb_z) sigma_z = PML%sigma_max * (dist_z/d)**PML%m

      ! The formulation fails if the sigmas are identical (and non-zero), and will be unreliable if they are too close. 
	  ! This issue is not discussed in the literature.
	  ! A work-around used here is to simply slightly change the values in an ad-hoc fashion to avoid this. 

      IF (ABS(sigma_x).GT.EPS) THEN
        IF (ABS((sigma_x-sigma_y)/sigma_x).LT.100.0*EPS) sigma_y = 0.99_SP*sigma_y
        IF (ABS((sigma_x-sigma_z)/sigma_x).LT.100.0*EPS) sigma_z = 1.01_SP*sigma_z
      END IF
      IF (ABS(sigma_y).GT.EPS) THEN
	    IF (ABS((sigma_y-sigma_z)/sigma_y).LT.100.0*EPS) sigma_y = 0.99_SP*sigma_y
      END IF

  	  SELECT CASE(zone_label)
	  ! Re-zero sigmas depending on zone. 
      CASE (1)
        sigma_y = 0
		sigma_z = 0
      CASE (2)
        sigma_x = 0
		sigma_z = 0
      CASE (3)
        sigma_x = 0
		sigma_y = 0
      CASE (4)
        sigma_z = 0
      CASE (5)
        sigma_y = 0
      CASE (6)
        sigma_y = 0
      CASE (7) 
	    CONTINUE ! All sigmas required.
	  END SELECT

	  ! Store the relevant value of sigma for subsequent re-use.
      elements(ielem)%sigma_x = sigma_x
      elements(ielem)%sigma_y = sigma_y
      elements(ielem)%sigma_z = sigma_z

	  L_p = 0.0_SP
	  L_p(1,1) = (sigma_y+sigma_z-sigma_x)/eps_0
	  L_p(2,2) = (sigma_x+sigma_z-sigma_y)/eps_0
	  L_p(3,3) = (sigma_x+sigma_y-sigma_z)/eps_0
	  L_q = 0.0_SP
	  L_q(1,1) = (sigma_x**2- sigma_x*(sigma_y+sigma_z) + sigma_y*sigma_z)/eps_0**2
	  L_q(2,2) = (sigma_y**2- sigma_y*(sigma_x+sigma_z) + sigma_x*sigma_z)/eps_0**2
	  L_q(3,3) = (sigma_z**2- sigma_z*(sigma_x+sigma_y) + sigma_x*sigma_y)/eps_0**2
	 
	  L_x = 0.0_SP
	  L_x(1,1) = - (sigma_x**2- sigma_x*(sigma_y+sigma_z) + sigma_y*sigma_z)/eps_0**2
	  L_y = 0.0_SP
	  L_y(2,2) = - (sigma_y**2- sigma_y*(sigma_x+sigma_z) + sigma_x*sigma_z)/eps_0**2
	  L_z = 0.0_SP
	  L_z(3,3) = - (sigma_z**2- sigma_z*(sigma_x+sigma_y) + sigma_x*sigma_y)/eps_0**2
	  
	  L_sx = 0.0_SP
	  L_sy = 0.0_SP
	  L_sz = 0.0_SP


  	  SELECT CASE(zone_label)
      CASE (1)
        L_sx(2,2) = (sigma_y-sigma_x)/(sigma_z-sigma_x) 
        L_sx(3,3) = (sigma_z-sigma_x)/(sigma_y-sigma_x) 
      CASE (2)
        L_sy(1,1) = (sigma_x-sigma_y)/(sigma_z-sigma_y) 
        L_sy(3,3) = (sigma_z-sigma_y)/(sigma_x-sigma_y) 
      CASE (3)
        L_sz(1,1) = (sigma_x-sigma_z)/(sigma_y-sigma_z) 
        L_sz(2,2) = (sigma_y-sigma_z)/(sigma_x-sigma_z) 
      CASE (4)
        L_sx(2,2) = (sigma_y-sigma_x)/(sigma_z-sigma_x) 
        L_sx(3,3) = (sigma_z-sigma_x)/(sigma_y-sigma_x) 
        L_sy(1,1) = (sigma_x-sigma_y)/(sigma_z-sigma_y) 
        L_sy(3,3) = (sigma_z-sigma_y)/(sigma_x-sigma_y) 
      CASE (5)
        L_sx(2,2) = (sigma_y-sigma_x)/(sigma_z-sigma_x) 
        L_sx(3,3) = (sigma_z-sigma_x)/(sigma_y-sigma_x) 
        L_sz(1,1) = (sigma_x-sigma_z)/(sigma_y-sigma_z) 
        L_sz(2,2) = (sigma_y-sigma_z)/(sigma_x-sigma_z) 
      CASE (6)
        L_sy(1,1) = (sigma_x-sigma_y)/(sigma_z-sigma_y) 
        L_sy(3,3) = (sigma_z-sigma_y)/(sigma_x-sigma_y) 
        L_sz(1,1) = (sigma_x-sigma_z)/(sigma_y-sigma_z) 
        L_sz(2,2) = (sigma_y-sigma_z)/(sigma_x-sigma_z) 
      CASE (7)
        L_sx(2,2) = (sigma_y-sigma_x)/(sigma_z-sigma_x) 
        L_sx(3,3) = (sigma_z-sigma_x)/(sigma_y-sigma_x) 
        L_sy(1,1) = (sigma_x-sigma_y)/(sigma_z-sigma_y) 
        L_sy(3,3) = (sigma_z-sigma_y)/(sigma_x-sigma_y) 
        L_sz(1,1) = (sigma_x-sigma_z)/(sigma_y-sigma_z) 
        L_sz(2,2) = (sigma_y-sigma_z)/(sigma_x-sigma_z) 
      CASE DEFAULT
        STOP 'IE in TD_MATRIX_FILL_APPROX: this should not happen.' 
      END SELECT
	  
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,L_p) 
      Te_p = REAL(Te)

      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,L_q)
      Te_q = REAL(Te)
	  
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,L_x) 
      Te_x = REAL(Te)
      
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,L_y)
      Te_y = REAL(Te)
      
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,L_z)
      Te_z = REAL(Te)

      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,-L_sx) 
      Se_x = REAL(Se)
      
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,-L_sy)
      Se_y = REAL(Se)
      
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,-L_sz)
      Se_z = REAL(Se)


! debugging code only! zero everything....
!te_p= 0 ;  te_q = 0; te_x=0; te_y=0; te_z=0
!se_x=0; se_y=0; se_z=0


    ELSE
      STOP 'IE: Unimplemented material type in TD_MATRIX_FILL.'
    END IF PML_TEST
    ! Cycle through all values of Se&Te (in PML if appropriate), add to system matrix according to dof number and
    ! the hierarchal order used. 
    ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE

      row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)

      COL_LOOP: DO nn = 1,ELEM_TET_MATRIX_SIZE

        column = LOCAL_TO_GLOBAL_INDEX_TET(ielem,nn)

        ! If both edges/faces are free and exist, then this value must be added to the system matrix.
		! See note above. "Imaginary" part of eps_r stores sigma for TD analysis.


        IF ((row.GT.0).AND.(column.GT.0)) THEN

		  IF (material_type(elements(ielem)%material).EQ.1) THEN 
		    ! Usual isotropic media
            rtemp1 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) + & 
    		         (1.0_SP/(2.0_SP*delta_t)) *sigma*Te_r(mm,nn) + & 
                     Newmark_beta*Se_r(mm,nn)/mu
		    rtemp2 = (2.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP - 2.0_SP*Newmark_beta)*Se_r(mm,nn)/mu
            rtemp3 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP/(2.0_SP*delta_t))*sigma* Te_r(mm,nn) + & 
                     Newmark_beta*Se_r(mm,nn)/mu
          ELSE IF(material_type(elements(ielem)%material).EQ.3) THEN
		    ! PML absorber
            C1 = ( Se_r(mm,nn) + & 
			    Se_x(mm,nn)*(1-EXP(-sigma_x*delta_t/EPS_0))+Se_y(mm,nn)*(1-EXP(-sigma_y*delta_t/EPS_0))+ &
				Se_z(mm,nn)*(1-EXP(-sigma_z*delta_t/EPS_0)) )/mu + & 
				epsilon*( Te_q(mm,nn) + & 
			    Te_x(mm,nn)*(1-EXP(-sigma_x*delta_t/EPS_0))+Te_y(mm,nn)*(1-EXP(-sigma_y*delta_t/EPS_0))+ & 
				Te_z(mm,nn)*(1-EXP(-sigma_z*delta_t/EPS_0)) )  

! Test for debugging.... above "correct". When exp damping removed, goes unstable faster...
!            C1 = ( Se_r(mm,nn) + & 
!			    Se_x(mm,nn)*(1)+Se_y(mm,nn)*(1)+ &
!				Se_z(mm,nn)*(1) )/mu + & 
!				epsilon*( Te_q(mm,nn) + & 
!			    Te_x(mm,nn)*(1)+Te_y(mm,nn)*(1)+ & 
!				Te_z(mm,nn)*(1) )  

            rtemp1 = (1.0_SP/delta_t**2) *epsilon* Te_r(mm,nn) + & 
    		         (1.0_SP/(2.0_SP*delta_t)) *epsilon*Te_p(mm,nn) + & 
                     Newmark_beta*C1
		    rtemp2 = (2.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP - 2.0_SP*Newmark_beta)*C1
            rtemp3 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP/(2.0_SP*delta_t))*epsilon*TE_p(mm,nn) + & 
                     Newmark_beta*C1
          ELSE
  	    	STOP 'IE: Unimplemented material type in TD_MATRIX_FILL.'
          END IF

          IF (.NOT.SPARSE) THEN
            A_mat(row,column) = A_mat(row,column) + rtemp1
            B_mat(row,column) = B_mat(row,column) + rtemp2
            C_mat(row,column) = C_mat(row,column) + rtemp3
          ELSE
            Asparse(CONVERTCOR(row,column)) = Asparse(CONVERTCOR(row,column)) + rtemp1
            Bsparse(CONVERTCOR(row,column)) = Bsparse(CONVERTCOR(row,column)) + rtemp2
            Csparse(CONVERTCOR(row,column)) = Csparse(CONVERTCOR(row,column)) + rtemp3
		  END IF
        END IF
      END DO COL_LOOP
    END DO ROW_LOOP
  END DO FEM_contrib

  
  ! Comment - next debugging section not tested properly....
  DEBUG_MATRICES: IF (DEBUG_SYSTEM_MATRIX.AND..NOT.SPARSE) THEN
    DO row=1,dof
      WRITE (FILEOUT,'(A,I4)') 'Row = ', row
      linelength = 10
      DO line = 1,CEILING(dof/REAL(linelength))
        DO col=(line-1)*linelength+1,line*linelength
          IF (col.GT.dof) THEN
            EXIT
          END IF
      WRITE(FILEOUT,*) A_mat
!          WRITE(FILEOUT,'(A,G10.4,1X,A,G10.4,A,1X)',ADVANCE='NO') & 
!              '(',REAL(A_mat_c(row,col)),'+i*',AIMAG(A_mat_c(row,col)),')'
        END DO
        WRITE (FILEOUT,'()') ! New line
      END DO
      WRITE (FILEOUT,'(/)') ! Skip a line
    END DO 
  
  END IF DEBUG_MATRICES

END SUBROUTINE TD_MATRIX_FILL_EXACT_PML
!*******************************************************************************

SUBROUTINE TD_MATRIX_FILL_APPROX_PML
  USE feminterface, ONLY: LOCAL_TO_GLOBAl_INDEX_TET, CONVERTCOR
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE basis_function
  USE frequency_data
  USE problem_info
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! This routine fills the FEM system matrices, by calling a number of 
! routines that first compute, and then assemble, the elemental matrices.
!*******************************************************************************
! Note that this routine calls the general purpose FREQUENCY DOMAIN 
! routines. These return potentially complex valued elemental matrices.
! For time-domain analyis, the real part alone is retained (note that
! the matrices cannot be complex-valued in the time-domain).
!
! Note that for time-domain analysis, the "imaginary" part of eps_r stores
! the electrical conductivity. (This is for coding convenience.)
!
!
! Written   DBD 25 March 2003. 
! Extended  DBD June 2003 to include an approximate PML treatment. 
!*******************************************************************************
  INTEGER(I4B) ielem,jface,row,column,col,line,linelength,mm,nn
  INTEGER(I4B) m,n,itemp
  LOGICAL(LGT) matrix_symmetry
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te ! Elemental FE matrices
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se_r,Te_r ! Elemental FE matrices 
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Te_r_PML1,Te_r_PML2 ! Ditto for PML materials. 
  REAL(SP) :: rtemp1,rtemp2,rtemp3
  REAL(SP) :: epsilon,sigma,mu
  REAL(SP), DIMENSION(3,3) :: cap_sigma,sigma_tensor,sigma_E_tensor,sigma_H_tensor
  REAL(SP) sigma_x, sigma_y, sigma_z
  INTEGER(I4B) zone_label
  REAL(SP), DIMENSION(3) :: elem_centre
  REAL(SP) d, dist_x, dist_y, dist_z


  ! (Re-)initialize arrays and matrix. Unnecessary at present.
  IF (SPARSE) THEN
    Asparse = 0.0_SP ! Array initializations.
    Bsparse = 0.0_SP 
    Csparse = 0.0_SP 
  ELSE 
    A_mat = 0.0_SP
    B_mat = 0.0_SP
    C_mat = 0.0_SP
  END IF

  ! Add the FE contribution to the system matrices:

  FEM_contrib: DO ielem = 1,num_elements

    ! Compute elemental S and T matrices
    CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.) ! unscaled, but complex
	Se_r = REAL(Se)
    Te_r = REAL(Te)
    PML_TEST: IF(material_type(elements(ielem)%material).EQ.1) THEN ! Target or homogeneous embedding.
      epsilon = eps_0*REAL(eps_r(elements(ielem)%material))
      sigma   = AIMAG(eps_r(elements(ielem)%material))
      mu      = mu_0*REAL(mu_r(elements(ielem)%material))
    ELSE IF(material_type(elements(ielem)%material).EQ.3) THEN ! PML
      zone_label = material_type(elements(ielem)%material) - MAX_MATERIALS
      epsilon = eps_r(HOMOG_MEDIUM)*eps_0
	  mu = mu_r(HOMOG_MEDIUM)*mu_0
	  sigma_E_tensor = 0.0_SP
	  elem_centre = ELEMENT_CENTRE(ielem)
	  d = PML%thickness
	  dist_x = MIN(ABS(elem_centre(1)-PML%x_min ),ABS(elem_centre(1)-PML%x_plus))
      dist_y = MIN(ABS(elem_centre(2)-PML%y_min ),ABS(elem_centre(2)-PML%y_plus))
      dist_z = MIN(ABS(elem_centre(3)-PML%z_min ),ABS(elem_centre(3)-PML%z_plus))
	  ! Note: values for sigma_x,y and z are computed irrespective of whether the points lie inside or
	  ! outside the PML layer. However, they are only USED if the element centres lie in the appropriate PML zones, 
	  IF (PML%absorb_x) sigma_x = PML%sigma_max * (dist_x/d)**PML%m
 	  IF (PML%absorb_y) sigma_y = PML%sigma_max * (dist_y/d)**PML%m
	  IF (PML%absorb_z) sigma_z = PML%sigma_max * (dist_z/d)**PML%m

  	  SELECT CASE(zone_label)
      CASE (1)
        sigma_E_tensor(1,1) = -sigma_x
        sigma_E_tensor(2,2) = +sigma_x
        sigma_E_tensor(3,3) = +sigma_x
      CASE (2)
        sigma_E_tensor(1,1) = +sigma_y
        sigma_E_tensor(2,2) = -sigma_y
        sigma_E_tensor(3,3) = +sigma_y
      CASE (3)
        sigma_E_tensor(1,1) = +sigma_z
		sigma_E_tensor(2,2) = +sigma_z
        sigma_E_tensor(3,3) = -sigma_z
      CASE (4)
        sigma_E_tensor(1,1) = -sigma_x+sigma_y
        sigma_E_tensor(2,2) = +sigma_x-sigma_y
        sigma_E_tensor(3,3) = +sigma_x+sigma_y
      CASE (5)
        sigma_E_tensor(1,1) = -sigma_x+sigma_z
        sigma_E_tensor(2,2) = +sigma_x+sigma_z
        sigma_E_tensor(3,3) = +sigma_x-sigma_z
      CASE (6)
        sigma_E_tensor(1,1) = +sigma_y+sigma_z
        sigma_E_tensor(2,2) = -sigma_y+sigma_z
        sigma_E_tensor(3,3) = +sigma_y-sigma_z
      CASE (7)
        sigma_E_tensor(1,1) = -sigma_x+sigma_y+sigma_z
        sigma_E_tensor(2,2) = +sigma_x-sigma_y+sigma_z
        sigma_E_tensor(3,3) = +sigma_x+sigma_y-sigma_z
      CASE DEFAULT
        STOP 'IE in TD_MATRIX_FILL_APPROX: this should not happen.' 
      END SELECT
	  sigma_H_tensor = mu/epsilon * sigma_E_tensor
	  sigma_tensor = sigma_E_tensor + epsilon*sigma_H_tensor/mu
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,sigma_tensor) 
      Te_r_PML1 = REAL(Te)
      sigma_tensor = MATMUL(sigma_H_tensor,sigma_E_tensor)
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.,sigma_tensor)
      Te_r_PML2 = REAL(Te)
    ELSE
      STOP 'IE: Unimplemented material type in TD_MATRIX_FILL.'
    END IF PML_TEST
    ! Cycle through all values of Se&Te (in PML if appropriate), add to system matrix according to dof number and
    ! the hierarchal order used. 
    ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE

      row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)

      COL_LOOP: DO nn = 1,ELEM_TET_MATRIX_SIZE

        column = LOCAL_TO_GLOBAL_INDEX_TET(ielem,nn)

        ! If both edges/faces are free and exist, then this value must be added to the system matrix.
		! See note above. "Imaginary" part of eps_r stores sigma for TD analysis.


        IF ((row.GT.0).AND.(column.GT.0)) THEN

		  IF (material_type(elements(ielem)%material).EQ.1) THEN 
		    ! Usual isotropic media
            rtemp1 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) + & 
    		         (1.0_SP/(2.0_SP*delta_t)) *sigma*Te_r(mm,nn) + & 
                     Newmark_beta*Se_r(mm,nn)/mu
		    rtemp2 = (2.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP - 2.0_SP*Newmark_beta)*Se_r(mm,nn)/mu
            rtemp3 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP/(2.0_SP*delta_t))*sigma* Te_r(mm,nn) + & 
                     Newmark_beta*Se_r(mm,nn)/mu
          ELSE IF(material_type(elements(ielem)%material).EQ.3) THEN
		    ! PML absorber
            rtemp1 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) + & 
    		         (1.0_SP/(2.0_SP*delta_t)) *Te_r_PML1(mm,nn) + & 
                     Newmark_beta*(Se_r(mm,nn)+Te_r_PML2(mm,nn))/mu
		    rtemp2 = (2.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP - 2.0_SP*Newmark_beta)*(Se_r(mm,nn)+Te_r_PML2(mm,nn))/mu
            rtemp3 = (1.0_SP/delta_t**2)*epsilon* Te_r(mm,nn) - & 
		             (1.0_SP/(2.0_SP*delta_t))*Te_r_PML1(mm,nn) + & 
                     Newmark_beta*(Se_r(mm,nn)+Te_r_PML2(mm,nn))/mu
          ELSE
  	    	STOP 'IE: Unimplemented material type in TD_MATRIX_FILL.'
          END IF

          IF (.NOT.SPARSE) THEN
            A_mat(row,column) = A_mat(row,column) + rtemp1
            B_mat(row,column) = B_mat(row,column) + rtemp2
            C_mat(row,column) = C_mat(row,column) + rtemp3
          ELSE
            Asparse(CONVERTCOR(row,column)) = Asparse(CONVERTCOR(row,column)) + rtemp1
            Bsparse(CONVERTCOR(row,column)) = Bsparse(CONVERTCOR(row,column)) + rtemp2
            Csparse(CONVERTCOR(row,column)) = Csparse(CONVERTCOR(row,column)) + rtemp3
		  END IF
        END IF
      END DO COL_LOOP
    END DO ROW_LOOP
  END DO FEM_contrib

  
  ! Comment - next debugging section not tested properly....
  DEBUG_MATRICES: IF (DEBUG_SYSTEM_MATRIX.AND..NOT.SPARSE) THEN
    DO row=1,dof
      WRITE (FILEOUT,'(A,I4)') 'Row = ', row
      linelength = 10
      DO line = 1,CEILING(dof/REAL(linelength))
        DO col=(line-1)*linelength+1,line*linelength
          IF (col.GT.dof) THEN
            EXIT
          END IF
      WRITE(FILEOUT,*) A_mat
!          WRITE(FILEOUT,'(A,G10.4,1X,A,G10.4,A,1X)',ADVANCE='NO') & 
!              '(',REAL(A_mat_c(row,col)),'+i*',AIMAG(A_mat_c(row,col)),')'
        END DO
        WRITE (FILEOUT,'()') ! New line
      END DO
      WRITE (FILEOUT,'(/)') ! Skip a line
    END DO 
  
  END IF DEBUG_MATRICES

END SUBROUTINE TD_MATRIX_FILL_APPROX_PML
!*******************************************************************************


SUBROUTINE TD_ABC_ASSEMBLY
  USE B_matrix, ONLY: B_MAKE_HIERARCHAL
  USE boundary_conditions
  USE feminterface, ONLY: CONVERTCOR, LOCAL_TO_GLOBAL_INDEX_TRI
  USE gw_sys, ONLY: K_Z_MN
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This subroutine assembles the contribution of the ABC's to the system matrix.
! The matrix <Bs> is the contribution of a specific ABC face to the system matrix,
! thus its rows and columns refer to local face quantities (ie. 3 edges and 1 face).
! These local face quantities must be converted to local element quantities (ie. 
! 6 edges and 4 faces) in order that they can be used to obtain the correct dof
! number from the indices renumbered. This conversion is trivial for the face
! quantity, because the local element face number of the port face is know. In the
! case of the edges, the index connecting the local tri edge numbers to local tet
! edge number is obtained simply be calling LOCAL_FACEEDGES and remembering that
! the there local tri edge numbers (1,2 and 3) correspond to the three local tet 
! edge numbers in ascending order (eg. 2, 3, 6).
!
! Note that this routine is used for both the total field and scattered field
! formulations. 
!
! This routine calls LOCAL_TO_GLOBAL_INDEX_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating between them. 
!
! Written DBD 26 March 2003. Based on GW_MAKE_PORT_AMATRIX.
! 
! References: J Jin, "The FEM in EM", 2nd edn. Wiley 2003.
!*******************************************************************************

  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) ::  Bs 
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) ::  Bs_re 
  LOGICAL(LGT) :: ABC_face_found
  INTEGER(I4B) :: ielem,face_num,jface,ABC_num
  INTEGER(I4B) :: mm,nn,row,column,indpos
  INTEGER(I4B), DIMENSION(3) :: localfaceedges
  REAL(SP), DIMENSION(3) :: normal    ! Outward directed unit normal on port.
  REAL(SP) :: Y  ! See [Jin 2nd edn,p531]

  ! Compute the elemental [B] matrix if this element has a face on an ABC. 
  ! Note that an element should NOT have more than one face on any one ABC,
  ! (this is not checked, however...)
  ! but (unlike a port) may have faces on more than one ABC (in corner regions, for instance).


  ! Cycle through all element faces and add the contributions of those 
  ! that lie on an ABC:
  ELEMENT_LOOP: DO ielem = 1,num_elements

    ! Initialize for this element:
    ABC_face_found = .FALSE.
    face_num = 0 ! default (to flag no ABC for assembly routines)

    FACE_LOOP: DO jface = 1,4
      ABC_TEST: IF (faces(elements(ielem)%faces(jface))%ABC) THEN

        ! Error check and record the port face's data:
! Check removed, it was incorrect...
!        IF (ABC_face_found) CALL ERROR_FEMFEKO(1,4704,int1=ielem)
        ABC_face_found = .TRUE.
        ABC_num = faces(elements(ielem)%faces(jface))%ABCnumber
        face_num = jface
        localfaceedges = LOCAL_FACEEDGES(face_num)
		Y = ABCs(ABC_num)%Yc

        ! Set up outward-directed normal:
        normal = ABCs(ABC_num)%normal

        ! Calculate the elemental port face matrix:
        CALL B_MAKE_HIERARCHAL(ielem,jface,normal,Bs)
		Bs_re = REAL(Bs)

        ! Now add <Bs> to the system matrix:
        DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row counter
          ! Assign the global row:
          row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,face_num,mm)
          DO nn = 1,ELEM_TRI_MATRIX_SIZE ! column counter
            ! Assign the global column:
            column = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,face_num,nn)
            ! Add this elemental matrix element to the system matrix
			! if both edges/faces are free and exist.
            IF ((row.GT.0).AND.(column.GT.0)) THEN ! 
              IF (SPARSE) THEN
                indpos = CONVERTCOR(row,column)
                Asparse(indpos) = Asparse(indpos) + Y*(1.0_SP/(2.0_SP*delta_t))*Bs_re(mm,nn)
                Csparse(indpos) = Csparse(indpos) - Y*(1.0_SP/(2.0_SP*delta_t))*Bs_re(mm,nn)
              ELSE
                A_mat(row,column)  = A_mat(row,column) + Y*(1.0_SP/(2.0_SP*delta_t))*Bs_re(mm,nn)
                C_mat(row,column)  = C_mat(row,column) - Y*(1.0_SP/(2.0_SP*delta_t))*Bs_re(mm,nn)
              END IF
            END IF
          END DO
        END DO
      END IF ABC_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE TD_ABC_ASSEMBLY
!*******************************************************************************


END  SUBROUTINE TD_SYSMAT



SUBROUTINE TD_MATRIX_ALLOCATE
  USE feminterface, ONLY: MATRIX_SPARSE_ALLOCATE
  USE geometry
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine allocates storage for the time domain wave analysis. Full and sparse 
! storage is now implemented.
! DBD. 25 March 2003.
!*******************************************************************************
  IF (BANDRENUM_STORE) THEN
    CALL ERROR_FEMFEKO(1,4503)
  END IF

  IF (SPARSE) THEN
    CALL MATRIX_SPARSE_ALLOCATE
  ELSE
    ALLOCATE(A_mat(dof,dof))
	ALLOCATE(B_mat(dof,dof))
    ALLOCATE(C_mat(dof,dof))
    ! Allocate storage for LAPACK pivots. Note - may already be allocated, if so
	! deallocate first. 
    IF (ALLOCATED(ipiv)) DEALLOCATE(ipiv)
      ALLOCATE(ipiv(dof)) 
  END IF

  ALLOCATE(b_vec(dof))

  ALLOCATE(u_nplus1(dof))
  ALLOCATE(u_n(dof))
  ALLOCATE(u_nmin1(dof))
  ALLOCATE(f_nplus1(dof))
  ALLOCATE(f_n(dof))
  ALLOCATE(f_nmin1(dof))
  IF(PML_PRESENT.AND.PML%FULL) THEN
    ALLOCATE(psi_x_nplus1(dof))
    ALLOCATE(psi_x_n(dof))
    ALLOCATE(psi_y_nplus1(dof))
    ALLOCATE(psi_y_n(dof))
    ALLOCATE(psi_z_nplus1(dof))
    ALLOCATE(psi_z_n(dof))
  END IF 
  

END SUBROUTINE TD_MATRIX_ALLOCATE
!*******************************************************************************


SUBROUTINE TD_TIMESTEP
!  INCLUDE 'CXML_INCLUDE.F90' ! Retained for future use of CXML routines. 
! Removed for SGI port
  USE feminterface, ONLY: MATVECPROD, LOCAL_TO_GLOBAL_INDEX_TET


  USE geometry
  USE matrix
  USE problem_info
  USE td_source
  USE output_error, ONLY: ERROR_FEMFEKO
  USE unit_numbers
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! This routine performs the timestepping using the Newmark Beta algorithm. 
! See J Jin, "The FE in EM", 2nd edn, Wiley 2002, p. 535.
! Information relating to the plane wave is passed  via module TD_source. 
!
!
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! First version DBD 27 March 2003. 
! Extended  DBD 27 April 2003: sparse matrix direct solver included. 
! Extended  DBD 29 April 2003: sparse matrix iterative solver added. 
! Extended  DBD 20 May 2003: scattered-only option added. 
! 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the time domain field data at points as specified in the NE
! cards. 

!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! LAPACK routine SGETF2 (plus dependencies) for factorizing matrix.
!
! If not pre-installed on your system, can be downloaded from
!   http://www.netlib.org/lapack/single
! Note that the BLAS routines are also required.
!
! See detailed comments in main program and 
! NOTE THE WARNINGS RE. DEFAULT PRECISION!
! Also:
! Compaq XML routines for sparse real-valued matrices. 
! Must be purchased from Compaq. 
!
!
!*******************************************************************************
!
!*******************************************************************************
  ! LAPACK data
  INTEGER(I4B) lda,ldb                      ! Dimensioning info for SGESV 
  INTEGER(I4B) info                         
  
  !Sparse matrix Compaq XML routines data:
  INTEGER(I4B) ierror,idum
  REAL(DP) dum_DP
  ! Note: minimum 50 required for iterative solvers, but 100 for direct.
  INTEGER(I4B), PARAMETER :: IPARAM_LENGTH=100,RPARAM_LENGTH=100
  INTEGER(I4B), DIMENSION(IPARAM_LENGTH) ::  iparam                                
  REAL(DP), DIMENSION(RPARAM_LENGTH) ::  rparam                                
  EXTERNAL MATVEC

  ! REAL(DP),     DIMENSION(:),   ALLOCATABLE :: vec_in_DP,vec_out_DP
    
  ! Other variables
  INTEGER(I4B) n_time                            
  INTEGER(I4B) row
  REAL(SP) time_taken
  INTEGER(I4B) smallest_its,largest_its,average_its

  INTEGER(I4B) ielem,mm



  IF (SPARSE) THEN
	ALLOCATE(b_vec_DP(dof))
    ALLOCATE(x_vec_DP(dof))
    ALLOCATE(temp_vec_real(dof))

! Following routines not presently used. 
!	IF(SOLVER_TYPE.EQ.0) THEN ! Set up parameters for Compaq XML direct solver.
!      iparam = 0      ! Array initialization
!      rparam = 0.0_DP 
!  	  iparam(1) = IPARAM_LENGTH ! Length of IPARAM
!      iparam(2) = RPARAM_LENGTH ! Length of RPARAM
!      iparam(3) = 2*dof         ! Size of int. work array
!      iparam(4) = 0             ! Size of real array (unused)
!      iparam(5) = FILEOUT       ! Unit number for errors etc.
!      iparam(6) = 2             ! Detailed information reporting.
!      iparam(7) = 0             ! Use defaults otherwise
!	  iparam(8) = 1             ! Profile-in storage (default)
!    ELSE ! Set up parameters for Compaq XML iterative 
!      CALL DITSOL_DEFAULTS(IPARAM,RPARAM)
!	  SELECT CASE (SOLVER_TYPE)
!	  CASE(2) ! CG solver, no preconditioning
!  	    ALLOCATE(rwrk(3*dof))
!        iparam(3) = 0           ! Size of int. work array - presumably dummy
!        iparam(4) = 3*dof       ! Size of real array as allocated above.
!        iparam(7) = 0           ! No preconditioning.
!      END SELECT 
!	  ! Set the following to required or non-default values.
!	  ! See Table 10-4 pg. 10-11 of CMXL reference guide for details. 
!  	  iparam(1) = IPARAM_LENGTH ! Length of IPARAM
!      iparam(2) = RPARAM_LENGTH ! Length of RPARAM
!      iparam(5) = FILEOUT       ! Unit number for errors etc.
!      iparam(6) = 2             ! Reasonable level of information reporting.
!      iparam(8) = 2             ! Normalized residual stopping criteria.
!	  rparam(1) = residual_norm ! As set by FM card.
!	END IF
  END IF

  WRITE(FILEOUT,'(//,16X,(A))')             'TIME DOMAIN FINITE ELEMENT ANALYSIS'
  WRITE(FILEOUT,'(/,1X,(A),T26,(A),G12.4))')  'NEWMARK METHOD USED','BETA=',newmark_beta
  WRITE(FILEOUT,'(1X,(A),T26,(A),G12.4)')   'Time step','DELTA_T=',delta_t
  WRITE(FILEOUT,'(1X,(A),T26,(A),I8)')      'Num. steps','NUM_TIMESTEPS=',num_timesteps


  ! Initialize vectors
  u_nmin1 = 0.0_SP
  u_n     = 0.0_SP
  f_nmin1 = 0.0_SP
  f_n     = 0.0_SP
  IF (PML_PRESENT.AND.PML%full) THEN
	psi_x_n = 0.0_SP
	psi_y_n = 0.0_SP
	psi_z_n = 0.0_SP
	psi_x_nplus1 = 0.0_SP
    psi_y_nplus1 = 0.0_SP
	psi_z_nplus1 = 0.0_SP
  END IF

  CALL TD_INC_FIELD_SETUP ! (t_pw,p_pw,e_pw,mag_pw)
  TIMESTEPPING: DO n_time = 1,num_timesteps
    timestep = n_time+1 ! Global parameter to communicate with FE assembly routines.
    ! Compute new forcing function (f_n_plus1)
	CALL TD_INCIDENT_SOURCE ! returns f_nplus1
    ! Solve for u_nplus1 using Newmark scheme - eqn. (12.33) [Jin2nd]

	! Now solve the linear system [A] {u_n+1} = {b_vec}
	! using the already factored matrix - either sparse or full. 
	! Iterative solvers not implemented at present.

    matrix_storage: IF (SPARSE) THEN
      ! Compute the new b vector
      b_vec = 0.0_SP
      CALL MATVECPROD((7),u_n,temp_vec_real,(dof)) 
      b_vec = b_vec+temp_vec_real
      CALL MATVECPROD((8),u_nmin1,temp_vec_real,(dof)) 
	  b_vec = b_vec-temp_vec_real
	  b_vec = b_vec - (Newmark_beta*f_nplus1 & 
	         +(1-2.0_SP*Newmark_beta)*f_n + Newmark_beta*f_nmin1) 
      b_vec_DP = DBLE(b_vec)

      sparse_matrix_solver: IF (SOLVER_TYPE.EQ.0) THEN


        !********************************************************************    
        ! The Compaq XML for sparse symmetric real-valued matrices is used - 
	    ! skyline storage, profile-in. Matrix is already factored; 
        ! iwrk must remain unchanged since factorization. Parameters set 
		! above at start of routine. 
        !********************************************************************

	    ! Other integer and real parameters not required or defaults used.
        STOP ' Supply DSSKYS or equilvalent. DBD hack 8 Feb 05'
!        CALL DSSKYS((dof),Asparse_skyline,IAUdiag,(num_nonzeros_sky),b_vec_DP,(dof),(1),iparam,rparam,iwrk,rwrk,ierror) 

 ! write(fileout,*)  '[x]: 'b_vec_DP

	    u_nplus1 = REAL(b_vec_DP) 
 
        IF (DEBUG_SYSTEM_MATRIX) THEN
          ! WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER
          WRITE(FILEOUT,'(//A,I6,A)') '****** [u_n+1] vector following XML DSSKYS solution at timestep ',n_time, & 
                         '**************'
          WRITE(FILEOUT,*) u_nplus1
        END IF 
 
        SELECT CASE(ierror) ! Check error status on exit.
        CASE(0)
          IF(ON_SCREEN_REPORTING) PRINT *,'Routine DSSKYS within TD_TIMESTEP exited successfully at timestep ',n_time
        CASE DEFAULT 
          CALL ERROR_FEMFEKO(1,4706,int1=info)
        END SELECT   	
	  ELSE 
!	    x_vec_DP = 0.0_DP ! Initial guess - on return, the solution.

        ! The CXML routine DITSOL_PCG gave a number of problems. The interface statements in 
		! Program Files\Microsoft Visual Studio\DF98\CXML\INCLUDE
		! declare as integer arguments that should be integer arrays (eg the 6th). 
		! MATVEC also does not correctly receive the vector to be multiplied.
		! This problem may be a library problem, and a local equivalent has been written
		! to provide similar functionality.
		
        ! CALL DITSOL_PCG(MATVEC,PCONDL,PCONDR,MSTOP,dum_dp,idum,x_vec_DP,b_vec_DP,(dof),dum_DP,idum,dum_DP,idum,&
		!                iparam,rparam,idum,rwrk,ierror)

!write(fileout,*)  '[b]: ',b_vec_DP
        CALL ITER_SOLVE_DP(time_taken,smallest_its,largest_its,average_its)
	    u_nplus1 = REAL(x_vec_DP) 
!write(fileout,*)  '[u_nplus1]',u_nplus1

	  ENDIF sparse_matrix_solver 
    ELSE 
	  full_matrix_solver: IF (SOLVER_TYPE.EQ.0) THEN
        b_vec = MATMUL(B_mat,u_n) - MATMUL(C_mat,u_nmin1) - (Newmark_beta*f_nplus1 & 
            +(1-2.0_SP*Newmark_beta)*f_n + Newmark_beta*f_nmin1) 
        !********************************************************************    
        ! The LAPACK solver for full, real-valued matrices is used
        !********************************************************************
        lda = dof ! Set leading dimensions for LAPACK routine for [A] 
        ldb = dof ! Set leading dimensions for LAPACK routine for  {b_vec}

        CALL SGETRS('N',(dof),(1),A_mat,lda,ipiv,b_vec,ldb,info)
        ! This is the LAPACK routine for real single precision matrices.
        ! It solves the linear system [A]{x} = {b}.
        ! On exit, {b_vec} is overwritten by the solution {x} 
        u_nplus1 = b_vec 
        SELECT CASE(info) ! Check error status on exit.
        CASE(:-1) ! < 0
          CALL ERROR_FEMFEKO(1,4703,int1=ABS(info))
        CASE(0)
          IF(ON_SCREEN_REPORTING) PRINT *,'Routine SGETRS within TD_TIMESTEP exited successfully at timestep ',n_time
        END SELECT  

        IF (DEBUG_SYSTEM_MATRIX) THEN
          ! WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER
          WRITE(FILEOUT,'(//A,I6,A)') '****** [u_n+1] vector following LAPACK SGETRS solution at timestep ',n_time, & 
                         '**************'
          WRITE(FILEOUT,*) u_nplus1
        END IF 
 	  ELSE 
	    STOP 'Iterative routines for full matrix storage not implemented in TD_TIMESTEP'
      END IF full_matrix_solver
    END IF matrix_storage

    ! Write near field data
	CALL TD_OUTPUT_FE_RESULTS


	! update vectors before next time step
    IF (PML_PRESENT.AND.PML%full) THEN
	  ! Update material time history.
      ELEMENT_LOOP: DO ielem = 1,num_elements
        IF (material_type(elements(ielem)%material).EQ.3) THEN 
          ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE 
            row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)
			IF (row.GT.0) THEN
              psi_x_nplus1(row) = EXP(-elements(ielem)%sigma_x*delta_t/EPS_0) *psi_x_n(row) + &
    		                     (1-EXP(-elements(ielem)%sigma_x*delta_t/EPS_0))*u_nplus1(row)
              psi_y_nplus1(row) = EXP(-elements(ielem)%sigma_y*delta_t/EPS_0) *psi_y_n(row) + &
	    	                     (1-EXP(-elements(ielem)%sigma_y*delta_t/EPS_0))*u_nplus1(row)
              psi_z_nplus1(row) = EXP(-elements(ielem)%sigma_z*delta_t/EPS_0) *psi_z_n(row) + &
		                         (1-EXP(-elements(ielem)%sigma_z*delta_t/EPS_0))*u_nplus1(row)
            END IF
		  END DO ROW_LOOP
	    END IF
      END DO ELEMENT_LOOP
	  psi_x_n = psi_x_nplus1
	  psi_y_n = psi_y_nplus1
	  psi_z_n = psi_z_nplus1
!write (fileout, *) 'psi_x_n, time step n=',n_time
!write (fileout, *) psi_x_n
    END IF

    ! Solution vectors. 
    u_nmin1  = u_n
	u_n      = u_nplus1
    f_nmin1  = f_n
	f_n      = f_nplus1


  END DO TIMESTEPPING
  WRITE(FILEOUT,'(//,20X,A)') 'STATISTICS OF ITERATIVE SOLVER'
  WRITE(FILEOUT,'(/,1X,A,2X,I12)') 'Smallest number of iterations required:         ',smallest_its
  WRITE(FILEOUT,'(/,1X,A,2X,I12)') 'Largest number of iterations required:          ',largest_its
  WRITE(FILEOUT,'(/,1X,A,2X,I12)') 'Average number of iterations required:          ',average_its

CONTAINS 

SUBROUTINE USER_MATVEC(job,iparam,rparam,dum_DP,idum,w_DP,u_DP,v_DP,n)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This is a dummy driver routine which provides an interface to the CXML CRS-matrix - vector routine
! Note that is is hardwired to the Asparse_DP matrix, with row and column indices
! row_ind and col_ind
!
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: job,n
  INTEGER(I4B), DIMENSION(*) ::  iparam                                
  REAL(DP), DIMENSION(*) ::  rparam 
  ! Following threer 
  REAL(DP), DIMENSION(n),INTENT(IN)   ::  w_DP 
  REAL(DP), DIMENSION(n), INTENT(IN)  ::  u_DP 
  REAL(DP), DIMENSION(n), INTENT(OUT) ::  v_DP
  REAL(DP), DIMENSION(n) ::  temp
  REAL(DP), DIMENSION(*)  :: dum_DP
  INTEGER(I4B), DIMENSION(*) :: idum
  v_DP = 0.0_DP
  STOP 'USER_MATVEC is a dummy routine.'
   
END SUBROUTINE USER_MATVEC



SUBROUTINE MSTOP(iparam,rparam,vec_in1_DP,vec_in2_DP,vec_in3_DP,vec_in4_DP,n)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This is a dummy routine provided for the CXML iterative solvers. 
! See CMXL reference guide, Table 10.3
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: n
  INTEGER(I4B), DIMENSION(:), INTENT(IN) ::  iparam                                
  REAL(DP), DIMENSION(:) ::  rparam 
  REAL(DP), DIMENSION(:), INTENT(IN) ::  vec_in1_DP, vec_in2_DP, vec_in3_DP, vec_in4_DP  
  
  STOP 'MSTOP is dummy routine and should not be called.' 

END SUBROUTINE MSTOP


SUBROUTINE PCONDL 
! *******************************************************************************
! This is a dummy routine provided for the CXML iterative solvers. 
! See CMXL reference guide, Table 10.3
!*******************************************************************************
  CONTINUE
END SUBROUTINE PCONDL


SUBROUTINE PCONDR
! *******************************************************************************
! This is a dummy routine provided for the CXML iterative solvers. 
! See CMXL reference guide, Table 10.3
!*******************************************************************************
  CONTINUE
END SUBROUTINE PCONDR


SUBROUTINE TD_INCIDENT_SOURCE
  USE basis_function
  USE boundary_conditions
  USE feminterface, ONLY: U_MAKE_HIERARCHAL, &
                          LOCAL_TO_GLOBAL_INDEX_TRI,LOCAL_TO_GLOBAL_INDEX_TET
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE scattering_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! This routine computes the forcing vector f_i (eqn. 12.10, Jin 2nd edn.) for 
! an incident field source on the exterior (ABC) boundaries. It is loosely based on 
! GW_MAKE_PORT_BVECTOR. 
!
! Note that the time step is set via the global parameter timestep, communicated via module TD_source.
! 
! This routine calls LOCAL_TO_GLOBAL_INDEX_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating them. 
!
! Note also that using the surface basis functions for the surface ingegration 
! brings an accompanying minus sign,
! taken into account with the waveguide analysis in eqn.(8.85) [Jin2nd}
! which then becomes positive when taken to the RHS in (8.89) [ibid]
! and can thus be neglected in that case, but must be taken into 
! account with TD analysis total field formulation in (12.3) [ibid] when constructing the f vector.
! Hence the minus sign in this routine in those cases. 
! This is not relevant to the scattered field formulation of (12.80). 
!
! For the scattered-field analysis, the following is assumed:
! 1. The inhomogeneous scatterer is the interior region. It may be sited anywhere within the mesh, 
!    but there MUST be at least one layer of tets with HOMOG_MEDIUM properties between the scatterer 
!    and the ABC. 
! 2. Materials with label HOMOG_MEDIUM comprise the exterior homogeneous region. 
! 3. The scatterer consists of materials with labels other than this. Note that 
!    these may include materials with the same properties as label this, but 
!    the incident field is also computed in these regions. 
!
! Note that in this case, the incident field is computed throughoutthe region 
! occupied by the inhomogneous scatterer, instead of entering via only a boundary 
! integral. 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Written DBD 26 March 2003. 
! Extended  DBD 20 May 2003: scattered-only option added. 
! 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the forcing function for the Newmark algorithm. 
!*******************************************************************************

  INTEGER(I4B) ielem,iface,row,ABC_num,idof,dummy,mm ,nn
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te ! Elemental FE matrices
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se_r,Te_r ! Elemental FE matrices 
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE) :: interpolate_coeff
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  Us 
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  Us_re 

  REAL(SP), DIMENSION(3) :: outward_normal
  REAL(SP) vol_int
  REAL(SP) temp1, temp2
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes
  INTEGER(I4B) jnode

  LOGICAL(LGT) ABC_face_found
  
  f_nplus1  = 0.0_SP ! Array re-initialization at each time step. 

  SCAT_FIELD_TEST: IF(.NOT.SCAT_FIELD) THEN
    ! Total field formulation. Incident field enters only on exterior boundary.
    ELEMENT_LOOP1: DO ielem = 1,num_elements
      ! Note that an element should NOT have more than one face on any  exterior boundary,
      ! but (unlike port analysis) may have faces on more than one exterior boundary (eg in corner regions).
      ABC_face_found = .FALSE.
      FACE_LOOP1: DO iface = 1,4
        ABC_TEST: IF (faces(elements(ielem)%faces(iface))%ABC) THEN ! This face is part of a ABC
        
          ! Record the relevant face information:
          ABC_face_found = .TRUE.
          ABC_num = faces(elements(ielem)%faces(iface))%ABCnumber

          ! Calculate the elemental contribution:
          CALL U_MAKE_HIERARCHAL(ielem,iface,MAX_ORDER(ielem,iface),ABC_num,& 
                                 ABCs(ABC_num)%normal,Us)
           ! This returns a complex number, but only the real part is needed here. 
	  	   ! (Imaginary part is zero by definition here). 
          Us_re = REAL(Us) 
          ! Add this face's contribution to <f_nplus1>:
          ROW_LOOP: DO idof = 1,ELEM_TRI_MATRIX_SIZE 
            row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,idof)
            IF (row.GT.0) THEN ! this is a dof
		      ! See comments above w.r.t. sign of Us_re
              f_nplus1(row) = f_nplus1(row) - Us_re(idof)
            END IF
          END DO ROW_LOOP
        END IF ABC_TEST
      END DO FACE_LOOP1
	END DO ELEMENT_LOOP1

  ELSE 
    ! Scattered field formulation. Incident field enters via volume integral over inhomogeneous scatterer 
	! as well as interface between inhomogeneous  and homogeneous regions.
	! Note that only istropic materials are treated here, otherwise PML region is incorrectly treated 
	! as part of the  inhomogeneous region.
    ELEMENT_LOOP2: DO ielem = 1,num_elements
      SCAT_TOTAL_TEST: IF (elements(ielem)%material.NE.HOMOG_MEDIUM.AND.material_type(elements(ielem)%material).EQ.1) THEN 
	    ! Add contributions from volume integral - this is the 2nd term in eqn. (12.80), Jin 2nd. 
	    ! Compute the usual S and T matrices
        CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,.true.) ! unscaled, but complex
        Se_r = REAL(Se)
        Te_r = REAL(Te)

	    ! Firstly, add contribution from incident field, interpolated using VBF's:
	    ! This is the first term in the integrand of the above. 
        CALL REAL_INTERPOLATE_FUNCTION(0,ielem,interpolate_coeff,(HIERARCHAL_ORDER),(MIXED_ORDER_FLAG))
        ROW_LOOP3: DO mm = 1,ELEM_TET_MATRIX_SIZE 
          row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)
          IF (row.GT.0) THEN 
            vol_int = 0.0_SP
            DO nn = 1,ELEM_TET_MATRIX_SIZE 
		      ! Note that Se, Te and interpolate_coeff are all zero for non-existent higher-ordercoefficients of a particular order,
			  !  hence no special treament is needed. Ditto other two integrals below. 
              vol_int = vol_int + interpolate_coeff(nn)*Se_r(mm,nn)/(mu_0*REAL(mu_r(elements(ielem)%material))) 
	  	    END DO
            f_nplus1(row) = f_nplus1(row) + vol_int		  
          END IF
        END DO ROW_LOOP3

        ! Secondly, add contribution from 2nd time derivative of the incident field, interpolated using VBF's:
        ! This is the second term in the integrand of the above. 
  	    CALL REAL_INTERPOLATE_FUNCTION(2,ielem,interpolate_coeff,(HIERARCHAL_ORDER),(MIXED_ORDER_FLAG))
        ! Now, compute the usual S and T matrices
        ROW_LOOP4: DO mm = 1,ELEM_TET_MATRIX_SIZE 
          row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)
          IF (row.GT.0) THEN 
    	    vol_int = 0.0_SP
            DO nn = 1,ELEM_TET_MATRIX_SIZE 
              vol_int = vol_int + interpolate_coeff(nn)*eps_0*REAL(eps_r(elements(ielem)%material))* Te_r(mm,nn) 
	        END DO
            f_nplus1(row) = f_nplus1(row) + vol_int
		  END IF
        END DO ROW_LOOP4
		
        ! Finally, add contribution from 1st time derivative of the incident field, interpolated using VBF's:
        ! This is the third term in the integrand of the above. 
        ! Note that "imaginary" part of eps_r stores sigma for TD analysis.
        CALL REAL_INTERPOLATE_FUNCTION(1,ielem,interpolate_coeff,(HIERARCHAL_ORDER),(MIXED_ORDER_FLAG))
        ! Now, compute the usual S and T matrices
        ROW_LOOP5: DO mm = 1,ELEM_TET_MATRIX_SIZE 
          row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)
          IF (row.GT.0) THEN 
  	        vol_int = 0.0_SP
            DO nn = 1,ELEM_TET_MATRIX_SIZE 
              vol_int = vol_int + interpolate_coeff(nn)* AIMAG(eps_r(elements(ielem)%material)) * Te_r(mm,nn) 
	        END DO
            f_nplus1(row) = f_nplus1(row) + vol_int
		  END IF
        END DO ROW_LOOP5

        FACE_LOOP2: DO iface = 1,4
		  FACE_TEST: IF (faces(elements(ielem)%faces(iface))%scat_tot_boundary) THEN 
		  
		    ! This face is part of S_s, the fictitious surface demarcating the homogeneous and inhomogenous regions,
		    ! and an outward directed normal to the element with also be outward normal to this demarcation interface.
		    ! Note that this test is only performed on elements inside this surface. (This also ensures that the face 
		    ! is only integrated over once, since the face is shared by two elements, one inside and one outside this surface). 
		  
    		outward_normal = ELEMENT_NORMAL_VECTOR(ielem,iface)

            CALL U_MAKE_HIERARCHAL(ielem,iface,MAX_ORDER(ielem,iface),dummy,outward_normal,Us)
            Us_re = REAL(Us)
            ! Add this face's contribution to <f_nplus1>. 
            ROW_LOOP2: DO idof = 1,ELEM_TRI_MATRIX_SIZE 
              row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,idof)
              IF (row.GT.0) THEN 
                f_nplus1(row) = f_nplus1(row) - Us_re(idof)
			  END IF
            END DO ROW_LOOP2
  		  END IF FACE_TEST
    	END DO FACE_LOOP2	  
	  END IF SCAT_TOTAL_TEST! No action needed if not on this surface.
      IF  (material_type(elements(ielem)%material).EQ.3) THEN
	    ! This is the PML region. The forcing function must be augmented by (effectively) the exponentially
		! weighted prior time history of the material.
        ROW_LOOP6: DO mm = 1,ELEM_TET_MATRIX_SIZE 
          row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)
          IF (row.GT.0) THEN 
            temp1 = elements(ielem)%Se_x*EXP(-elements(ielem)%sigma_x*delta_t/EPS_0)*psi_x_n(row)+&
                    elements(ielem)%Se_y*EXP(-elements(ielem)%sigma_y*delta_t/EPS_0)*psi_y_n(row)+&
  	        		elements(ielem)%Se_z*EXP(-elements(ielem)%sigma_z*delta_t/EPS_0)*psi_z_n(row)
            temp2 = elements(ielem)%Te_x*EXP(-elements(ielem)%sigma_x*delta_t/EPS_0)*psi_x_n(row)+&
                    elements(ielem)%Te_y*EXP(-elements(ielem)%sigma_y*delta_t/EPS_0)*psi_y_n(row)+&
			        elements(ielem)%Te_z*EXP(-elements(ielem)%sigma_z*delta_t/EPS_0)*psi_z_n(row)
            f_nplus1(row) = f_nplus1(row) + temp1 + temp2
		  END IF
        END DO ROW_LOOP6
	  END IF			 
    END DO ELEMENT_LOOP2
  
  END IF SCAT_FIELD_TEST


END SUBROUTINE TD_INCIDENT_SOURCE


!*******************************************************************************



SUBROUTINE TD_OUTPUT_FE_RESULTS
  USE feminterface, ONLY: OUTPUT_FE_RESULTS
  USE near_field_data
!*******************************************************************************
! Outputs near field results for the specified analysis.
!*******************************************************************************
  INTEGER(I4B) this_analysis      ! analysis number to be processed  
  INTEGER(I4B) fec_count          ! Counters
  LOGICAL(LGT), SAVE :: first_call = .true. 

  this_analysis = 1 ! Set analysis number (to conform with FEMFEKO convention)
!  CALL DATE_AND_TIME (values=time1)

  ! Process all cards belonging to this analysis:    
  DO fec_count = 1,NUM_FE_cards
    IF (this_analysis.EQ.FEdata(fec_count)%analysis_no)  THEN ! this card is part of the analysis
      ! Confirm reading the card to be processed:
	  IF (first_call) THEN
        WRITE (FILEOUT,'(/,1X,A)') 'read from memory:'   
        WRITE (FILEOUT,'(1X,A,5X,5(1X,I5),6(1X,E12.5))') 'FE', FEdata(fec_count)%feltyp, FEdata(fec_count)%n_x, &
                                                   FEdata(fec_count)%n_y, FEdata(fec_count)%n_z, &
                                                   FEdata(fec_count)%felkor, FEdata(fec_count)%x0, &
                                                   FEdata(fec_count)%y0, FEdata(fec_count)%z0, &
                                                   FEdata(fec_count)%delta_x, FEdata(fec_count)%delta_y, &
                                                   FEdata(fec_count)%delta_z
        first_call = .false.
	  END IF
      CALL OUTPUT_FE_RESULTS(FEdata(fec_count)%feltyp,FEdata(fec_count)%n_x,       &
                             FEdata(fec_count)%n_y, FEdata(fec_count)%n_z,         &
                             FEdata(fec_count)%felkor, FEdata(fec_count)%x0,       &
                             FEdata(fec_count)%y0, FEdata(fec_count)%z0,           &
                             FEdata(fec_count)%delta_x, FEdata(fec_count)%delta_y, &
                             FEdata(fec_count)%delta_z)
    END IF
  END DO

!  CALL DATE_AND_TIME (values=time2)
!  FE_time = FE_time + TIME_DIFFERENCE(time1,time2)

END SUBROUTINE TD_OUTPUT_FE_RESULTS
!*******************************************************************************


END SUBROUTINE TD_TIMESTEP


SUBROUTINE MATVEC(job,iparam,rparam,dum_DP,idum,w_DP,u_DP,v_DP,n)
  USE matrix
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This is a driver routine which provides an interface to the CXML CRS-matrix - vector routine
! Note that is is hardwired to the Asparse_DP matrix, with row and column indices
! row_ind and col_ind
!
! This routine is NOT working correctly! 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: job,n
  INTEGER(I4B), DIMENSION(*) ::  iparam                                
  REAL(DP), DIMENSION(*) ::  rparam 
  REAL(DP), DIMENSION(n),INTENT(IN)   ::  w_DP 
  REAL(DP), DIMENSION(n), INTENT(IN)  ::  u_DP 
  REAL(DP), DIMENSION(n), INTENT(OUT) ::  v_DP
  REAL(DP), DIMENSION(n) ::  temp
  REAL(DP), DIMENSION(*)  :: dum_DP
  INTEGER(I4B), DIMENSION(*) :: idum

  STOP 'This routine is presently not working correctly, use USER_MATVEC instead.'

  ! This is where the code is having problems - input vector u_DP is undefined. ????
  ! Replacing it with MATVECPROD, which DOES work, is problematic....

temp = w_DP
temp = u_DP
temp = v_DP



  SELECT CASE(job)
  CASE (0)
    STOP 'Supply DMATVEC_GENR or equivalent. DBD hack 8 Feb 05'
!    CALL DMATVEC_GENR(job,Asparse_DP,row_ind,col_ind,(row_ind(n+1)-1),w_DP,u_DP,v_DP,n)
  CASE DEFAULT
    STOP 'Error in MATVEC, unimplemented job type'
  END SELECT
END SUBROUTINE MATVEC

