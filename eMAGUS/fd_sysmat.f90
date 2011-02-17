SUBROUTINE FD_SCAT_SYSMAT
  USE gw_sys, ONLY: GW_MAKE_PORT_AMATRIX
  USE feminterface, ONLY: MAKE_FE_AMATRIX
  USE frequency_data
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs the system matrix for the frequency domain solver.
!
! Boundary conditions are assumed as ABC on the external boundary. 
! The source is an external incident plane wave.
! It is imposed via the approach in Section 9.3, J Jin, "The FE in EM", 2nd edn, Wiley 2002. 
!
! This routine uses many of the standard system matrix routines, as well as the conventions 
! of the rest of the code. 
! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 16 Dec 2003 by DBD. 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the frequency domain field solution.
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! None.
!
!*******************************************************************************


  ! (Re-)initialize arrays and matrix. (Re-initialization needed for 
  ! frequency loop).
  IF (SPARSE) THEN
    Asparse_c = ZERO_C ! Array initializations.
  ELSE 
    A_mat_c = ZERO_C 
  END IF
  x_vec_c = ZERO_C 

  ! Add the FE (S&T matrices) contribution to the system matrix:
  CALL MAKE_FE_AMATRIX(k0) ! This is the standard FE frequency domain assembly routine. 
                           ! and assembles the S and T elemental matrices into the global matrix.

  ! Add contribution of ABC's to the system matrix:
  CALL FD_SCAT_ABC_ASSEMBLY(k0)

!*******************************************************************************

CONTAINS


SUBROUTINE FD_SCAT_ABC_ASSEMBLY(k)
  USE B_matrix, ONLY: B_MAKE_HIERARCHAL
  USE boundary_conditions
  USE feminterface, ONLY: CONVERTCOR, LOCAL_TO_GLOBAL_INDEX_TRI
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
! Written DBD 16 December 2003. Based on TD_ABC_ASSEMBLY (which in turn was based on GW_MAKE_PORT_AMATRIX).
! Extended DBD 21 July 2004, to include a partial implementation of the 2nd order ABC,
! and extended again by DBD on 29 Sept 2004 to include the full implementation hereof.
!
! Extended DBD 27 Jan 2005 to support non-sperihcal ABC's with differing effective radii.
! 
! This routine supports both the total and scattered field formulations. 
! 
! References: J Jin, "The FEM in EM", 2nd edn. Wiley 2003.
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k  ! wavenumber in external medium 
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) ::  BsABC1, BsABC2, BsABC3 
  LOGICAL(LGT) :: ABC_face_found
  INTEGER(I4B) :: ielem,face_num,jface,ABC_num
  INTEGER(I4B) :: mm,nn,row,column,indpos
  INTEGER(I4B), DIMENSION(3) :: localfaceedges
  REAL(SP), DIMENSION(3) :: normal    ! Outward directed unit normal on boundary.
  REAL(SP) :: Y     ! See [Jin 2nd edn,p531]
  COMPLEX(SPC) :: beta  ! See [Jin 2nd edn,p354]
  REAL(SP) :: ABC_Radius ! Radius of ABC. 


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
        ABC_face_found = .TRUE.
        face_num = jface
        localfaceedges = LOCAL_FACEEDGES(face_num)
		Y = SQRT(eps_0/mu_0) ! Presently only matched to free space, real.

        ! Find outward-directed normal:
        normal = ELEMENT_NORMAL_VECTOR(ielem,face_num)

        ! Calculate the elemental face matrix for the 1st order ABC and the additional terms required by the 
		! the 2nd order ABC and add to the system matrix:
 
        ABC_Radius = ABCs(1)%r
		IF (ABC_Radius.LT.0.0_SP) THEN
		  ABC_Radius = SQRT(DOT_PRODUCT(FACE_CENTRE(ielem,jface),FACE_CENTRE(ielem,jface)))
	    END IF
        CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC1)
        IF (ABCs(1)%type.EQ.-2) THEN
          CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC2,1)
          BsABC3 = 0.0_SPC ! Ditto 
          beta = 1.0_SP/(2.0_SP*j*k + 2.0_SP/ABC_Radius) ! 1/(2jk+2/r)
        ELSE IF (ABCs(1)%type.EQ.2) THEN
          CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC2,1)
          CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC3,2)
          beta = 1.0_SP/(2.0_SP*j*k + 2.0_SP/ABC_Radius) ! 1/(2jk+2/r)
        ELSE IF (ABCs(1)%type.EQ.3) THEN ! this is experimental code....
          CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC2,1)
          CALL B_MAKE_HIERARCHAL(ielem,jface,normal,BsABC3,3)
          beta = 1.0_SP/(2.0_SP*j*k + 2.0_SP/ABC_Radius) ! 1/(2jk+2/r)
        ELSE
          BsABC2 = 0.0_SPC ! Initialize to zero, so that if the ABC is 1st order, this term can be added without effect.
          BsABC3 = 0.0_SPC ! Ditto 
		  beta = 0.0_SPC   ! Ditto
        END IF

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
                Asparse_c(indpos) = Asparse_c(indpos) + j*k*BsABC1(mm,nn) + beta*BsABC2(mm,nn) - beta*BsABC3(mm,nn) 
              ELSE
                A_mat_c(row,column)  = A_mat_c(row,column) + j*k*BsABC1(mm,nn) + beta*BsABC2(mm,nn) - beta*BsABC3(mm,nn) 
              END IF
            END IF
          END DO
        END DO


      END IF ABC_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE FD_SCAT_ABC_ASSEMBLY

!*******************************************************************************


END SUBROUTINE FD_SCAT_SYSMAT


SUBROUTINE FD_SCAT_BVECTOR_SCAT(k0)
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TET, LOCAL_TO_GLOBAL_INDEX_PRE_TET
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine incorporates the effects of prescribed degrees of freedom. The routine is very 
! similar to MAKE_FE_AMATRIX, except that only dof's corresponding to Dirichlet BC's
! contribute. 
!
! The source is an external incident plane wave.
! It is imposed via the approach in Section 9.3, J Jin, "The FE in EM", 2nd edn, Wiley 2002. 
!
! This routine uses many of the standard system matrix routines, as well as the conventions 
! of the rest of the code. 
! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 07 Aug 2004 by DBD. 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
! The Dirichlet BC's must already have been set.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the frequency domain field solution.
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! None.
!
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k0    ! wavenumber at which assembly must take place
  
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te ! Elemental FE matrices
  INTEGER(I4B) :: ielem,row,column,mm,nn

  FEM_contrib: DO ielem = 1,num_elements

    ! Compute elemental S and T matrices
    CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te) ! this incorporates epsilon and mu

    ! Cycle through all values of Se&Te, add to RHS vector according to dof number and
    ! the hierarchal order used:
    ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE

      row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm) 

      COL_LOOP: DO nn = 1,ELEM_TET_MATRIX_SIZE

        column = LOCAL_TO_GLOBAL_INDEX_PRE_TET(ielem,nn) ! This give the "column" corresponding to the Sfp etc
		                                                     ! entries. 

        ! If the dof corresponding to the row is free, and that corresponding to the column is prescribed,
		! the contribution must be added to the RHS vector and multiplied by the prescribed dof's, 
        IF ((row.GT.0).AND.(column.GT.0)) THEN
          b_vec_c(row) = b_vec_c(row) - (Se(mm,nn)+ (k0**2)*Te(mm,nn))*x_pre_vec_c(column)
        END IF
      END DO COL_LOOP
    END DO ROW_LOOP
  END DO FEM_contrib

!write (fileout,'(1X,///)')
!write (fileout,*) 'b_vec_c',b_vec_c

! test code only...
!b_vec_c = 0.5*b_vec_c ! Neither 0.5 nor 2 help here...


END SUBROUTINE FD_SCAT_BVECTOR_SCAT


!*******************************************************************************

SUBROUTINE FD_SCAT_INC_SCAT(freq)

  USE basis_function
  USE boundary_conditions
  USE feminterface, ONLY: U_MAKE_HIERARCHAL, &
                          LOCAL_TO_GLOBAL_INDEX_PRE_TRI
  USE fd_scat_source
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  IMPLICIT NONE
!*******************************************************************************
! This routine computes the forcing vector for 
! an incident field source on the scattered/total field boundary. 
! It is loosely based on TD_INCIDENT_SOURCE. However, in the scattered field
! frequency-domain formulation, the scattered field enters as Dirichlet BC's
! in the volumetric integral and a surface integral over the scatterer. This routine
! returns the last. 
! 
! This routine calls LOCAL_TO_GLOBAL_INDEX_PRE_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating between them. 
! ROUTINE NOT PRESENTLY USED.
!
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Written DBD 26 July 2004. 
! 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the forcing function. 
!*******************************************************************************

  REAL(SP), INTENT(IN) :: freq
  INTEGER(I4B) ielem,iface,row,idof,dummy
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  Us 

  REAL(SP), DIMENSION(3) :: outward_normal

  ELEMENT_LOOP: DO ielem = 1,num_elements
    ELEMENT_TEST: IF (elements(ielem)%material.EQ.HOMOG_MEDIUM) THEN 		  
	  ! Each edge and face must only contribute once to the surface integrals; here, this is implemented
	  ! by ONLY elements in the homogeneous "embedding" contributing.
	  FACE_LOOP: DO iface = 1,4
        FACE_TEST: IF (faces(elements(ielem)%faces(iface))%scat_tot_boundary) THEN 		  
	      ! This face is part of S_s, the fictitious surface demarcating the homogeneous and inhomogenous regions,
	      ! (an outward directed normal to the element actually points inwards).
          ! Note that this test is effectively only performed on elements outside this surface
		  ! since the inhomogenous region is presently assumed to be entirely PEC and the field is 
		  ! computed ONLY on this surface. 
	      outward_normal = - ELEMENT_NORMAL_VECTOR(ielem,iface)
          CALL U_MAKE_HIERARCHAL(ielem,iface,MAX_ORDER(ielem,iface),dummy,outward_normal,Us)
            ! Add this face's contribution to forcing vector, if it is a prescribed dof. 
            ROW_LOOP: DO idof = 1,ELEM_TRI_MATRIX_SIZE 
              row = LOCAL_TO_GLOBAL_INDEX_PRE_TRI(ielem,iface,idof)
              IF (row.GT.0) THEN 
                b_vec_c(row) = b_vec_c(row) - j*2*pi*freq*mu_0*Us(idof)
!write(fileout,*) 'ielem =',ielem,'  iface=',iface,'  dof=',row,'  b_vec_c(row)=',b_vec_c(row)
	          END IF
            END DO ROW_LOOP
          END IF FACE_TEST
        END DO FACE_LOOP	  
     END IF ELEMENT_TEST
   END DO ELEMENT_LOOP
 
!write(fileout,*) 'b_vec_c',b_vec_c
END SUBROUTINE FD_SCAT_INC_SCAT


!*******************************************************************************


SUBROUTINE FD_SCAT_BVECTOR_TOT(k)
  USE basis_function
  USE boundary_conditions
  USE fd_scat_source
  USE feminterface, ONLY: U_MAKE_HIERARCHAL, LOCAL_TO_GLOBAL_INDEX_TRI
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  IMPLICIT NONE
!*******************************************************************************
! Adds the contribution of the incident field at the outer ABC to the excitation vector
! for the total field formulation.
!
! This routine calls LOCAL_TO_GLOBAL_INDEX_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating them. 
!
! Based on GW_MAKE_PORT_BVECTOR. 

! Note that as in GW_MAKE_PORT_BVECTOR, using the surface basis functions for the surface ingegration 
! brings an accompanying minus sign, 
! taken into account with the waveguide analysis in eqn.(8.85) [Jin2nd}
! which then becomes positive when taken to the RHS in (8.89) [ibid]; the same happens in this case with the 
! contribution from the last term in (9.61) [ibid].
!
! Written by DBD 16 Dec 2003. 
! Extended DBD 30 Nov 2004, to include effects of 2nd order ABC. 
!
! Extended DBD 27 Jan 2005 to support non-sperihcal ABC's with differing effective radii.
!
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k  ! wavenumber in external medium 
  INTEGER(I4B) ielem,iface,row,idof
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  Us
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  UsABC1, UsABC2, UsABC3 
  COMPLEX(SPC) Uinc ! Incident field function in (9.62), [Jin 2nd]
  LOGICAL(LGT) :: ABC_face_found
  REAL(SP), DIMENSION(3) :: normal    ! Outward directed unit normal on boundary.
  INTEGER(I4B) ABC_num ! Dummy parameter, not used at present.
  COMPLEX(SPC) :: beta  ! See [Jin 2nd edn,p354]
  REAL(SP) :: ABC_Radius
  CALL FD_INC_FIELD_SETUP ! Set up field parameters.
  
  ABC_num  = 0 ! Not actually used, 

  ELEMENT_LOOP: DO ielem = 1,num_elements
    ABC_face_found = .FALSE.
    FACE_LOOP: DO iface = 1,4
      ABC_TEST: IF (faces(elements(ielem)%faces(iface))%ABC) THEN
        ! Record the relevant post face information:
        ABC_face_found = .TRUE.

        ! Find outward-directed normal:
        normal = ELEMENT_NORMAL_VECTOR(ielem,iface)

        ! Calculate the elemental contribution:
        CALL U_MAKE_HIERARCHAL(ielem,iface,MAX_ORDER(ielem,iface),ABC_num,& 
                               normal,Us)
        ! Add this face's contribution to <b_vec_c>:
        ROW_LOOP: DO idof = 1,ELEM_TRI_MATRIX_SIZE 
          row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,idof)
          IF (row.GT.0) THEN ! this is a dof
             b_vec_c(row) = b_vec_c(row) + j*k*Us(idof)
          END IF
          ABC_Radius = ABCs(1)%r
		  IF (ABC_Radius.LT.0.0_SP) THEN
		    ABC_Radius = SQRT(DOT_PRODUCT(FACE_CENTRE(ielem,iface),FACE_CENTRE(ielem,iface)))
	      END IF

! 2nd order ABC: effect of  curl-curl term and surface div term. 
          IF (ABCs(1)%type.EQ.-2) THEN
            CALL U_MAKE_ABC2(ielem,iface,normal,UsABC2,1)
            UsABC3 = 0.0_SPC ! Ditto 
            beta = 1.0_SP/(2.0_SP*j*k + 2.0_SP/ABC_Radius) ! 1/(2jk+2/r)
          ELSE IF (ABCs(1)%type.EQ.2) THEN
            CALL U_MAKE_ABC2(ielem,iface,normal,UsABC2,1)
            CALL U_MAKE_ABC2(ielem,iface,normal,UsABC3,2)
            beta = 1.0_SP/(2.0_SP*j*k + 2.0_SP/ABC_Radius) ! 1/(2jk+2/r)
          ELSE IF (ABCs(1)%type.EQ.3) THEN ! 
		   STOP 'Option not supported in FD_SCAT_BVECTOR_TOT'
		  ELSE
            UsABC2 = 0.0_SPC ! Initialize to zero, so that if the ABC is 1st order, this term can be added without effect.
            UsABC3 = 0.0_SPC ! Ditto 
		    beta = 0.0_SPC   ! Ditto
          END IF
          IF (row.GT.0) THEN ! this is a dof, add effects of 2nd order ABC as above.
             b_vec_c(row) = b_vec_c(row) + beta*UsABC2(idof) - beta*UsABC3(idof)  
          END IF
! DBD 29 Nov 04

        END DO ROW_LOOP
      END IF ABC_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP
! write (fileout,*) b_vec_c

CONTAINS 

SUBROUTINE U_MAKE_ABC2(elem,local_face,normal,Us,int_request)
  USE geometry
  USE math_tools, ONLY: CROSS_PRODUCT,VECTOR_LENGTH
  USE nrtype
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs the additional terms required for the 
! 2nd order ABC applied to the total field formulation. It is based loosely on B_MAKE_HIERARCHAL.. 
!
! Currently implemented are CT/LN, LT/LN, and LT/QN . 
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! D B Davidson.
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the element number, order and element geometry. The latter
! is passed via MODULE geometry.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the elemental matrix contributing to the RHS vector.
! 
! The numbering convention of nodes, edges and faces follows [S&P] in this
! sub-routine [Table II,S&P]. 
!
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face
  REAL(SP), DIMENSION(3), INTENT(IN) :: normal
  COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: Us
  INTEGER(I4B), INTENT(IN) :: int_request
  REAL(SP),DIMENSION(3) :: ell,u_m
  REAL(SP),DIMENSION(3) :: xi_hat,eta_hat,zeta_hat ! Local rectangular coordinate system
  INTEGER(I4B) :: num_qpts ! number of quadrature points
  INTEGER(I4B) :: iquad,iedge,ifunc,mm,nn,row,col ! counters
  REAL(SP),DIMENSION(4,4) :: vertmat
  REAL(SP),DIMENSION(4) :: r_vec,r_vec_xip,r_vec_xim,r_vec_etap,r_vec_etam,s_vec
  INTEGER(I4B),DIMENSION(3) :: facenodes
  INTEGER(I4B), DIMENSION(3) :: face_edges
  INTEGER(I4B),DIMENSION(2) :: edgenodes
  REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: vbf_values,vbf_values_xip, vbf_values_xim, & 
                                                vbf_values_etap,vbf_values_etam,vbf_tan

  COMPLEX(SPC), DIMENSION(3) :: E_inc_xip,E_inc_xim,E_inc_etap,E_inc_etam
  REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE) :: vbf_norm_values, div_vbf_values
  COMPLEX(SPC) div_E_inc,temp_sca
  COMPLEX(SPC),DIMENSION(ELEM_TET_MATRIX_SIZE) :: Us_large
  REAL(SP) delta ! Delta used for central difference approximation
 

  ! Check that normal is valid:
  IF ( ABS(SQRT(DOT_PRODUCT(normal,normal))-1.0_SP).GE.EPS) THEN
    STOP 'IE: U_MAKE_ABC2 called with invalid unit normal.'
  END IF
  ! Check that a defined operation has been requested.

  SELECT CASE(int_request) 
  CASE(1:2) ! 
    CONTINUE 
  CASE DEFAULT
	STOP 'IE: Error in U_MAKE_ABC2, undefined type of integral requested'
  END SELECT

  num_qpts = 7  ! set the number of surface quadrature points - 5th order complete
  delta = 0.02_SP*SUM(T_LENGTHS(elem))/6.0_SP ! Sets the delta used to compute the central difference to 1/50 of the 
                                                ! average edge lengths.

  ! Establish vertices matrix for repeated use:
  CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vertmat)
  facenodes = LOCAL_FACENODES(local_face)

  Us_large = 0.0 ! matrix assignment
  s_vec = 0.0 ! vector assignment
  QUAD_LOOP: DO iquad = 1,num_qpts
    
    ! Calculate the xyz coords of current quad point:
    s_vec(facenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
    r_vec = MATMUL(vertmat,s_vec)

    ! Evaluate basis functions (or curl thereof):
    SELECT CASE(int_request)
	CASE(1) ! Normal component of curl of vbf.
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),1,vbf_values) 
      ! DOT with the normal:        
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	    vbf_norm_values(ifunc) = DOT_PRODUCT(normal,vbf_values(ifunc,1:3))
      END DO
	  ! Find normal component of curl of incident field (computed via H field with correct scaling)
      temp_sca = -j*k*c_0*mu_0*DOT_PRODUCT(CMPLX(normal,0.0_SP),FD_INC_FIELD(r_vec(1),r_vec(2),r_vec(3),'H',1/Z_zero))       ! k c = omega
      ! Add to <Us_large>:
      DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
        Us_large(mm) = Us_large(mm) + &
		               quad_tri_rules(num_qpts)%rule(iquad,4) * &
		               vbf_norm_values(mm)*temp_sca
      END DO
	CASE(2) ! Surface divergence term. Surface divergence is computed numerically using 
	        ! central differences of basis functions AND incident field. 
			! Local coordinate system is xi,eta,zeta, with xi chosen
			! abitrarily in direction of first edge, and zeta aligned with normal.
      face_edges = LOCAL_FACEEDGES(local_face)
      edgenodes = LOCAL_EDGENODES(face_edges(1))
      ell(:) = EDGE_UNIT_VECTOR(edgenodes(1),edgenodes(2),elem)
      xi_hat = ell
      zeta_hat = normal      
      eta_hat = CROSS_PRODUCT(zeta_hat,xi_hat)
	  r_vec_xip(1:3) = r_vec(1:3) + xi_hat*delta   ! r vector at xi plus delta
	  r_vec_xim(1:3) = r_vec(1:3) - xi_hat*delta   ! ditto minus delta
	  r_vec_etap(1:3) = r_vec(1:3) + eta_hat*delta ! etc
	  r_vec_etam(1:3) = r_vec(1:3) - eta_hat*delta
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_xip(1:3),0,vbf_values_xip) 
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_xim(1:3),0,vbf_values_xim) 
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_etap(1:3),0,vbf_values_etap) 
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_etam(1:3),0,vbf_values_etam)
	  E_inc_xip = FD_INC_FIELD(r_vec_xip(1), r_vec_xip(2), r_vec_xip(3),'E') 
	  E_inc_xim = FD_INC_FIELD(r_vec_xim(1), r_vec_xim(2), r_vec_xim(3),'E') 
	  E_inc_etap= FD_INC_FIELD(r_vec_etap(1),r_vec_etap(2),r_vec_etap(3),'E') 
	  E_inc_etam= FD_INC_FIELD(r_vec_etam(1),r_vec_etam(2),r_vec_etam(3),'E') 
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
        ! Find surface divergence using central difference approximation for the two in-plane
    	! components
 	    div_vbf_values(ifunc) = DOT_PRODUCT((vbf_values_xip(ifunc,1:3) - vbf_values_xim(ifunc,1:3)),xi_hat) + & 
                                DOT_PRODUCT((vbf_values_etap(ifunc,1:3) - vbf_values_etam(ifunc,1:3)),eta_hat)
		div_vbf_values(ifunc) = div_vbf_values(ifunc)/(2.0_SP*delta)  
        ! CAUTION HERE: E_inc is complex, DOT_PRODUCT conjugates first arugment, hence order swopped (since normal is has no imaginary
		! component). 
        div_E_inc             = DOT_PRODUCT(CMPLX(xi_hat),(E_inc_xip(1:3) - E_inc_xim(1:3))) + & 
                                DOT_PRODUCT(CMPLX(eta_hat),(E_inc_etap(1:3) - E_inc_etam(1:3)))
        div_E_inc              = div_E_inc/(2.0_SP*delta)  

	  END DO
      ! Add to <Bs_large>:
      DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
          Us_large(mm) = Us_large(mm) + &
		                    quad_tri_rules(num_qpts)%rule(iquad,4) * &
		                    div_vbf_values(mm)*div_E_inc
      END DO
	END SELECT
  END DO QUAD_LOOP


  ! Scale by face area (final step in quadrature):
  Us_large = FACE_AREA(elem,local_face) * Us_large

  ! Extract <Us> fron <Us_large>:
    DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row loop
      row = LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face,mm)
        Us(mm) = Us_large(row)
    END DO

END SUBROUTINE U_MAKE_ABC2
!*******************************************************************************



END SUBROUTINE FD_SCAT_BVECTOR_TOT
!*******************************************************************************


SUBROUTINE FD_SCAT_MATRIX_ALLOCATE
  USE feminterface, ONLY: MATRIX_SPARSE_ALLOCATE
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine allocates storage for the guided wave analysis. Full and sparse 
! storage is now implemented.
! Original DBD. Made external - MMB 7 Feb 2001
!*******************************************************************************
  IF (BANDRENUM_STORE) THEN
    CALL ERROR_FEMFEKO(1,4503)
  END IF

  IF (SPARSE) THEN
    CALL MATRIX_SPARSE_ALLOCATE
  ELSE
    ALLOCATE(A_mat_c(dof,dof))
    ! Allocate some storage for LAPACK pivots and workspace
    ALLOCATE(ipiv(dof)) ! dof = dim in LAPACK notation
    !lwork = 5*dof ! guess block size for initialization
!    READ(INFILE,NML=LAPACKDATA)   
!    ALLOCATE(work_c(lwork) )
  END IF

  ALLOCATE(x_vec_c(dof))
  ALLOCATE(b_vec_c(dof))

END SUBROUTINE FD_SCAT_MATRIX_ALLOCATE
!*******************************************************************************


SUBROUTINE FD_SCAT_PRE_X_VEC
  USE basis_function
  USE boundary_conditions
  USE fd_scat_source
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_PRE_TET
  USE matrix
  USE problem_info
  USE geometry
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This routine computes the prescribed dof's for the frequency domain scattered field
! formulation.
! Note that this routine assumes that the formulation produces a symmetrical
! system matrix! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! D B Davidson 
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 07 Aug 2004 by DBD. 
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
! The Dirichlet BC's must already have been set.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the frequency domain field solution.
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! None.
!
!*******************************************************************************
  
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE) ::  interpolate_coeff ! Coefficients of basis functions
  INTEGER(I4B) :: ielem,row,column,mm,nn

  CALL FD_INC_FIELD_SETUP ! Set up incident field parameters.

  ALLOCATE(x_pre_vec_c(pre_dof))

  DO ielem = 1,num_elements

    ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE
      row = LOCAL_TO_GLOBAL_INDEX_PRE_TET(ielem,mm) ! 
      ! If the dof corresponding to the row (column) is prescribed, it is part of the vector of prescribed dof's, 
      IF (row.GT.0) THEN  
        CALL COMPLEX_INTERPOLATE_FUNCTION((0),ielem,interpolate_coeff,(HIERARCHAL_ORDER),(MIXED_ORDER_FLAG))
	  	x_pre_vec_c(row)  = - interpolate_coeff(mm) ! Note that in a subsequent element, the same dof overwrites
		                                            ! this - but with the same information, so there is no problem.
												    ! UNLIKE matrix assembly, adjoining elements do NOT sum contributions - 
												    ! this is simply a vector of prescribed values, one per prescribed dof. 
													! The negative sign is required since the prescribed dof's are 
													! (tangential) scattered fields, and COMPLEX_INTERPOLATE_FUNCTION
													! returns the incident fields. 
      END IF
    END DO ROW_LOOP
  END DO 
!write(fileout,*) 'x_pre_vec_c',x_pre_vec_c
END SUBROUTINE FD_SCAT_PRE_X_VEC
