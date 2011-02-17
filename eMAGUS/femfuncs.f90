! Last changed:
! 04 Mar 2002. U_MAKE_HIERARCHAL extended to QT/QN.
!              Various local to global and local to local 
!              numbering routines extended to support up to QT/QN
!              elements. DBD.
! 28 Feb 2002. DBD. Hard dimensions (8 and 20) replaced with parameters
! in several places to support higher-order and polynomial complete elements. 
! 21 Feb 2002. DBD: Routine ASSIGN_FACE_AND_EDGE_ORDERS



SUBROUTINE MAKE_FE_AMATRIX(k0)
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE feminterface, ONLY: CONVERTCOR, LOCAL_TO_GLOBAL_INDEX_TET
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! Calculates the FE contribution to the A-matrix.
! Stores to a full matrix (<A_mat_c>) or a sparse matrix (<Asparse_c>).
! This routine calls LOCAL_TO_GLOBAL_INDEX_TET, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating them. 
! MMB moved from CBAA_SYSMAT to an external routine. 2001-10-01.
! Extended to support polymial complete elements, DBD, 2002-02-28.
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k0    ! wavenumber at which assembly must take place
  
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) ::  Se,Te ! Elemental FE matrices
  INTEGER(I4B) :: ielem,row,column,mm,nn,iface
  LOGICAL(LGT) :: curvilinear_flag
  FEM_contrib: DO ielem = 1,num_elements

    ! Compute elemental S and T matrices
!    CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te) ! this incorporates epsilon and mu
    IF(ALL_CURVILINEAR) THEN
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,element_curvilinear=.TRUE.) 
	  ! this incorporates epsilon and mu
   ELSE IF(CURVILINEAR) THEN
      curvilinear_flag=.FALSE.
      DO iface = 1,4
         IF (faces(elements(ielem)%faces(iface))%curvilinear) THEN
            curvilinear_flag =.TRUE. 
            EXIT ! Found a face on the scatterer boundary. 
         END IF
! print *,'untested code in MAKE_FE_AMATRIX'
      END DO
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te,element_curvilinear=curvilinear_flag) 
   ELSE 
      CALL S_AND_T_MAKE_HIERARCHAL(ielem,MAX_ORDER(ielem),MIXED_ORDER(ielem),Se,Te) 
      ! this incorporates epsilon and mu
   END IF
    ! Cycle through all values of Se&Te, add to system matrix according to dof number and
    ! the hierarchal order used:
    ROW_LOOP: DO mm = 1,ELEM_TET_MATRIX_SIZE

      row = LOCAL_TO_GLOBAL_INDEX_TET(ielem,mm)

      COL_LOOP: DO nn = 1,ELEM_TET_MATRIX_SIZE

        column = LOCAL_TO_GLOBAL_INDEX_TET(ielem,nn)

        ! If both edges/faces are free and exist, then this value must be added to the system matrix.
        IF ((row.GT.0).AND.(column.GT.0)) THEN
          IF (.NOT.SPARSE) THEN
            A_mat_c(row,column) = A_mat_c(row,column) + Se(mm,nn) &
                                  - (k0**2)*Te(mm,nn)
          ELSE
            Asparse_c(CONVERTCOR(row,column)) = Asparse_c(CONVERTCOR(row,column)) &
                                                + Se(mm,nn) - (k0**2)*Te(mm,nn)
          END IF
        END IF

      END DO COL_LOOP
    END DO ROW_LOOP
  END DO FEM_contrib

END SUBROUTINE MAKE_FE_AMATRIX
!*******************************************************************************


FUNCTION LOCAL_TO_GLOBAL_INDEX_TET(elem,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This fuction returns the system matrix index (dof number), provider the local 
! elemental contribution matrix index.
! Returns 0 if there is no associated global matrix index.
! The range of <local> is:
! QT/QN => 1:30
! DBD 2 Mar 2002
!
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_TET

  SELECT CASE(local)
  CASE(1:6)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_e1(elements(elem)%edges(local)) ! E1 functions
  CASE(7:12)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_e2(elements(elem)%edges(local-6)) ! E2 functions
  CASE(13:16)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_f1(elements(elem)%faces(local-12)) ! F2 functions
  CASE(17:20)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_f2(elements(elem)%faces(local-16)) ! F2 functions
  CASE(21:26)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_e3(elements(elem)%edges(local-20)) ! E3 functions
  CASE(27:30)
    LOCAL_TO_GLOBAL_INDEX_TET = renumbered_f3(elements(elem)%faces(local-26)) ! F3 functions
  CASE DEFAULT
    STOP 'IE: out-of-range local value in LOCAL_TO_GLOBAL_INDEX_TET.'    
  END SELECT

END FUNCTION LOCAL_TO_GLOBAL_INDEX_TET
!*******************************************************************************

FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TET(elem,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This fuction returns the system matrix index (PRESCRIBED dof number), given the local 
! elemental contribution matrix index.
! Returns 0 if there is no associated global matrix index.
! The range of <local> is:
! QT/QN => 1:30
! DBD 08 Aug 2004
! NB! NB! This routine is for PRESCRIBED dof's. Use LOCAL_TO_GLOBAL_INDEX_TET for free dof's.
!
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_PRE_TET

  SELECT CASE(local)
  CASE(1:6)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_e1(elements(elem)%edges(local)) ! E1 functions
  CASE(7:12)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_e2(elements(elem)%edges(local-6)) ! E2 functions
  CASE(13:16)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_f1(elements(elem)%faces(local-12)) ! F2 functions
  CASE(17:20)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_f2(elements(elem)%faces(local-16)) ! F2 functions
  CASE(21:26)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_e3(elements(elem)%edges(local-20)) ! E3 functions
  CASE(27:30)
    LOCAL_TO_GLOBAL_INDEX_PRE_TET = renumbered_pre_f3(elements(elem)%faces(local-26)) ! F3 functions
  CASE DEFAULT
    STOP 'IE: out-of-range local value in LOCAL_TO_GLOBAL_INDEX_PRE_TET.'    
  END SELECT

END FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TET
!*******************************************************************************



FUNCTION LOCAL_TO_GLOBAL_INDEX_TRI(elem,local_face,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This fuction returns the system matrix index (dof number), provided the local 
! facial contribution matrix index.
! Returns 0 if there is no associated global matrix index.
! The range of <local> is:
! LT/QN => 1:12
! DBB 2 Mar 2002.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_TRI

  INTEGER(I4B), DIMENSION(3) :: faceedges

  ! Find the local edges of the face:
  faceedges = LOCAL_FACEEDGES(local_face)

  SELECT CASE(local)
  CASE(1:3)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_e1(elements(elem)%edges(faceedges(local))) ! E1 functions
  CASE (4:6)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_e2(elements(elem)%edges(faceedges(local-3))) ! E2 functions         
  CASE (7)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_f1(elements(elem)%faces(local_face)) ! F1 functions         
  CASE (8)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_f2(elements(elem)%faces(local_face)) ! F2 functions         
  CASE (9:11)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_e3(elements(elem)%edges(faceedges(local-8))) ! E3 functions         
  CASE (12)
    LOCAL_TO_GLOBAL_INDEX_TRI = renumbered_f3(elements(elem)%faces(local_face)) ! F3 functions         
  CASE DEFAULT
    STOP 'IE: out-of-range local value in LOCAL_TO_GLOBAL_INDEX_TRI.'
  END SELECT

END FUNCTION LOCAL_TO_GLOBAL_INDEX_TRI
!*******************************************************************************


FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TRI(elem,local_face,local)
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This fuction returns the system matrix index (prescribed dof number), provided the local 
! facial contribution matrix index.
! Returns 0 if there is no associated global matrix index.
! The range of <local> is:
! LT/QN => 1:12
! DBB 07 Aug 2004.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face,local
  INTEGER(I4B) :: LOCAL_TO_GLOBAL_INDEX_PRE_TRI

  INTEGER(I4B), DIMENSION(3) :: faceedges

  ! Find the local edges of the face:
  faceedges = LOCAL_FACEEDGES(local_face)

  SELECT CASE(local)
  CASE(1:3)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_e1(elements(elem)%edges(faceedges(local))) ! E1 functions
  CASE (4:6)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_e2(elements(elem)%edges(faceedges(local-3))) ! E2 functions         
  CASE (7)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_f1(elements(elem)%faces(local_face)) ! F1 functions         
  CASE (8)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_f2(elements(elem)%faces(local_face)) ! F2 functions         
  CASE (9:11)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_e3(elements(elem)%edges(faceedges(local-8))) ! E3 functions         
  CASE (12)
    LOCAL_TO_GLOBAL_INDEX_PRE_TRI = renumbered_pre_f3(elements(elem)%faces(local_face)) ! F3 functions         
  CASE DEFAULT
    STOP 'IE: out-of-range local value in LOCAL_TO_GLOBAL_INDEX_PRE_TRI.'
  END SELECT

END FUNCTION LOCAL_TO_GLOBAL_INDEX_PRE_TRI
!*******************************************************************************


! Last changed 27 March 2003. Option to compute incident field for TD analysis
! also added
! Last changed: 10 Apr 2001. 
! New function CONVERT_2_LOCAL_PORT_COORD added. 
! Function K_Z_10 replaced by K_Z_MN
! Function E_10 extended; GAUSS_QUAD_VBF changed; GAUSS_QUAD_CFIELD ditto. 
! 20 March 2001. DBD. Function VECTOR_LENGTH added (lines 231-245). 
! 6 March 2001. DBD. Function VBF_Approx corrected.  
! 17 Feb 2001 DBD - new functions VBF,  VBF_S and NABLA_VBF added. 
! Functions FUNCTION GAUSS_QUAD_VBF and FUNCTION GAUSS_QUAD_CFIELD
! changed to use these. 
! Routines E_FIELD_CMPLX and E_FIELD_REAL were corrected. 


! Changed DBD 27 March 2003
FUNCTION GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,&
                        VBF_type,i1,i2,normal,i3,ell)
! End Changed DBD 27 March 2003
  
  USE basis_function, ONLY: VBF
  USE boundary_conditions
  USE gw_sys, ONLY : E_MN
  USE geometry
!  USE geometry, ONLY: XYZ_COORDINATES, GRADIENT_LAMBDA, FACE_AREA, &
!                      LOCAL_FACENODES
  USE gw_data
  USE math_tools, ONLY: CROSS_PRODUCT
  USE problem_info
  USE quad_tables
  USE TD_source
  USE scattering_analysis_data
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! Real-valued Gaussian quadrature over a triangle for
! the product of a vector basis function of type VBF_type and the incident field. 
! This may be the incident waveguide TE_mn mode, as required by 
! eqn. (8.87) Jin, (8.88 in Jin 2nd edn, 2002) "The Finite Element Method in Electromagnetics", Wiley 1993,
! extended to multi-mode analysis, 
! or the incident field in (12.10), Jin 2nd edn for total field analysis
! or as in the 4th term in (12.80) (Jin 2nd), for scattered field analysis.
! 
! 1, 3, 4 and 6 point rules are available. 
! References:
! Savage & Peterson, "Quadrature rules for numerical integration over 
! triangles and tetrahedra", IEEE AP Magazine, June 1996, pp. 100-102.
! See also Cowper, "Gaussian quadrature formulas for triangles", 
! Int. Jnl. Num. Meth. Eng., 1973, pp.405-408.
!
! May 2000: First version. DBD.
! 13 July 2000: Extended to support LT/QN elements. DBD.
! 16 Feb  2001: Changed to call the external function VBF to compute
!               the vector basis functions. 
! 28 March 2001: Extended to support arbitrarily-orientated ports. DBD.
! 27 March 2003: Extended to support other types of incident fields. 
!                Comment statements corrected (routine is valid fcr abritary
!                triangle orientation). 
! 21 May 2003:   Extended to support the scattered field formulation of Jin 2nd p.554.
! 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: port_or_ABC_num,element_num,& 
                              local_face_num,VBF_type
  INTEGER(I4B), INTENT(IN):: i1, i2            ! Local tet nodes making up 
                                               ! local triangle edges.
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal ! Outward unit normal on surface
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i3     ! As for i1 & i2, remaining node;
                                               ! only needed for face-based VBF's.
  REAL(SP), INTENT(IN), OPTIONAL :: ell        ! Edge length


! Added DBD 27 March 2003
  INTEGER(I4B) port_num,ABC_num                ! Number of port or ABC. 
! End added DBD 27 March 2003
  REAL(SP) :: GAUSS_QUAD_VBF
  REAL(SP), DIMENSION(3) :: Omega_12, Omega_13, Omega_23
  REAL(SP), DIMENSION(3) :: tempcoord
  INTEGER(I4B), DIMENSION(3) :: tri_nodes       ! Triangle nodes corresponding
                                               ! to tet nodes.
  INTEGER(I4B) n1,n2,n3                         ! Ditto.
  REAL(SP) :: x,y,z                            ! 3D coordinates. 
  REAL(SP), DIMENSION(3) :: e1_i1_i2           ! Vector basis function (VBF) e1.
  REAL(SP), DIMENSION(3) :: e2_i1_i2           ! Vector basis function (VBF) e2.
  REAL(SP), DIMENSION(3) :: f1                 ! Vector basis function (VBF) f1.
  REAL(SP), DIMENSION(3) :: f2                 ! Vector basis function (VBF) f2.
  REAL(SP), DIMENSION(3) :: n_X_VBF            ! n X  relevant VBF. 
  REAL(SP), DIMENSION(3) :: E_MN_X_nrml        ! Incident field mode cross n. 
  REAL(SP), DIMENSION(3) :: U_inc              ! Incident field at boundary. 
  REAL(SP), DIMENSION(3) :: U_inc_X_nrml       ! Incident field at boundary cross n. 
  REAL(SP), DIMENSION(3) :: temp_vec1,temp_vec2,temp_vec3 
                                               ! Temporary vectors. 
  INTEGER(I4B) :: ii,inode                     ! Misc. counter.
  REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
  REAL(SP), DIMENSION(4,3) :: nabla_lambda     ! Gradients of simplex
                                               ! coordinates. Constant
                                               ! within element.
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes  ! local nodes of the face
  REAL(SP) Y_char                              ! Characteristic impedance of homogeneous medium.


  ! Test for correct data:
  IF (SCALE_BY_EDGE_LENGTH.AND.ELEMENT_TYPE.EQ.1.AND..NOT.PRESENT(ell) & 
     .AND.(VBF_type.EQ.1.OR.VBF_type.EQ.2) ) THEN
    STOP 'IE: Edge length required for e1 or e2 scaling in GAUSS_QUAD_VBF with this basis function.'
  END IF 
  IF ((VBF_type.EQ.3.OR.VBF_type.EQ.4).AND..NOT.PRESENT(i3)) THEN
    STOP 'IE: Must supply third node for f1 or f2 in GAUSS_QUAD_VBF.'
  END IF
! Added DBD 27 March 2003

  IF(GW_ANALYSIS) THEN
    port_num = port_or_ABC_num
  ELSE IF (TD_ANALYSIS) THEN
    IF(.NOT.SCAT_FIELD) THEN 
      ABC_num = port_or_ABC_num
	ELSE 
	  Y_char = REAL(eps_r(HOMOG_MEDIUM)/mu_r(HOMOG_MEDIUM))/Z_zero ! Real-valued characteristic admittance of 
	                                         ! homogeneous background medium. 
	END IF
  END IF
! End added DBD 27 March 2003

  ! Get the (normalized) gradients of the simplex coordinates.
  nabla_lambda =  GRADIENT_LAMBDA(element_num,.true.)

  GAUSS_QUAD_VBF = 0.0 ! Initialize
  tempfacenodes = LOCAL_FACENODES(local_face_num)

  ! Perform quadrature:
  DO ii = 1,gauss_points

    ! Associate the correct 3D simplex coordinates with the 2D integration:
    lambda = 0.0
    DO inode = 1,3
      lambda(tempfacenodes(inode)) = quad_tri_rules(gauss_points)%rule(ii,inode)
    END DO

    ! Find values for x,y and z corresponding to sample points (all coordinates used). 
	tempcoord(1:3) = XYZ_COORDINATES(element_num,lambda)
    x = tempcoord(1)
    y = tempcoord(2)
    z = tempcoord(3)

    ! Find ^n cross VBF:
    IF (.NOT.PRESENT(i3)) THEN
      n_X_VBF =CROSS_PRODUCT(normal,VBF(VBF_type,lambda,nabla_lambda,i1,i2))
    ELSE 
      n_X_VBF =CROSS_PRODUCT(normal,VBF(VBF_type,lambda,nabla_lambda,i1,i2,i3))
    END IF

    ! Scale by edge lengths if appropriate:
    IF(ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH.AND. &
      (VBF_type.EQ.1.OR.VBF_type.EQ.2)) THEN
      n_X_VBF = n_X_VBF * ell
    END IF

! Changed DBD 27 March 2003
    IF(GW_ANALYSIS) THEN

      ! Find ^n cross WG mode:
      E_MN_X_nrml = CROSS_PRODUCT(E_MN(x,y,z,port_num,mode_m,mode_n),normal)

      ! Add weighted function value. Note - sign of dS carried in integrand.
      GAUSS_QUAD_VBF = GAUSS_QUAD_VBF +                          &
                       quad_tri_rules(gauss_points)%rule(ii,4) * & 
                       DOT_PRODUCT(n_X_VBF,E_MN_X_nrml)
    ELSE IF (TD_ANALYSIS) THEN
	  IF(.NOT.SCAT_FIELD) THEN
        temp_vec1 = CROSS_PRODUCT(normal,TD_D_BY_DT_INC_FIELD(x,y,z,'H'))
        temp_vec2 = CROSS_PRODUCT(normal,TD_D_BY_DT_INC_FIELD(x,y,z,'E'))
	    temp_vec3 = CROSS_PRODUCT(normal,temp_vec2)
        U_inc = ABCs(ABC_num)%Yc*(temp_vec1+temp_vec3)
        U_inc_X_nrml = CROSS_PRODUCT(U_inc,normal)
        ! Add weighted function value. Note - sign of dS carried in integrand.
        GAUSS_QUAD_VBF = GAUSS_QUAD_VBF +                          &
                         quad_tri_rules(gauss_points)%rule(ii,4) * & 
                         DOT_PRODUCT(n_X_VBF,U_inc_X_nrml)
      ELSE 
        GAUSS_QUAD_VBF = GAUSS_QUAD_VBF +                          &
                         quad_tri_rules(gauss_points)%rule(ii,4) * & 
                         DOT_PRODUCT(n_X_VBF,TD_D_BY_DT_INC_FIELD(x,y,z,'H',Y_char))
      END IF
    END IF				     
! End changed DBD 27 March 2003

  END DO ! Quadrature loop


  ! Scale by area of triangle (above results are normalized for unit area). 
  GAUSS_QUAD_VBF = GAUSS_QUAD_VBF * FACE_AREA(element_num,local_face_num)

END FUNCTION GAUSS_QUAD_VBF
!*******************************************************************************



FUNCTION GAUSS_QUAD_VBF_SPC(ABC_num,element_num,local_face_num,&
                        VBF_type,i1,i2,normal,i3,ell)
  USE basis_function, ONLY: VBF
  USE boundary_conditions
  USE geometry
  USE gw_data
  USE math_tools, ONLY: CROSS_PRODUCT
  USE problem_info
  USE quad_tables
  USE FD_scat_source
  IMPLICIT NONE
!*******************************************************************************
! Complex-valued Gaussian quadrature over a triangle for
! the product of a vector basis function of type VBF_type and the incident field. 
! Code based on GAUSS_QUAD_VBF, but does not support guided wave analysis.
! Presently, incident field as in (9.60, Jin 2nd edn) 
! "The Finite Element Method in Electromagnetics", Wiley 2002
! implemented - i.e. 1st order ABC only - for total field formulation. 
! 
!
! NB! ASSUMES EXTERNAL REGION IS FREE SPACE!!!
! 
! 1, 3, 4 and 6 point rules are available. 
! References:
! Savage & Peterson, "Quadrature rules for numerical integration over 
! triangles and tetrahedra", IEEE AP Magazine, June 1996, pp. 100-102.
! See also Cowper, "Gaussian quadrature formulas for triangles", 
! Int. Jnl. Num. Meth. Eng., 1973, pp.405-408.
!
! 16 Dec 2003 : First version. DBD.
! 26 Jul 2004 : Extended to support scattered field formulation (assuming PEC scatterer).
! 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: ABC_num,element_num,& 
                              local_face_num,VBF_type ! ABC_num not used at present.
  INTEGER(I4B), INTENT(IN):: i1, i2            ! Local tet nodes making up 
                                               ! local triangle edges.
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal ! Outward unit normal on surface
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i3     ! As for i1 & i2, remaining node;
                                               ! only needed for face-based VBF's.
  REAL(SP), INTENT(IN), OPTIONAL :: ell        ! Edge length


  COMPLEX(SPC) :: GAUSS_QUAD_VBF_SPC
  REAL(SP), DIMENSION(3) :: Omega_12, Omega_13, Omega_23
  REAL(SP), DIMENSION(3) :: tempcoord
  INTEGER(I4B), DIMENSION(3) :: tri_nodes       ! Triangle nodes corresponding
                                               ! to tet nodes.
  INTEGER(I4B) n1,n2,n3                         ! Ditto.
  REAL(SP) :: x,y,z                            ! 3D coordinates. 
  REAL(SP), DIMENSION(3) :: e1_i1_i2           ! Vector basis function (VBF) e1.
  REAL(SP), DIMENSION(3) :: e2_i1_i2           ! Vector basis function (VBF) e2.
  REAL(SP), DIMENSION(3) :: f1                 ! Vector basis function (VBF) f1.
  REAL(SP), DIMENSION(3) :: f2                 ! Vector basis function (VBF) f2.
  REAL(SP), DIMENSION(3) :: n_X_VBF            ! n X  relevant VBF. 
  REAL(SP) eta_0                               ! Wave impedance of free space
  COMPLEX(SPC), DIMENSION(3) :: U_inc              ! Incident field at boundary. 
  COMPLEX(SPC), DIMENSION(3) :: U_inc_X_nrml       ! Incident field at boundary cross n. 
  COMPLEX(SPC), DIMENSION(3) :: temp_vec1,temp_vec2,temp_vec3 
                                               ! Temporary vectors. 
  INTEGER(I4B) :: ii,inode                     ! Misc. counter.
  REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
  REAL(SP), DIMENSION(4,3) :: nabla_lambda     ! Gradients of simplex
                                               ! coordinates. Constant
                                               ! within element.
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes  ! local nodes of the face
!  REAL(SP) Y_char                              ! Characteristic impedance of homogeneous medium.


  ! Test for correct data:
  IF (SCALE_BY_EDGE_LENGTH.AND.ELEMENT_TYPE.EQ.1.AND..NOT.PRESENT(ell) & 
     .AND.(VBF_type.EQ.1.OR.VBF_type.EQ.2) ) THEN
    STOP 'IE: Edge length required for e1 or e2 scaling in GAUSS_QUAD_VBF with this basis function.'
  END IF 
  IF ((VBF_type.EQ.3.OR.VBF_type.EQ.4).AND..NOT.PRESENT(i3)) THEN
    STOP 'IE: Must supply third node for f1 or f2 in GAUSS_QUAD_VBF.'
  END IF
! Added DBD 27 March 2003

  IF(GW_ANALYSIS.OR.TD_ANALYSIS) THEN
    STOP 'Routine GAUSS_QUAD_VBF_SPC does not support GW_ANALYSIS or TD_ANALYSIS'
  END IF
  ! Get the (normalized) gradients of the simplex coordinates.
  nabla_lambda =  GRADIENT_LAMBDA(element_num,.true.)

  GAUSS_QUAD_VBF_SPC = 0.0 ! Initialize
  tempfacenodes = LOCAL_FACENODES(local_face_num)

  ! Perform quadrature:
  DO ii = 1,gauss_points

    ! Associate the correct 3D simplex coordinates with the 2D integration:
    lambda = 0.0
    DO inode = 1,3
      lambda(tempfacenodes(inode)) = quad_tri_rules(gauss_points)%rule(ii,inode)
    END DO

    ! Find values for x,y and z corresponding to sample points (all coordinates used). 
	tempcoord(1:3) = XYZ_COORDINATES(element_num,lambda)
    x = tempcoord(1)
    y = tempcoord(2)
    z = tempcoord(3)

    ! Find ^n cross VBF:
    IF (.NOT.PRESENT(i3)) THEN
      n_X_VBF =CROSS_PRODUCT(normal,VBF(VBF_type,lambda,nabla_lambda,i1,i2))
    ELSE 
      n_X_VBF =CROSS_PRODUCT(normal,VBF(VBF_type,lambda,nabla_lambda,i1,i2,i3))
    END IF

    ! Scale by edge lengths if appropriate:
    IF(ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH.AND. &
      (VBF_type.EQ.1.OR.VBF_type.EQ.2)) THEN
      n_X_VBF = n_X_VBF * ell
    END IF
	IF(.NOT.SCAT_FIELD) THEN
      ! Note - routine FD_INC_FIELD called to return H without Y scaling, hence eta in temp_vec1
      ! is effectively cancelled out.
      temp_vec1 = -CROSS_PRODUCT(CMPLX(normal,0.0_SP),FD_INC_FIELD(x,y,z,'H'))
      temp_vec2 = CROSS_PRODUCT(CMPLX(normal,0.0_SP),FD_INC_FIELD(x,y,z,'E'))
	  temp_vec3 = CROSS_PRODUCT(CMPLX(normal,0.0_SP),temp_vec2)
      U_inc = (temp_vec1+temp_vec3)
      U_inc_X_nrml = CROSS_PRODUCT(U_inc,CMPLX(normal,0.0_SP))
      ! Add weighted function value. Note - sign of dS carried in integrand.
      GAUSS_QUAD_VBF_SPC = GAUSS_QUAD_VBF_SPC +                          &
                      quad_tri_rules(gauss_points)%rule(ii,4) * & 
                      DOT_PRODUCT(n_X_VBF,U_inc_X_nrml)
    ELSE 
	  temp_vec1 = FD_INC_FIELD(x,y,z,'H',(1/Z_zero))
      temp_vec2 = CROSS_PRODUCT(CMPLX(normal,0.0_SP),FD_INC_FIELD(x,y,z,'E'))
      ! Add weighted function value. (Scaling done in calling routine)
      GAUSS_QUAD_VBF_SPC = GAUSS_QUAD_VBF_SPC +                          &
                      quad_tri_rules(gauss_points)%rule(ii,4) * & 
                      DOT_PRODUCT(CONJG(temp_vec1),temp_vec2)
					  ! DOT_PRODUCT conjugates first argument - undesired here, so conjugate first to correct!
	END IF
  END DO ! Quadrature loop

  ! Scale by area of triangle (above results are normalized for unit area). 
  GAUSS_QUAD_VBF_SPC = GAUSS_QUAD_VBF_SPC * FACE_AREA(element_num,local_face_num)

END FUNCTION GAUSS_QUAD_VBF_SPC
!*******************************************************************************




FUNCTION GAUSS_QUAD_CFIELD(port_num,element_num,local_face_num,&
                           e1,e2,e3,f1,f2,f3,ell)
  USE boundary_conditions
  USE gw_sys, ONLY : E_MN
  USE geometry, ONLY: XYZ_COORDINATES, GRADIENT_LAMBDA, FACE_AREA, &
                      LOCAL_FACENODES
  USE gw_data
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! Gaussian quadrature over a triangle in a plane of constant z for
! the product of the interpolated Complex interpolated FIELD using these 
! basis funtions and  the incident waveguide TE_mn mode, as required by 
! eqns. (8.91) and (8.92), Jin, "The Finite Element Method in Electromagnetics",
! Wiley 1993 and extended for multi-mode cylindrical (rectangular) waveguide. 
!
! References:
! Savage & Peterson, "Quadrature rules for numerical integration over 
! triangles and tetrahedra", IEEE AP Magazine, June 1996, pp. 100-102.
! See also Cowper, "Gaussian quadrature formulas for triangles", 
! Int. Jnl. Num. Meth. Eng., 1973, pp.405-408.
!
! May 2000: First version. DBD.
! 13 July 2000: Extended to also support LT/QN elements. DBD.
! 28 March 2001: Extended to support arbitrarily-orientated ports. DBD.
!  2 March 2002: Extended to also support LT/LN, QT/QN elements. DBD.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: port_num, element_num, & 
                              local_face_num
  COMPLEX(SPC), DIMENSION(3), INTENT(IN) :: e1
  COMPLEX(SPC), DIMENSION(3), INTENT(IN), OPTIONAL :: e2,e3
  COMPLEX(SPC), INTENT(IN), OPTIONAL :: f1,f2,f3
  REAL(SP), DIMENSION(3), INTENT(IN),OPTIONAL  :: ell
  COMPLEX(SPC) :: GAUSS_QUAD_CFIELD

  REAL(SP), DIMENSION(3) :: tempcoord
  REAL(SP) :: x,y,z                            ! 3D cordinates. 
  COMPLEX(SPC), DIMENSION(3) :: E_field        ! E field on element.
  INTEGER(I4B) :: ii,inode                     ! Misc. counter.
  REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
  REAL(SP), DIMENSION(4,3) :: nabla_lambda     ! Gradients of simplex
                                               ! coordinates. Constant
                                               ! within element.
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes  ! local nodes of the face

  IF (SCALE_BY_EDGE_LENGTH.AND..NOT.PRESENT(ell)) THEN
    STOP 'IE: Edge length required if scaling requested in GAUSS_QUAD_CFIELD.' 
  END IF 
 
  GAUSS_QUAD_CFIELD = ZERO_C ! Initialize
  tempfacenodes = LOCAL_FACENODES(local_face_num)

  ! Get the (normalized) gradients of the simplex coordinates.
  nabla_lambda =  GRADIENT_LAMBDA(element_num,.true.)

  ! Perform quadrature:
  DO ii = 1,gauss_points
    ! Associate the correct 3D simplex coordinates with the 2D integration:
    lambda = 0.0
    DO inode = 1,3
      lambda(tempfacenodes(inode)) = quad_tri_rules(gauss_points)%rule(ii,inode)
    END DO

    ! Find values for x corresponding to sample points (y and z also
    ! found but not used).
    tempcoord(1:3) = XYZ_COORDINATES(element_num,lambda)
    x = tempcoord(1)
    y = tempcoord(2)
    z = tempcoord(3)

    ! Calculate the value <E_field>:
    CALL VBFApprox

    ! Add weighted function value. 
    GAUSS_QUAD_CFIELD = GAUSS_QUAD_CFIELD +                           &
      quad_tri_rules(gauss_points)%rule(ii,4) *                       &
      DOT_PRODUCT(CONJG(E_field),E_MN(x,y,z,port_num,mode_m,mode_n))
    ! Note - DOT_PRODUCT conjugates the first argument, which is undesired, 
    ! hence the further  conjugate to correct this.
  END DO

  ! Scale by area of triangle (above results are normalized for unit area). 
  GAUSS_QUAD_CFIELD = GAUSS_QUAD_CFIELD * FACE_AREA(element_num,local_face_num)

CONTAINS 

SUBROUTINE VBFApprox
  USE basis_function, ONLY : VBF
  USE geometry, ONLY : elements, LOCAL_TRI_EDGENODES, LOCAL_FACENODES, & 
                                 MAX_ORDER, MIXED_ORDER
                          
  IMPLICIT NONE
!*******************************************************************************
! This internal subroutine returns the approximated field on the triangular
! face for hierarchal elements at the required quaudature point.
! (Note that this is the FIELD, not the normal cross-producted with the field, 
! thus it is the sum of the appropriate (3D)
! vector basis functions, not the surface vector basis functions. 
! May 2000: First version. DBD.
! 13 July 2000: Extended to support LT/QN elements. DBD.
! Re-written DBD 17 Feb 2001, calling external function VBF (thus removing
! the explicit basis function definitions from this routine) and using
! function LOCAL_TRI_EDGENODES to find the local triangular edge nodes. 
! Corrected DBD 6 March 2001: e2,f1 and f2 functions no longer evaluated
! for CT/LN elements (could have caused unpredictable results). 
! Extended DBD 2 Mar 2002: support added for up to QT/QN elements. 
!*******************************************************************************
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE,3) :: func_values      ! (Unscaled) basis function values
  INTEGER(I4B) :: count                        ! Counter
  INTEGER(I4B), DIMENSION(2) :: tempnodes
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes

  IF (ELEMENT_TYPE.NE.1.AND.SCALE_BY_EDGE_LENGTH) THEN
    STOP 'IE in VBFApprox. Scaling by edge length only supported for S&P elements'
  END IF 

  ! Initialize func_values to zero. 
  func_values = 0.0_SP

  ! Evaluate the three E1 functions:
  DO count = 1,3
    tempnodes = LOCAL_TRI_EDGENODES(count,local_face_num)
    func_values(count,1:3) = VBF(1,lambda,nabla_lambda,tempnodes(1),tempnodes(2))
    IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
      func_values(count,1:3) = ell(count)*func_values(count,1:3)
    END IF
  END DO
  
  ! LT/LN and higher order terms.
  IF (MAX_ORDER(element_num,local_face_num).GE.2.OR.&
      .NOT.MIXED_ORDER(element_num,local_face_num).AND.&
      MAX_ORDER(element_num,local_face_num).GE.1) THEN
    ! Evaluate the three E2 functions:
    DO count = 1,3
      tempnodes = LOCAL_TRI_EDGENODES(count,local_face_num)
      func_values(count+3,1:3) = VBF(2,lambda,nabla_lambda,tempnodes(1),tempnodes(2))
      IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
        func_values(count+3,1:3) = ell(count)*func_values(count+3,1:3)
      END IF
    END DO
  END IF

  ! LT/QN and higher order terms.
  IF (MAX_ORDER(element_num,local_face_num).GE.2) THEN
    ! Evaluate the single F1 function:
    tempfacenodes = LOCAL_FACENODES(local_face_num)
    func_values(7,1:3) = VBF(3,lambda,nabla_lambda,& 
                            tempfacenodes(1),tempfacenodes(2),tempfacenodes(3))

    ! Evaluate the single F2 function:
    tempfacenodes = LOCAL_FACENODES(local_face_num)
    func_values(8,1:3) = VBF(4,lambda,nabla_lambda,& 
                            tempfacenodes(1),tempfacenodes(2),tempfacenodes(3))

  END IF

  ! QT/QN and higher order terms.
  IF (MAX_ORDER(element_num,local_face_num).GE.3.OR.&
      .NOT.MIXED_ORDER(element_num,local_face_num).AND. &
	  MAX_ORDER(element_num,local_face_num).GE.2) THEN
    ! Evaluate the three E3 functions:
    DO count = 1,3
      tempnodes = LOCAL_TRI_EDGENODES(count,local_face_num)
      func_values(count+8,1:3) = VBF(5,lambda,nabla_lambda,tempnodes(1),tempnodes(2))
    END DO
    ! Evaluate the single F3 function:
    tempfacenodes = LOCAL_FACENODES(local_face_num)
    func_values(12,1:3) = VBF(6,lambda,nabla_lambda,& 
                            tempfacenodes(1),tempfacenodes(2),tempfacenodes(3))
  END IF

  IF (MAX_ORDER(element_num,local_face_num).GE.4.OR.&
      .NOT.MIXED_ORDER(element_num,local_face_num).AND.&
	 MAX_ORDER(element_num,local_face_num).GE.3) THEN
	STOP 'IE: Unimplemted hierarchal order in VBFApprox.'
  END IF

  ! Compute E field. Note that func_values was initialized to zero.
  ! Similarly, the coefficients e1-3 and f1-3 were initialized 
  ! as zero.
  ! Thus it is not necessary to explicitly differentiate between 
  ! various hierarchal orders here. 
  
  E_field(1:3) = e1(1)*func_values(1,1:3) + e1(2)*func_values(2,1:3) + &
                 e1(3)*func_values(3,1:3) + &
                 e2(1)*func_values(4,1:3) + e2(2)*func_values(5,1:3) + &
                 e2(3)*func_values(6,1:3) + &
                 f1*func_values(7,1:3) + & 
                 f2*func_values(8,1:3) + &
				 e3(1)*func_values(9,1:3) + e3(2)*func_values(10,1:3) + &
                 e3(3)*func_values(11,1:3) + & 
				 f3*func_values(12,1:3)
END SUBROUTINE VBFApprox

END FUNCTION GAUSS_QUAD_CFIELD
!******************************************************************************


SUBROUTINE MESH_INFO_WRITE 
  USE geometry
  USE matrix
  USE problem_info
  USE quad_tables
  USE unit_numbers
  use scattering_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! Write out some mesh information.
!*******************************************************************************
  INTEGER(I4B) iedge                              ! counters
  INTEGER(I4B), DIMENSION(2) :: tempnodes         ! temp storage


  WRITE(FILEOUT,'(//,20X,A,/)') 'MESH INFORMATION'
  WRITE(FILEOUT,'(1X,A,I8)')    'Number of nodes:       ',num_nodes
  WRITE(FILEOUT,'(1X,A,I8)')    'Number of elements:    ',num_elements
  WRITE(FILEOUT,'(1X,A,I8)')    'Number of edges:       ',num_edges
  WRITE(FILEOUT,'(1X,A,I8)')    'Number of faces:       ',num_faces
  ! Find and print average edge length in mesh:
  avg_edge_length  = 0.0_SP
  DO iedge = 1,num_edges
    tempnodes(:) = edges(iedge)%nodes
    avg_edge_length = avg_edge_length + T_LENGTH(tempnodes(1),tempnodes(2))
  END DO
  avg_edge_length = avg_edge_length/num_edges
  WRITE(FILEOUT,'(1X,A,G10.4)') 'Average edge length:   ',avg_edge_length


  WRITE(FILEOUT,'(//,20X,A,/)') 'MATRIX INFORMATION'
  WRITE(FILEOUT,'(A,I8)')  'No. of degrees of freedom (free edges & faces):  ', dof
  IF (FD_SCAT_ANALYSIS.AND.SCAT_FIELD) THEN
    WRITE(FILEOUT,'(A,I8)')  'No. of non-zero prescribed "degrees of freedom" (prescribed edges & faces):  ', pre_dof
  END IF

  SELECT CASE (ELEMENT_TYPE)
  CASE(1)
     WRITE(FILEOUT,'(A)')   'Formulation used for S&T:                        Savage and Peterson'
  CASE(2)
     WRITE(FILEOUT,'(A)')   'Formulation used for S&T:                        Andersen and Volakis'
  CASE(3)
     WRITE(FILEOUT,'(A)')   'Formulation used for S&T:                        Webb99'
  CASE DEFAULT
     STOP 'IE: Error in MESH_INFO_WRITE. Unimplemented element type' 
  END SELECT
  IF(SCALE_BY_EDGE_LENGTH) THEN
    WRITE(FILEOUT,'(A)')   'S and T matrices scaled by edge lengths:         Yes'
  ELSE
    WRITE(FILEOUT,'(A)')   'S and T matrices scaled by edge lengths:         No'
  END IF

!    WRITE(FILEOUT,'(I3,A)') gauss_points,&
!      '-point rule used for quadrature on simplexes'


END SUBROUTINE MESH_INFO_WRITE
!*******************************************************************************


SUBROUTINE DIRECT_SOLVE(flag,time_taken)
  USE math_tools, ONLY: TIME_DIFFERENCE
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! This routine has multiple functions, depending on <flag>:
! flag=0: LU decomposition of the matrix <A_mat_c>.
! flag=1: Solving <x_vec_c> by back substitution. (A_mat_c*x_vec_c = b_vec_c)
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: flag        ! 0=LU decomp.,1=back subs.
  REAL(SP) :: time_taken     ! time taken by this routine

!  REAL(SP) :: time_taken     ! time taken by the routine
  INTEGER(I4B) :: info                                 ! Info from LAPACK factorization.
  INTEGER(I4B), DIMENSION(8) :: time_start,time_finish ! temp variables for timing all events

  CALL DATE_AND_TIME (values=time_start)

  SELECT CASE (flag)

  CASE (0) ! LU decomposition of A_mat_c:
    IF (ALLOCATED(ipiv)) DEALLOCATE(ipiv,STAT=info)
    ALLOCATE(ipiv(dof))
    ! Uses LAPACK routine SGETRF to LU factorize (single precision)
    CALL CGETRF (dof,dof,A_mat_c,dof,ipiv,info)  
    IF(Info.NE.0) CALL ERROR_FEMFEKO(1,4302,int1=info) ! Check error condition

  CASE (1) ! Solve x_vec_c by back substitution:
    CALL CGETRS('N',dof,1,A_mat_c,dof,ipiv,b_vec_c,dof,info) ! LAPACK routine
    IF(Info.NE.0) CALL ERROR_FEMFEKO(1,4303,int1=info) ! Check error condition
    x_vec_c = b_vec_c

  END SELECT

  CALL DATE_AND_TIME (values=time_finish)
  time_taken = TIME_DIFFERENCE(time_start,time_finish)

  ! Write to the output file:
  SELECT CASE (flag)
  CASE (0) ! LU decomposition of A_mat_c:
    WRITE (FILEOUT,'(//,20X,A)')   'DIRECT SOLUTION OF THE SYSTEM MATRIX EQUATION: LU DECOMPOSITION'
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  LU decomposition'
    WRITE (FILEOUT,'(1X,A)')       'Subroutine(s) used:              CGETRF (LAPACK)'
    WRITE (FILEOUT,'(1X,A,F12.3)') 'Time taken by LU decomp.(sec):  ', time_taken
  CASE (1) ! Solve x_vec_c by back substitution:
    WRITE (FILEOUT,'(//,20X,A)')   'DIRECT SOLUTION OF THE SYSTEM MATRIX EQUATION: BACK SUBSTITUTION'
    WRITE (FILEOUT,'(/,1X,A)')     'Technique used:                  LU back substitution'
    WRITE (FILEOUT,'(1X,A)')       'Subroutine(s) used:              CGETRS (LAPACK)'
    WRITE (FILEOUT,'(1X,A,F12.3)') 'Time taken by back subs.(sec):  ', time_taken
  END SELECT

END SUBROUTINE DIRECT_SOLVE
!*******************************************************************************


SUBROUTINE NUMBER_DOF
  USE eigen_analysis_data
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
!*******************************************************************************
!*** Establish and (re)-number degrees of freedom                           ****
!*******************************************************************************
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! DB Davidson
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 20 March 2000 by DBD.
! 21 March 00: initialization of renumbered_e1,e2,f1 and f2 changed to integer 0 
!              (not wrong, but confusing). DBD.
! 13 May 00  : DEBUG_DOF moved to end of this routine. DBD.
! 16 Nov 2000: MMB: Added the FIND_FREEDOM routine and its dependants for a
!              unified identification process of edges and faces (identified as
!              free/CBAA_aperture etc.)
! 14 Dec 2000: MMB: Added code in subroutine FACE_NUMBERING_LIST that ensures that
!              in the CBAA_ANALYSIS case, the aperture faces are numbered first,
!              together with the aperture edges.
! 8 June 2001: DBD. Subroutine DEBUG_DOF extended to report PEC boundary
!              condition
!              Subroutine PEC_QUADRILATERAL_SEARCH corrected to 
!              also flag faces (around line 1010)
! 28 Feb 2002: DBD. Support added for polynomial complete elements LT/LN and
!              QT/QN  elements.
! 17 Mar 2002: MMB. Corrected the dof assignment in the case of the f3 
!              dofs. They were (incorrectly) only assigned in the case of full order
!              elements of order >= 2.
!*******************************************************************************
  IMPLICIT NONE
  INTEGER(I4B) i_loop,pos
  INTEGER(I4B) iedge,ielem,jface,kedge,kface,jdof             ! Counters.
  INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE :: edge_facelink
 
  CALL DOF_CALC ! internal routine to verify the number of dof's

  ! Check that there are a non-zero number of degrees of freedom
  IF (dof.LE.0) THEN
    CALL ERROR_FEMFEKO(1,4039)
  END IF

  CALL FACE_NUMBERING_LIST ! Make the list edge_facelink for subsequent use.

  ! Number the free edges and faces for later use in assembling S and T.
  ! Note that the d.o.f. number is unique - a renumbered d.o.f. is 
  ! either an edge (type e1 or e2) or face (f1 or f2).
  ! For the CT/LN case, the numbers in renumbered_e1 will be
  ! consecutive.
  ! For the LT/QN renumbered d.o.f.'s, 
  ! the numbering for each type will have "missing" numbers; 
  ! these have of course been assigned to the other d.o.f.'s. 

  ! Allocate and initialize degree of freedom variables
  ! Note that the allocation is conservative: num_edges and 
  ! num_faces must be >= the respective d.o.f.'s
  ALLOCATE (renumbered_e1(num_edges)) ! CT/LN  
  ALLOCATE (renumbered_e2(num_edges)) ! LT/LN and LT/QN
  ALLOCATE (renumbered_e3(num_edges)) ! QT/QN
  ALLOCATE (renumbered_f1(num_faces)) ! LT/QN
  ALLOCATE (renumbered_f2(num_faces)) ! LT/QN
  ALLOCATE (renumbered_f3(num_faces)) ! QT/QN
  ALLOCATE (dof_type(3*num_edges+3*num_faces))    
  renumbered_e1 = 0 ! Array initialization (0 implies prescribed or non-existent)
  renumbered_e2 = 0 ! Ditto
  renumbered_e3 = 0 ! Ditto
  renumbered_f1 = 0 ! Ditto
  renumbered_f2 = 0 ! Ditto
  renumbered_f3 = 0 ! Ditto
  dof_type      = 0 ! Ditto

  jdof  = 0 ! initialise

  ! Re-written to handle mixed/complete elements of varying orders, DBD. 
  DOF_ASSIGN_EDGE_LOOP: DO iedge = 1,num_edges
    IF (edges(iedge)%free) THEN ! D.o.f. associated with free edge
      ! Always do the following: CT/LN 
	  jdof  = jdof  + 1 
      renumbered_e1(iedge) = jdof
      dof_type(jdof)  = 1
      IF(edges(iedge)%order.GE.2.OR.&
	  	.NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.1) THEN
        ! LT/LN and LT/QN
	    jdof  = jdof  + 1
        renumbered_e2(iedge) = jdof  
        dof_type(jdof)  = 2                
      END IF
      IF(edges(iedge)%order.GE.3.OR.&
	    .NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.2) THEN
        ! QT/QN and QT/CuN
        jdof  = jdof  + 1
        renumbered_e3(iedge) = jdof  
        dof_type(jdof)  = 5
	  END IF                
    END IF

    ! Now assign the dofs assiciated with the faces that are associated
    ! with the current edge (i.e. iegde):
    CALL count_nonzeros(SIZE(edge_facelink(iedge,:)),edge_facelink(iedge,:),pos)  
    DOF_ASSIGN_FACE_LOOP: DO i_loop = 1,pos  

      kface = edge_facelink(iedge,i_loop)
      IF (kface.NE.0) THEN ! i.e. free face and this edge indeed 
                           ! linked to a face (not all edges are
                           ! linked to faces for this renumbering scheme)   

        ! Note: CT/LN and LT/LN elements have no face-based dof's
        IF (faces(kface)%order.GE.2) THEN
          jdof  = jdof  + 1 
          renumbered_f1(kface) = jdof  
          dof_type(jdof) = 3                              
          jdof  = jdof  + 1 
          renumbered_f2(kface) = jdof
          dof_type(jdof) = 4
        END IF
        IF (((.NOT.faces(kface)%mixed).AND.(faces(kface)%order.EQ.2)) &
		  .OR.((faces(kface)%order.GE.3))) THEN
		  jdof  = jdof  + 1 
          renumbered_f3(kface) = jdof
          dof_type(jdof) = 6
        END IF
      END IF
    END DO DOF_ASSIGN_FACE_LOOP
  END DO DOF_ASSIGN_EDGE_LOOP

  ! Consistency check following renumbering.
  IF(jdof.NE.dof) THEN
    CALL ERROR_FEMFEKO(1,4040)
  END IF

  IF (DEBUG_SYSTEM_MATRIX) THEN
    CALL DEBUG_DOF
  END IF



  DEALLOCATE(edge_facelink)

CONTAINS

!*******************************************************************************

  SUBROUTINE DEBUG_DOF
    IMPLICIT NONE
!*******************************************************************************
! List degrees of freedom - for debugging only
! Extended DBD 8 July 2001
! Not extended to include elements beyond LT/QN. 
!*******************************************************************************
    INTEGER(I4B) iedge,iface,jdof,kedge,kface

    IF (HIERARCHAL_ORDER.GE.2.AND..NOT.MIXED_ORDER_FLAG) THEN
	  STOP 'IE: DEBUG_PRE_DOF: Only implemented to LT/QN'
	END IF 

    WRITE(FILEOUT,'(//A)') '****** DEGREES OF FREEDOM **************'
    WRITE(FILEOUT,'(/A)') 'Requested by flag DEBUG_SYSTEM_MATRIX'

    IF(GW_ANALYSIS) THEN
      WRITE(FILEOUT,'(/A/)') 'Free/prescribed/port edges and faces:'
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A,2X,A,1X,A,2X,A/)') & 
        'Edge','Node 1', 'Node 2','Free (T or F)', 'Port (T or F)',&
         'Port #', 'PEC (T or F)'
      DO kedge = 1,num_edges
        WRITE(FILEOUT,'(I8,3X,2(I8,2X),11X,L1,11X,L1,11X,I2,5X,L1)') & 
          kedge, edges(kedge)%nodes(:), edges(kedge)%free, edges(kedge)%port, &
          edges(kedge)%portnumber, edges(kedge)%PEC
      END DO 
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A,4X,A,2X,A,1X,A,2X,A/)') & 
         'Face','Node 1', 'Node 2','Node 3','Free (T or F)', 'Port (T or F)',&
         'Port #', 'PEC (T or F)'
      DO kface = 1,num_faces
        WRITE(FILEOUT,'(I8,3X,3(I8,2X),11X,L1,11X,L1,11X,I2,5X,L1)') & 
          kface, faces(kface)%nodes, faces(kface)%free, faces(kface)%port, &
          faces(kface)%portnumber,faces(kface)%PEC
      END DO 
    ELSE IF (FD_SCAT_ANALYSIS) THEN 
	  WRITE(FILEOUT,'(/A/)') 'Free/prescribed edges and faces:'
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A,4X,A,4X,A/)') 'Edge','Node 1', & 
                                      'Node 2','Free (True or False)','Scat/tot boundary (T/F)','Dirichlet'
      DO kedge = 1,num_edges
        WRITE(FILEOUT,'(I8,3X,2(I8,2X),9X,L1,19X,L1,24X,L1)') kedge, edges(kedge)%nodes(:), &
                                             edges(kedge)%free, edges(kedge)%scat_tot_boundary,edges(kedge)%Dirichlet
      END DO 
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A,4X,A,4X,A,4X,A/)') 'Face','Node 1', & 
                                      'Node 2','Node 3','Free (True or False)','Scat/tot boundary (T/F)','Dirichlet'
      DO kface = 1,num_faces
        WRITE(FILEOUT,'(I8,3X,3(I8,2X),9X,L1,19X,L1,24X,L1)') kface, faces(kface)%nodes, &
                                            faces(kface)%free, faces(kface)%scat_tot_boundary, faces(kface)%Dirichlet
      END DO 
    ELSE
      WRITE(FILEOUT,'(/A/)') 'Free/prescribed edges and faces:'
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A/)') 'Edge','Node 1', & 
                                      'Node 2','Free (True or False)'
      DO kedge = 1,num_edges
        WRITE(FILEOUT,'(I8,3X,2(I8,2X),9X,L1)') kedge, edges(kedge)%nodes(:), &
                                             edges(kedge)%free
      END DO 
      WRITE(FILEOUT,'(/5X,A,7X,A,4X,A,4X,A,4X,A/)') 'Face','Node 1', & 
                                      'Node 2','Node 3','Free (True or False)'
      DO kface = 1,num_faces
        WRITE(FILEOUT,'(I8,3X,3(I8,2X),9X,L1)') kface, faces(kface)%nodes, &
                                            faces(kface)%free
      END DO 
    END IF
    WRITE(FILEOUT,'(/A)') 'Degrees of freedom by type:'
    WRITE(FILEOUT,'(/A,7X,A,4X,A/)') 'Type','D.o.f.','Old global edge/face#'               

    DO iedge = 1,num_edges      
      IF (renumbered_e1(iedge).NE.0) THEN ! Always
        WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'e1',renumbered_e1(iedge),iedge
      END IF
      IF(edges(iedge)%order.GE.2.OR.&
	  	.NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.1) THEN
        IF (renumbered_e2(iedge).NE.0) THEN ! LT/LN and higher 
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'e2',renumbered_e2(iedge),iedge
		END IF
      END IF
	  ! Not implemented yet for higher order elements. 	 
    END DO 

    DO iface = 1,num_faces
      SELECT CASE(HIERARCHAL_ORDER)
      CASE(1) ! CT/LN     
        ! No action
      CASE(2) ! LT/QN
        IF (renumbered_f1(iface).NE.0) THEN
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'f1',renumbered_f1(iface),iface
        END IF
        IF (renumbered_f2(iface).NE.0) THEN
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'f2',renumbered_f2(iface),iface
        END IF
      CASE DEFAULT
        STOP 'Error in routine DEBUG_DOF - unimplemented hierarchal order'
      END SELECT
    END DO

    WRITE(FILEOUT,'(//A)') 'Degrees of freedom in sequence'
    WRITE(FILEOUT,'(/A,7X,A,4X,A/)') 'Type','D.o.f.'    

    DO jdof = 1,dof 
      SELECT CASE(dof_type(jdof))
      CASE(1) ! e1     
        WRITE(FILEOUT,'(A,3X,I8)') 'e1',jdof              
      CASE(2) ! e2     
        WRITE(FILEOUT,'(A,3X,I8)') 'e2',jdof  
      CASE(3) ! e3     
        WRITE(FILEOUT,'(A,3X,I8)') 'f1',jdof  
      CASE(4) ! e4     
        WRITE(FILEOUT,'(A,3X,I8)') 'f2',jdof
      CASE DEFAULT
        STOP 'IE: Internal error in routine EIGEN_SYSMAT - incorrect d.o.f. type'
      END SELECT
    END DO         

  END SUBROUTINE DEBUG_DOF
!*******************************************************************************

  SUBROUTINE DOF_CALC
    USE nrtype
    IMPLICIT NONE
!*******************************************************************************
! Calculates the number of d.o.f.'s for verification purposes.
! This routine assumes MAX(HIERARCHAL_ORDER) = 2
! and supports both mixed order and polynomial complete elements.
! It does not assume that ALL elements in the mesh have the same order and 
! mixed/complete nature.
!*******************************************************************************
    INTEGER(I4B) :: ielem,iedge,iface
    INTEGER(I4B), DIMENSION(2) :: temp_edge_nodes

    ! Initial (max) size of dof: (QN/QN elements)
      dof = 3*num_edges + 3*num_faces

    ! Subtract the effect of edges that are not free or not 
    ! of the highest order:
    DO iedge = 1,num_edges
      IF (.NOT.edges(iedge)%free) THEN
        dof = dof - 3 
      ELSE
        SELECT CASE(edges(iedge)%order)
        CASE(1) 
  		  IF(edges(iedge)%mixed) THEN ! CT/LN
            dof = dof - 2
		  ELSE ! LT/LN
            dof = dof - 1
		  END IF
        CASE(2) 
		  IF(edges(iedge)%mixed) THEN ! LT/QN
            dof = dof - 1
		  ELSE
            CONTINUE ! QT/QN
		  END IF
        CASE DEFAULT
          STOP 'IE: Invalid hierarchal order in routine DOF_CALC.'
        END SELECT      
      END IF       
    END DO

    ! Subtract the effect of faces that are not free or not
    ! of the highest order:
    DO iface = 1,num_faces
      IF (.NOT.faces(iface)%free) THEN
        dof = dof - 3
      ELSE
        SELECT CASE(faces(iface)%order)
        CASE(1) ! CT/LN and LT/LN, no face-based dof's          
          dof = dof - 3 ! 3 dof's associated with a free face
		                ! for QN/QN elements..
        CASE(2) ! LT/QN
		  IF (faces(iface)%mixed) THEN ! LT/QN: remove one QN/QN
		                               ! dof.
            dof = dof - 1
		  ELSE 
            CONTINUE ! QT/QN
		  END IF
        CASE DEFAULT
          STOP 'IE: Invalid hierarchal order in routine DOF_CALC.'
        END SELECT      
      END IF       
    END DO

  END SUBROUTINE DOF_CALC
!*******************************************************************************

  SUBROUTINE FACE_NUMBERING_LIST
    USE problem_info
    IMPLICIT NONE
!*******************************************************************************
! Set up an index for the d.o.f.'s associate with the faces. 
! With CT/LN elements, the d.o.f.'s are associated only with the edges.
! The d.o.f. can be numbered by stepping through the edges in the mesh, element
! by element. With LT/QN elements, there are additional d.o.f.'s associated
! with edges AND faces. Although the face-based d.o.f.'s could also be
! numbered in a subsequent loop over the elements, face by face, 
! after the edge-based d.o.f.'s have been numbered, this will result in a 
! very large matrix bandwidth. 
! 
! Instead, each face is linked to an edge (which forms part of
! the face.)  The edge is not required to be free, and several faces 
! may typically be linked to a single edge.
! The face then assumes a d.o.f. number relatively close to that of the edge.   
!
! Two lists edge_facelink  and face_edgelink are created for this procedure. 

! The list edge_facelink assigns a list of faces to a single edge.  
! Row one, for example, lists faces linked to edge one. 
! Note that this is NOT a complete list of faces sharing edge one; 
! it is however sufficient for the renumbering process. 
! Note further that some edges might not have any
! faces linked to them at all. ONLY free faces are linked to edges 
! (which do not necessarily have to be free)
! 
! The list face_edgelink contains essentially the same information but 
! listed per face; all FREE faces will be linked to one and only one edge
! in this scheme, and otherwise will have a (default) zero entry.
!
! 14 Dec 2000: MMB: Added <secondary_lookup_table> in order that, in case of a
!                   LT/QN, CBAA_analysis NO FACE NOT BELONGING TO THE CBAA
!                   APERTURE WILL BE ASOCIATED WITH A CBAA APERTURE EDGE. This
!                   is done in order that the CBAA aperture dof's are assigned
!                   first, so that they lead to a fully populated upper left 
!                   portion of the system matrix.
!*******************************************************************************
    INTEGER(I4B), DIMENSION(4) :: lookup_table
    INTEGER(I4B), DIMENSION(4) :: secondary_lookup_table
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: face_edgelink
    INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE :: temp_edge_facelink
    INTEGER(I4B) tempedge,tempface                 
    INTEGER(I4B) EST_DIM,TMP_DIM ! dimensioning info

    ! Initialize some variables used for the LT/QN elements, used in 
    ! internal routine FACE_NUMBERING_LIST.
    EST_DIM = 2
    ALLOCATE(edge_facelink(num_edges,EST_DIM)) 
                                        ! EST_DIM is an estimate for the 
                                        ! maximum number of faces connected
                                        ! to a single edge; it is increased
                                        ! if needed. 
    ALLOCATE(face_edgelink(num_faces))

    edge_facelink = 0 ! array assignment
    face_edgelink = 0 ! array assignment

    ! Lookup_table(face#) assigns (arbitrarily) one out of three valid 
    ! possibilities for an edge contained by face# - and not already assigned
    ! to another face. See [Table II, Savage and Peterson, IEEE T-AP June 96]
    ! for numbering conventions.
    lookup_table(1) = 1 
    lookup_table(2) = 5 
    lookup_table(3) = 2 
    lookup_table(4) = 6  
    secondary_lookup_table(1) = 2
    secondary_lookup_table(2) = 1
    secondary_lookup_table(3) = 3 
    secondary_lookup_table(4) = 4

    DO ielem = 1,num_elements 
      DO kface = 1,4 ! local face counter
        tempface = elements(ielem)%faces(kface) ! global face           

        IF (face_edgelink(tempface).EQ.0) THEN 
          ! This face has not been processed before
          IF (faces(tempface)%free) THEN ! This is free face 

            tempedge = elements(ielem)%edges(lookup_table(kface)) ! global edge

            ! In the CBAA case, the edge may not belong to the CBAA aperture if the 
            ! face does not belong to it (MMB):
            IF ((.NOT.faces(tempface)%CBAA_aperture).AND.edges(tempedge)%CBAA_aperture) THEN
              ! this new edge is not in the aperture, because a non-aperture face cannot
              ! have more than one edge within the aperture:
              tempedge = elements(ielem)%edges(secondary_lookup_table(kface)) ! global edge
            END IF

            face_edgelink(tempface) = tempedge
            CALL count_nonzeros(SIZE(edge_facelink(tempedge,:)),& 
                 edge_facelink(tempedge,:),pos)
            IF (pos+1.GT.EST_DIM) THEN
              ! Increase array size from original estimate - double, for convenience
              TMP_DIM = EST_DIM*2 
              ! Copy into temporary variable, first re-initializing.
              ALLOCATE(temp_edge_facelink(num_edges,TMP_DIM))
              temp_edge_facelink = 0 ! array assignment
              temp_edge_facelink(:,1:EST_DIM) = edge_facelink(:,1:EST_DIM)
              ! Now re-allocate extra memory for edge_facelink and copy back.
              DEALLOCATE(edge_facelink)          
              ALLOCATE(edge_facelink(num_edges,TMP_DIM))
              edge_facelink = temp_edge_facelink ! array assignment, same size
              DEALLOCATE(temp_edge_facelink)     ! Cleanup temporary variable
              EST_DIM = TMP_DIM                  ! Increase array size 
            END IF
            edge_facelink(tempedge,pos+1) = tempface            

          END IF 
        END IF
      END DO
    END DO 
    IF (DEBUG_NUMBER_DOF) THEN
      WRITE(TESTFILE,*) 'num_faces : ', num_faces
      WRITE(TESTFILE,*) 'face_edgelink : ',face_edgelink
    END IF
    DEALLOCATE(face_edgelink) 

  END SUBROUTINE FACE_NUMBERING_LIST
!*******************************************************************************

END SUBROUTINE NUMBER_DOF

!*******************************************************************************


SUBROUTINE NUMBER_PRE_DOF
  USE eigen_analysis_data
  USE geometry
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  USE scattering_analysis_data
!*******************************************************************************
!***  (Re)-number prescribed "degrees of freedom"                           ****
!***  Since the prescribed dof's enter only in the RHS vector,              ****
!***  it is not important to minimize bandwidth, and they can  be           ****
!***  numbered consecutively for each type of dof.                          ****
!*******************************************************************************
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! DB Davidson
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version 6 August 2004 by DBD. Based on NUMBER_DOF.
!*******************************************************************************
  IMPLICIT NONE
  INTEGER(I4B) iedge,kface,jdof             ! Counters.
 

  ! Number the prescibed edges and faces for later use in assembling S and T.
  ! Note that the d.o.f. number is unique - a renumbered d.o.f. is 
  ! either an edge (type e1 or e2) or face (f1 or f2).
  ! Since the prescribed dof's enter only in the RHS vector, 
  ! numbers in each type of dof are consecutive, UNLIKE the free dof's.
  
  ! Allocate and initialize degree of freedom variables
  ! Note that the allocation is conservative: num_edges and 
  ! num_faces must be >= the respective d.o.f.'s
  ALLOCATE (renumbered_pre_e1(num_edges)) ! CT/LN  
  ALLOCATE (renumbered_pre_e2(num_edges)) ! LT/LN and LT/QN
  ALLOCATE (renumbered_pre_e3(num_edges)) ! QT/QN
  ALLOCATE (renumbered_pre_f1(num_faces)) ! LT/QN
  ALLOCATE (renumbered_pre_f2(num_faces)) ! LT/QN
  ALLOCATE (renumbered_pre_f3(num_faces)) ! QT/QN
  ALLOCATE (pre_dof_type(3*num_edges+3*num_faces))    
  renumbered_pre_e1 = 0 ! Array initialization (0 implies non-existent)
  renumbered_pre_e2 = 0 ! Ditto
  renumbered_pre_e3 = 0 ! Ditto
  renumbered_pre_f1 = 0 ! Ditto
  renumbered_pre_f2 = 0 ! Ditto
  renumbered_pre_f3 = 0 ! Ditto
  pre_dof_type      = 0 ! Ditto

  jdof  = 0 ! initialise

  PRE_ASSIGN_EDGE_LOOP: DO iedge = 1,num_edges
    IF (edges(iedge)%Dirichlet) THEN ! D.o.f. associated with prescribed edge
      ! Always do the following: CT/LN 
	  jdof  = jdof  + 1 
      renumbered_pre_e1(iedge) = jdof
      pre_dof_type(jdof)  = 1
      IF(edges(iedge)%order.GE.2.OR.&
	  	.NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.1) THEN
        ! LT/LN and LT/QN
	    jdof  = jdof  + 1
        renumbered_pre_e2(iedge) = jdof  
        pre_dof_type(jdof)  = 2                
      END IF
      IF(edges(iedge)%order.GE.3.OR.&
	    .NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.2) THEN
        ! QT/QN and QT/CuN
        jdof  = jdof  + 1
        renumbered_pre_e3(iedge) = jdof  
        pre_dof_type(jdof)  = 5
	  END IF                
    END IF
  END DO PRE_ASSIGN_EDGE_LOOP

  PRE_ASSIGN_FACE_LOOP: DO kface = 1,num_faces
    ! Now assign the dofs associated with the faces 
    IF (faces(kface)%Dirichlet) THEN ! D.o.f. associated with prescribed face
        ! Note: CT/LN and LT/LN elements have no face-based dof's
        IF (faces(kface)%order.GE.2) THEN
          jdof  = jdof  + 1 
          renumbered_pre_f1(kface) = jdof  
          pre_dof_type(jdof) = 3                              
          jdof  = jdof  + 1 
          renumbered_pre_f2(kface) = jdof
          pre_dof_type(jdof) = 4
        END IF
        IF (((.NOT.faces(kface)%mixed).AND.(faces(kface)%order.EQ.2)) &
		  .OR.((faces(kface)%order.GE.3))) THEN
		  jdof  = jdof  + 1 
          renumbered_pre_f3(kface) = jdof
          pre_dof_type(jdof) = 6
        END IF
      END IF
  END DO PRE_ASSIGN_FACE_LOOP

  pre_dof = jdof
  
  IF (DEBUG_SYSTEM_MATRIX) THEN
    CALL DEBUG_PRE_DOF
  END IF

  
  ! Following section of code is NOT for production! It tests the interpolation function 
  ! by essentially setting ALL dof's as prescribed....
  IF (TEST_INTERPOLATE_FIELD) THEN
    jdof  = 0 ! initialise
    DO iedge = 1,num_edges
      jdof  = jdof  + 1 
      renumbered_pre_e1(iedge) = jdof
      pre_dof_type(jdof)  = 1
      IF(edges(iedge)%order.GE.2.OR.&
	  	.NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.1) THEN
	    jdof  = jdof  + 1
        renumbered_pre_e2(iedge) = jdof  
        pre_dof_type(jdof)  = 2                
      END IF
      IF(edges(iedge)%order.GE.3.OR.&
	    .NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.2) THEN
        ! QT/QN and QT/CuN
        jdof  = jdof  + 1
        renumbered_pre_e3(iedge) = jdof  
        pre_dof_type(jdof)  = 5
	  END IF                
    END DO 
  
    DO kface = 1,num_faces
    ! Now assign the dofs associated with the faces 
        ! Note: CT/LN and LT/LN elements have no face-based dof's
      IF (faces(kface)%order.GE.2) THEN
          jdof  = jdof  + 1 
          renumbered_pre_f1(kface) = jdof  
          pre_dof_type(jdof) = 3                              
          jdof  = jdof  + 1 
          renumbered_pre_f2(kface) = jdof
          pre_dof_type(jdof) = 4
      END IF
      IF (((.NOT.faces(kface)%mixed).AND.(faces(kface)%order.EQ.2)) &
             .OR.((faces(kface)%order.GE.3))) THEN
		  jdof  = jdof  + 1 
          renumbered_pre_f3(kface) = jdof
          pre_dof_type(jdof) = 6
      END IF
    END DO 
    pre_dof = jdof
	! Now set all the dof's to NOT FREE:
	dof= 0 
    renumbered_e1(:) = 0
    renumbered_e2(:) = 0
    renumbered_e3(:) = 0
    renumbered_f1(:) = 0
    renumbered_f2(:) = 0
    renumbered_f3(:) = 0
    

  END IF


CONTAINS

  SUBROUTINE DEBUG_PRE_DOF
    IMPLICIT NONE
!*******************************************************************************
! List prescribed degrees of freedom - for debugging only
! DBD 6 August 2004
!*******************************************************************************
    INTEGER(I4B) iedge,iface,jdof,kedge,kface

    IF (HIERARCHAL_ORDER.GE.2.AND..NOT.MIXED_ORDER_FLAG) THEN
	  STOP 'IE: DEBUG_PRE_DOF: Only implemented to LT/QN'
	END IF 

    WRITE(FILEOUT,'(//A)') '****** PRESCRIBED "DEGREES OF FREEDOM" **************'
    WRITE(FILEOUT,'(/A)') 'Requested by flag DEBUG_SYSTEM_MATRIX'

    WRITE(FILEOUT,'(/A)') 'Prescribed degrees of freedom by type:'
    WRITE(FILEOUT,'(/A,7X,A,4X,A/)') 'Type','Prescribed D.o.f.','Old global edge/face#'               

    DO iedge = 1,num_edges
      IF (renumbered_pre_e1(iedge).NE.0) THEN
        WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'e1',renumbered_pre_e1(iedge),iedge
      END IF
      IF(edges(iedge)%order.GE.2.OR.&
	  	.NOT.edges(iedge)%mixed.AND.edges(iedge)%order.GE.1) THEN
        IF (renumbered_pre_e2(iedge).NE.0) THEN ! LT/LN and higher 
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'e2',renumbered_pre_e2(iedge),iedge
		END IF
      END IF
    END DO  
    DO iface = 1,num_faces
      SELECT CASE(HIERARCHAL_ORDER)
      CASE(1) ! CT/LN     
        ! No action
      CASE(2) ! LT/QN
        IF (renumbered_pre_f1(iface).NE.0) THEN
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'f1',renumbered_pre_f1(iface),iface
        END IF
        IF (renumbered_pre_f2(iface).NE.0) THEN
          WRITE(FILEOUT,'(A,3X,I8,2X,I8)') 'f2',renumbered_pre_f2(iface),iface
        END IF
      CASE DEFAULT
        STOP 'Error in routine DEBUG_PRE_DOF - unimplemented hierarchal order'
      END SELECT
    END DO

    WRITE(FILEOUT,'(//A)') 'Prescribed "Degrees of freedom" in sequence'
    WRITE(FILEOUT,'(/A,7X,A,4X,A/)') 'Type','D.o.f.'    

    DO jdof = 1,pre_dof
      SELECT CASE(pre_dof_type(jdof))
      CASE(1) ! e1     
        WRITE(FILEOUT,'(A,3X,I8)') 'e1',jdof              
      CASE(2) ! e2     
        WRITE(FILEOUT,'(A,3X,I8)') 'e2',jdof  
      CASE(3) ! e3     
        WRITE(FILEOUT,'(A,3X,I8)') 'f1',jdof  
      CASE(4) ! e4     
        WRITE(FILEOUT,'(A,3X,I8)') 'f2',jdof
      CASE DEFAULT
        STOP 'IE: Internal error in routine DEBUG_PRE - incorrect d.o.f. type'
      END SELECT
    END DO         

  END SUBROUTINE DEBUG_PRE_DOF
!*******************************************************************************

END SUBROUTINE NUMBER_PRE_DOF

!*******************************************************************************

! Changed DBD 31 March 2003
 SUBROUTINE U_MAKE_HIERARCHAL(element_num,local_face_num,element_order,port_or_ABC_num,&
                              normal,Us)
! End changed DBD 31 March 2003
  USE geometry
  USE nrtype
  USE feminterface, ONLY: GAUSS_QUAD_VBF,GAUSS_QUAD_VBF_SPC 
  USE math_tools, ONLY: CROSS_PRODUCT
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs the U matrix for the surface integration of 
! the incident field weighted by the vector basis functions  over the port
! surfaces, for the given element and local face number. 
! [DBD correction to documentation:] Currently implemented up to complete 2nd order. 
! The [Us] matrix returned is the [b^s] matrix in the notation of: 
! J-M Jin, "The Finite element method in electromagnetics", Wiley 1993, 
! eqn. (8.87) p. 265. 
!
! The same code is also used to compute the forcing vector f_i 
! (eqn.12.10,p.531, Jin 2nd edn) for time-domain analysis, for either the 
! total field formulation (in which case the integral is over the outer boundary)
! or for the scattered field formulation (in which case the integral is over the 
! fictitious internal boundary demarcated the homogeneous exterior from inhomogeneous 
! scatterer. The necessary discrimination is done within GAUSS_QUAD_VBF.
!
! The same code is also used to compute the contribution of the incident 
! field to the boundary integral for frequency-domain open region problems.
! Both the scattered and total field formulations are supported. Again, the  
! necessary discrimination is done within GAUSS_QUAD_VBF.
! (A complex valued version of GAUSS_QUAD_VBF is called in this case).
! 
! The formulation used to compute the required dot products  of the vector 
! basis functions is that of:
! [S&P] Savage and Peterson, "Higher-order vector finite elements for tetrahedral 
! cells", IEEE MTT, June 96, pp. 874-879. 
!
! with extensions: 
! [DBD] DB Davidson, Lab Notebook FEM Vol III, p 39-41, March 2000.  
!
! The S&P elements can either be computed using closed-form expressions, or 
! using cubature.  
!
! The code also offers the Andersen & Volakis elements, implemented using
! cubature: LS Andersen,JL Volakis, "Hierarhical tangential vector finite 
! elements for tetrahedra", IEEE M&GW Lttrs, March 1998, p.127--129.
!
! The surface vector basis functions are treated as in B_MAKE_HIERARCHAL_ANALYTIC.
!
! [Correction by DBD 26 March 2003] 
! The (outward directed) unit normal vector is general.
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! D B Davidson.
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! 24 May 2000: First release. DBD.
! 13 Jul 2000: Extended to included LT/QN elements. DBD.
! 16 Feb 2001  Cubature and A&V elements added. DBD. 
! 03 Mar 2002: Support for mixed order, and up to QT/QN elements added. DBD.
! 17 Mar 2002: Corrected the e2 contribution to the LT/LN element case,
!              as specified by DBD via email. MMB. 
! 31 Mar 2003: k0 removed from argument list (unused). Minor changed 
!              to permit this routine to also handle time-domain anaylysis.
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
! Output data is the elemental [Us] vector - {b^s} in [eqn.(8.86),p. 265.,Jin ].
! The full system matrix assembly is done elsewhere, taking the BC's into 
! account. 
! The storage convention used is the following:
!
! | Ue1|
! | Ue2|
! | Uf1|
! | Uf2|
! | Ue3|
! | Uf3|
!
! where e1, e2, f1 and f2 refer to the edge and face based functions.
! See [S&P] and [DBD]. 
!
! Note the [S&P] include the length of the edge in the H_0(curl) element
! and it is (optionally) included here. See notes in S_AND_T_MAKE_HIERARCHAL.
! 
! The numbering convention of nodes, edges and faces follows [S&P] in this
! sub-routine [Table II,S&P]. 
!
!*******************************************************************************
! Changed DBD 31 March 2003  
  INTEGER(I4B), INTENT(IN) :: element_num, element_order, & 
                              local_face_num, port_or_ABC_num
! Changed DBD 29 March 2001 and 31 March 2003
!  REAL(SP), INTENT(IN) ::  k0 !,x0
! Changed DBD 29 March 2001
!  REAL(SP), INTENT(IN), OPTIONAL ::  y0
  REAL(SP), INTENT(IN), DIMENSION(3) :: normal ! Outward unit normal on surface

  COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: Us
  COMPLEX(SPC), DIMENSION(3) :: Ue1,Ue2,Ue3
  COMPLEX(SPC) :: Uf1,Uf2,Uf3
  INTEGER(I4B) :: iedge              ! Misc counters.
  INTEGER(I4B) :: i1,i2,i3           ! Misc
  INTEGER(I4B), DIMENSION(2) :: i_edge_nodes 
  INTEGER(I4B), DIMENSION(3) :: temp_global_face_nodes
  INTEGER(I4B), DIMENSION(3) :: temp_local_face_nodes
  REAL(SP), DIMENSION(3) :: ell      ! edge lengths

  IF ( ABS(SQRT(DOT_PRODUCT(normal,normal))-1.0_SP).GE.EPS) THEN
    STOP 'IE: U_MAKE_HIERARCHAL called with invalid unit normal.'
  END IF

  Us = ZERO_C ! NB - this is a precaution, to ensure that 
              ! non-existent higher order terms are zero
			  ! for lower-order elements. 
  ! Find the lengths of the element edges
  temp_global_face_nodes(1:3) = GLOBAL_FACENODES(element_num,local_face_num)
  ell(1) = T_LENGTH(temp_global_face_nodes(1),temp_global_face_nodes(2))
  ell(2) = T_LENGTH(temp_global_face_nodes(1),temp_global_face_nodes(3))
  ell(3) = T_LENGTH(temp_global_face_nodes(2),temp_global_face_nodes(3))
  DO iedge=1,3
    i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
    i1 = i_edge_nodes(1)
    i2 = i_edge_nodes(2)
	IF(FD_SCAT_ANALYSIS) THEN
      Ue1(iedge) = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                            1,i1,i2,normal,ell=ell(iedge))  
    ELSE 
      Ue1(iedge) = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                            1,i1,i2,normal,ell=ell(iedge))
    END IF  
  END DO
  Us(1:3) = Ue1(1:3)

  IF(MAX_ORDER(element_num,local_face_num).GE.2.OR.&
    .NOT.MIXED_ORDER(element_num,local_face_num).AND.& 
    MAX_ORDER(element_num,local_face_num).GE.1) THEN
    DO iedge=1,3
      i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
      i1 = i_edge_nodes(1)
      i2 = i_edge_nodes(2)
      IF(FD_SCAT_ANALYSIS) THEN
        Ue2(iedge) = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                                2,i1,i2,normal,ell=ell(iedge))  
      ELSE 
	    Ue2(iedge) = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                                2,i1,i2,normal,ell=ell(iedge))  
      END IF
    END DO
    Us(4:6) = Ue2(1:3)   ! <- this is the correction!
  END IF

  IF(MAX_ORDER(element_num,local_face_num).GE.2) THEN
    temp_local_face_nodes(1:3) = LOCAL_FACENODES(local_face_num)
    i1 = temp_local_face_nodes(1)
    i2 = temp_local_face_nodes(2)
    i3 = temp_local_face_nodes(3)
    IF(FD_SCAT_ANALYSIS) THEN
      Uf1 = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                            3,i1,i2,normal,i3=i3)  
      Uf2 = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                            4,i1,i2,normal,i3=i3)  
    ELSE 
      Uf1 = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                            3,i1,i2,normal,i3=i3)  
      Uf2 = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                            4,i1,i2,normal,i3=i3)  
    END IF
	Us(7) = Uf1
    Us(8) = Uf2
  END IF

  IF(MAX_ORDER(element_num,local_face_num).GE.3.OR.&
    .NOT.MIXED_ORDER(element_num,local_face_num).AND.& 
	MAX_ORDER(element_num,local_face_num).GE.2) THEN
    DO iedge=1,3
      i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
      i1 = i_edge_nodes(1)
      i2 = i_edge_nodes(2)
      IF(FD_SCAT_ANALYSIS) THEN
        Ue3(iedge) = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                              5,i1,i2,normal,ell=ell(iedge))
      ELSE   
        Ue3(iedge) = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                              5,i1,i2,normal,ell=ell(iedge))  

      END IF
    END DO
    i1 = temp_local_face_nodes(1)
    i2 = temp_local_face_nodes(2)
    i3 = temp_local_face_nodes(3)
    IF(FD_SCAT_ANALYSIS) THEN
      Uf3 = GAUSS_QUAD_VBF_SPC(port_or_ABC_num,element_num,local_face_num,& 
                            6,i1,i2,normal,i3=i3)  
    ELSE 
	   Uf3 = GAUSS_QUAD_VBF(port_or_ABC_num,element_num,local_face_num,& 
                            6,i1,i2,normal,i3=i3)  
    END IF
    Us(9:11) = Ue3(1:3)
	Us(12) = Uf3
  END IF

END SUBROUTINE U_MAKE_HIERARCHAL  
!*******************************************************************************



SUBROUTINE ASSIGN_FACE_AND_EDGE_ORDERS(init_flag)
  USE geometry
  USE problem_info
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Assigns values to the variables <faces%order> and <edges%order> according to
! the values <elements%order>, i.e. the hierarchal orders of the elements.
! <init_flag> indicates whether this routine must initialise the values
! <elements%order> as well.
! Also assigns values to the variables <faces%mixed> and <edges%mixed> according to
! the values <elements%mixed>
!
! Created: MMB 2001-10-08.
! Extended: DBD 2002-02-23 to include mixed-order flag on edges and faces.
! 2002-03-17: Corrections to mixed assignment, 'a 2-mixed entity overrides
!             a 1-full entity'. MMB.
!*******************************************************************************
  LOGICAL(LGT), INTENT(IN) :: init_flag

  INTEGER(I4B) :: ielem,iface,globalface,iedge
  INTEGER(I4B), DIMENSION(3) :: globaledges

  ! Initialise the element orders if requested:
  IF (init_flag) THEN 
    elements%order = HIERARCHAL_ORDER
	elements%mixed = MIXED_ORDER_FLAG ! Added DBD 21 Feb 02
  END IF

  ! Initialise the variables to be assigned, with the lowest order:
  faces%order = 1      ! matrix assignment
  edges%order = 1      !      ditto.
  faces%mixed = .TRUE. ! matrix assignment (Additions: DBD 21 Feb 02)
  edges%mixed = .TRUE. !      ditto.

  ! Assign the orders of the faces and edges according to the orders
  ! of the elements they belong to. The order of a face/edge shared by more 
  ! than one element will be that of the highest order element. 
  ! Similarly, polynomial completeness over-rides mixed order.
  ELEMENT_LOOP: DO ielem = 1,num_elements
    FACE_LOOP: DO iface = 1,4

      globalface = elements(ielem)%faces(iface)

      ! Assign the face order and mixed/complete status:
      IF (elements(ielem)%order.GT.faces(globalface)%order) THEN
        faces(globalface)%order = elements(ielem)%order
        faces(globalface)%mixed = elements(ielem)%mixed
      ELSE IF (elements(ielem)%order.EQ.faces(globalface)%order) THEN
        faces(globalface)%mixed = faces(globalface)%mixed.AND.elements(ielem)%mixed
      END IF
      
      globaledges = LOCAL_FACEEDGES(iface)
      globaledges = elements(ielem)%edges(globaledges)

      EDGE_LOOP: DO iedge = 1,3
        
        ! Assign the edge order and mixed/complete status:
        IF (elements(ielem)%order.GT.edges(globaledges(iedge))%order) THEN
          edges(globaledges(iedge))%order = elements(ielem)%order
          edges(globaledges(iedge))%mixed = elements(ielem)%mixed
        ELSE IF (elements(ielem)%order.EQ.edges(globaledges(iedge))%order) THEN
          edges(globaledges(iedge))%mixed = edges(globaledges(iedge))%mixed.AND.elements(ielem)%mixed
        END IF
      
      END DO EDGE_LOOP
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE ASSIGN_FACE_AND_EDGE_ORDERS
!*******************************************************************************

