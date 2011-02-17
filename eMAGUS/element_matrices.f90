MODULE B_matrix
  USE nrtype
  USE basis_function
  IMPLICIT NONE

CONTAINS

  SUBROUTINE B_MAKE_HIERARCHAL(elem,local_face,normal,Bs,int_request)
    USE geometry
    USE math_tools, ONLY: CROSS_PRODUCT,VECTOR_LENGTH
    USE quad_tables
    IMPLICIT NONE
    !*******************************************************************************
    ! SUBROUTINE DESCRIPTION
    !*******************************************************************************
    ! This subroutine constructs the B matrix for the surface integration of 
    ! the vector basis functions required over the ports in the 
    ! guided wave analysis for the given element and local face number. 
    ! It is performs the same function as B_MAKE_HIERARCHAL_ANALYTIC. 
    ! It has now replaced the analytical integration in that routine. 
    ! This routine calculates the elemental matrix defined by equation 8.86 
    ! [Jin], 8.87 [Jin 2nd], Quadrature is used.
    ! This routine may also be used without modification for the FE time domain analysis, 
    ! and with an extension, for 2nd order ABC computation. 
    !
    ! It may also be used to evalute the surface integral appearing in the first
    ! order ABC, eqn. 9.61 [Jin 2nd], and by setting the optional flag
    ! int_request to 1, to compute the 2nd surface integral appearing
    ! in the second order ABC, eqn. 9.71 [Jin 2nd]. (Leaving this flag out
    ! defaults to the former), and with int_request set to 2, both the 2nd 
    ! and 3rd integrals are evaluated.
    !
    ! Currently implemented are CT/LN, LT/LN, LT/QN and QT/QN elements (the last NOT for ABCs!). 
    ! (Various higher order elements are supported, see S_AND_T_MAKE_HIERARCHAL
    ! for more detail,. 
    ! The [Bs] matrix returned  is the [B^s] matrix in the notation of: 
    ! [Jin ]J-M Jin, "The Finite element method in electromagnetics", Wiley 1993, 
    ! eqn. (8.86) p. 265.
    ! This routine returns the matrix without the \gamma scaling factor. 
    !
    ! The formulation used to compute the required dot products  of the vector 
    ! basis functions is an extension of:
    ! [S&P] Savage and Peterson, "Higher-order vector finite elements for tetrahedral 
    ! cells", IEEE MTT, June 96, pp. 874-879. 
    !
    ! with extensions: 
    !
    ! [DBD] DB Davidson, Lab Notebook FEM Vol III, p 39-41, March 2000.
    !
    ! Note that the vector basis functions in this routine are actually the 
    ! SURFACE basis functions, computed as:
    ! \vec{S}^s_i = \hat{N} \times \vec{N}^s_i
    ! The 3D simplex coordinates, and the gradients thereof, are used throughout 
    ! this routine, although the surface basis functions are written in terms 
    ! of 2D simplex coordinates. As far as the simplex coordinates go, one of them
    ! is zero; as regards the gradient of the simplex coordinates, it may be shown
    ! that taking \hat{N} \times the (3D) gradient (ie the component in the
    ! plane of the triangle) yields the correct 2D result.
    !
    !*******************************************************************************
    ! AUTHOR
    !*******************************************************************************
    ! D B Davidson.
    !
    !*******************************************************************************
    ! Last revised:
    !*******************************************************************************
    ! 12 Feb 2002: Created by MMB from B_MAKE_HIERARCHAL_ANALYTIC
    ! 19 Feb 2002: Changed to default routine. MMB
    ! 28 Feb 2002: Extended to cater for possible elements beyond LT/QN. DBD.
    ! 07 Mar 2002: Documentation updated. DBD.
    ! 01 Apr 2003: Ditto. DBD.
    ! 20 Jul 2004: Extended to compute 2nd surface integral appearing
    !              in 2nd order ABC.
    ! 29 Sep 2004: Ditto for 3rd integral.
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
    ! Output data is the elemental B matrix, [eqn.(8.86),p. 265.,Jin ]
    ! The full system matrix assembly is done elsewhere, taking the BC's into 
    ! account. 
    !
    ! The storage convention used is the following:
    !
    ! | Be1e1 Be1e2 Be1f1 Be1f2 Be1e3 Be1f3 |
    ! | Be2e1 Be2e2 Be2f1 Be2f2 Be2e3 Be2f3 |
    ! | Bf1e1 Bf1e2 Bf1f1 Bf1f2 Bf1e3 Bf1f3 |
    ! | Bf2e1 Bf2e2 Bf2f1 Bf2f2 Bf2e3 Bf2f3 |
    ! | Be3e1 Be3e2 Be3f1 Be3f2 Be3e3 Be3f3 |
    ! | Bf3e1 Bf3e2 Bf3f1 Bf3f2 Bf3e3 Bf3f3 |
    !
    ! where e1, e2, e3, f1, f2 and f3 refer to the edge and face based functions.
    ! See [S&P], [DBD] and [Webb99].
    !
    ! Note the [S&P] include the length of the edge in the H_0(curl) element
    ! and it (optionally) is included here. See notes in S_AND_T_MAKE_HIERARCHAL.
    ! 
    ! The numbering convention of nodes, edges and faces follows [S&P] in this
    ! sub-routine [Table II,S&P]. 
    !
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: elem,local_face
    REAL(SP), DIMENSION(3), INTENT(IN) :: normal
    COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE), INTENT(OUT) :: Bs
    INTEGER(I4B), INTENT(IN), OPTIONAL :: int_request
    INTEGER(I4B) ::  int_value
    REAL(SP),DIMENSION(3) :: ell,u_m
    REAL(SP),DIMENSION(3) :: xi_hat,eta_hat,zeta_hat ! Local rectangular coordinate system
    INTEGER(I4B) :: num_qpts,num_qpts_line ! number of quadrature points
    INTEGER(I4B) :: iquad,iedge,ifunc,mm,nn,row,col ! counters
    REAL(SP),DIMENSION(4,4) :: vertmat
    REAL(SP),DIMENSION(4) :: r_vec,r_vec_xip,r_vec_xim,r_vec_etap,r_vec_etam,s_vec
    INTEGER(I4B),DIMENSION(3) :: facenodes
    INTEGER(I4B), DIMENSION(3) :: face_edges
    INTEGER(I4B),DIMENSION(2) :: edgenodes
    REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: vbf_values,vbf_values_xip, vbf_values_xim, & 
         vbf_values_etap,vbf_values_etam,vbf_tan
    REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE) :: vbf_norm_values, div_vbf_values
    REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Bs_large
    REAL(SP) delta ! Delta used for central difference approximation


    ! Added DBD 20 Jul 04
    IF (PRESENT(int_request)) THEN
       int_value = int_request
    ELSE
       int_value = 0
    END IF
    ! End added DBD 20 Jul 04

    ! Check that normal is valid:
    IF ( ABS(SQRT(DOT_PRODUCT(normal,normal))-1.0_SP).GE.EPS) THEN
       STOP 'IE: B_MAKE_HIERARCHAL called with invalid unit normal.'
    END IF
    ! Check that a defined operation has been requested.

    ! Added DBD 18 Nov 04
    SELECT CASE(int_value) 
    CASE(0:3) ! "Default" case.
       CONTINUE 
    CASE DEFAULT
       STOP 'IE: Error in B_MAKE_HIERARCHAL, undefined type of integral requested'
    END SELECT
    ! End added DBD 18 Nov 04


    num_qpts = 7  ! set the number of surface quadrature points - 5th order complete
    delta = 0.02_SP*SUM(T_LENGTHS(elem))/6.0_SP ! Sets the delta used to compute the central difference to 1/50 of the 
    ! average edge lengths.

    ! Establish vertices matrix for repeated use:
    CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vertmat)
    facenodes = LOCAL_FACENODES(local_face)

    Bs_large = 0.0 ! matrix assignment

    s_vec = 0.0 ! vector assignment
    QUAD_LOOP: DO iquad = 1,num_qpts

       ! Calculate the xyz coords of current quad point:
       s_vec(facenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
       r_vec = MATMUL(vertmat,s_vec)

       ! Evaluate basis functions (or curl thereof):
       ! Changed DBD 20 Jul 04 and 29 Sep 04
       SELECT CASE(int_value)
       CASE(0) ! "Default" case. 
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values)
          ! Cross with the normal:        
          DO ifunc = 1,ELEM_TET_MATRIX_SIZE
             vbf_values(ifunc,1:3) = CROSS_PRODUCT(normal,vbf_values(ifunc,1:3))
          END DO
          ! Add to <Bs_large>:
          DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
             DO nn = 1,ELEM_TET_MATRIX_SIZE ! column loop
                Bs_large(mm,nn) = Bs_large(mm,nn) + &
                     quad_tri_rules(num_qpts)%rule(iquad,4) * &
                     SUM(vbf_values(mm,1:3)*vbf_values(nn,1:3))
             END DO
          END DO
       CASE(1) ! Normal component of curl of vbf.
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),1,vbf_values) 
          ! DOT with the normal:        
          DO ifunc = 1,ELEM_TET_MATRIX_SIZE
             vbf_norm_values(ifunc) = DOT_PRODUCT(normal,vbf_values(ifunc,1:3))
          END DO
          ! Add to <Bs_large>:
          DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
             DO nn = 1,ELEM_TET_MATRIX_SIZE ! column loop
                Bs_large(mm,nn) = Bs_large(mm,nn) + &
                     quad_tri_rules(num_qpts)%rule(iquad,4) * &
                     vbf_norm_values(mm)*vbf_norm_values(nn)
             END DO
          END DO
       CASE(2) ! Surface divergence term. Surface divergence is computed numerically using 
          ! central differences. Local coordinate system is xi,eta,zeta, with xi chosen
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
          ! CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values) ! Not actually needed
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_xip(1:3),0,vbf_values_xip) 
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_xim(1:3),0,vbf_values_xim) 
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_etap(1:3),0,vbf_values_etap) 
          CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec_etam(1:3),0,vbf_values_etam) 
          DO ifunc = 1,ELEM_TET_MATRIX_SIZE
             ! Find surface divergence using central difference approximation for the two in-plane
             ! components
             div_vbf_values(ifunc) = DOT_PRODUCT((vbf_values_xip(ifunc,1:3) - vbf_values_xim(ifunc,1:3)),xi_hat) + & 
                  DOT_PRODUCT((vbf_values_etap(ifunc,1:3) - vbf_values_etam(ifunc,1:3)),eta_hat)
             div_vbf_values(ifunc) = div_vbf_values(ifunc)/(2.0_SP*delta)  
          END DO
          ! Add to <Bs_large>:
          DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
             DO nn = 1,ELEM_TET_MATRIX_SIZE ! column loop
                Bs_large(mm,nn) = Bs_large(mm,nn) + &
                     quad_tri_rules(num_qpts)%rule(iquad,4) * &
                     div_vbf_values(mm)*div_vbf_values(nn)
             END DO
          END DO
       END SELECT
    END DO QUAD_LOOP


!!! Following is really dead code, since the implementation is apparently 
!!! incorrect (it implements the surface div squared term  via the square of the line 
!!! integral, which does not appear to be correct); however, the results are
!!! plausible. Thus it has been left for possible re-use in future investigations.


    num_qpts_line = 3  ! set the number of line quadrature points - 5th order complete
    SELECT CASE(int_value)
       ! Changed DBD 18 Nov 04
    CASE(3) ! Surface divergence of vbf.
       s_vec = 0.0 ! vector assignment
       face_edges = LOCAL_FACEEDGES(local_face)
       DO iedge = 1,3         ! This algorithm works on an edge-by-edge basis.
          edgenodes = LOCAL_EDGENODES(face_edges(iedge))
          ! Get edge unit vectors.
          ell(:) = EDGE_UNIT_VECTOR(edgenodes(1),edgenodes(2),elem)
          ! Compute u_m vectors on edge.
          u_m(:) = CROSS_PRODUCT(ell,normal)
          ! Integrate in-plane (tangential) component of appropriate vbf along each edge
          LINE_QUAD_LOOP: DO iquad = 1,num_qpts_line
             ! Calculate the xyz coords of current quad point:
             s_vec(edgenodes) = quad_line_rules(num_qpts_line)%rule(iquad,1:2)
             r_vec = MATMUL(vertmat,s_vec)
             CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values)
             ! Get tangential components in the plane of the face ([Jin, p.354]):        
             DO ifunc = 1,ELEM_TET_MATRIX_SIZE
                vbf_tan(ifunc,:) = - CROSS_PRODUCT(normal,CROSS_PRODUCT(normal,vbf_values(ifunc,:)))
                vbf_norm_values(ifunc) = DOT_PRODUCT(vbf_tan(ifunc,:),u_m)
             END DO

             ! Add to <Bs_large>:
             DO mm = 1,ELEM_TET_MATRIX_SIZE ! row loop
                DO nn = 1,ELEM_TET_MATRIX_SIZE ! column loop
                   Bs_large(mm,nn) = Bs_large(mm,nn) + &
                        quad_line_rules(num_qpts_line)%rule(iquad,3) * &
                        vbf_norm_values(mm)*vbf_norm_values(nn)
                END DO
             END DO
             Bs_large = - Bs_large ! Test code, works better with minus here. 
          END DO LINE_QUAD_LOOP

       END DO
    END SELECT
    ! End changed DBD 18 Nov 04


    ! Scale by face area (final step in quadrature):
    Bs_large = FACE_AREA(elem,local_face) * Bs_large

    ! Extract <Bs> fron <Bs_large>:
    DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row loop
       row = LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face,mm)
       DO nn = 1,ELEM_TRI_MATRIX_SIZE ! column loop
          col = LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face,nn)
          Bs(mm,nn) = Bs_large(row,col)
       END DO
    END DO

  END SUBROUTINE B_MAKE_HIERARCHAL
  !*******************************************************************************


  ! Following group of routines will probably eventually be removed. 
  ! Retained for speed tests to compare quadrature to analytical 
  ! evaluation. 

  SUBROUTINE B_MAKE_HIERARCHAL_ANALYTIC(element_num,local_face_num,element_order,&
       element_mixed_order,normal,Bs)
    USE geometry
    USE math_tools, ONLY: CROSS_PRODUCT
    USE problem_info
    USE quad_tables
    IMPLICIT NONE
    !*******************************************************************************
    ! SUBROUTINE DESCRIPTION
    !*******************************************************************************
    ! This subroutine constructs the B matrix for the surface integration of 
    ! the vector basis functions required over the ports in the 
    ! guided wave analysis for the given element and local face number. 
    ! It is now retained for verification purposes only; B_MAKE_HIERARCHAL performs
    ! the same functionality using quadrature. 
    ! Currently implemented are CT/LN and the Savage LT/QN elements. 
    ! The [Bs] matrix returned  is the [B^s] matrix in the notation of: 
    ! [Jin ]J-M Jin, "The Finite element method in electromagnetics", Wiley 1993, 
    ! eqn. (8.86) p. 265. 
    !
    ! The formulation used to compute the required dot products  of the vector 
    ! basis functions is an extension of:
    ! [S&P] Savage and Peterson, "Higher-order vector finite elements for tetrahedral 
    ! cells", IEEE MTT, June 96, pp. 874-879. 
    !
    ! with extensions: 
    !
    ! [DBD] DB Davidson, Lab Notebook FEM Vol III, p 39-41, March 2000.
    !
    ! The S&P elements can either be computed using closed-form expressions, or 
    ! using cubature.  
    !
    ! The code also offers the Andersen & Volakis elements, implemented using
    ! cubature: LS Andersen,JL Volakis, "Hierarhical tangential vector finite 
    ! elements for tetrahedra", IEEE M&GW Lttrs, March 1998, p.127--129.
    !
    ! Note that the vector basis functions in this routine are actually the 
    ! SURFACE basis functions, computed as:
    ! \vec{S}^s_i = \hat{N} \times \vec{N}^s_i
    ! The 3D simplex coordinates, and the gradients thereof, are used throughout 
    ! this routine, although the surface basis functions are written in terms 
    ! of 2D simplex coordinates. As far as the simplex coordinates go, one of them
    ! is zero; as regards the gradient of the simplex coordinates, it may be shown
    ! that taking \hat{N} \times the (3D) gradient (ie the component in the
    ! plane of the triangle) yields the correct 2D result.
    !
    ! Although ports are presently constrained to lie in  planes of constant z 
    ! in GW_SYSMAT, this routine does not make this assumption. The (outward
    ! directed) unit normal vector is general.
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
    ! 17 Jul 2000: LT/QN element support added. DBD.
    ! 16 Feb 2001: Cubature and A&V elements added. DBD.
    !*******************************************************************************
    ! INPUT 
    !*******************************************************************************
    ! Input data is the element number, order and element geometry. The latter
    ! is passed via MODULE geometry.
    !
    !*******************************************************************************
    ! OUTPUT 
    !*******************************************************************************
    ! Output data is the elemental B matrix, [eqn.(8.86),p. 265.,Jin ]
    ! The full system matrix assembly is done elsewhere, taking the BC's into 
    ! account. 
    !
    ! The storage convention used is the following:
    !
    ! | Be1e1 Be1e2 Be1f1 Be1f2 |
    ! | Be2e1 Be2e2 Be2f1 Be2f2 |
    ! | Bf1e1 Bf1e2 Bf1f1 Bf1f2 |
    ! | Bf2e1 Bf2e2 Bf2f1 Bf2f2 |
    !
    ! where e1, e2, f1 and f2 refer to the edge and face based functions.
    ! See [S&P] and [DBD]. 
    !
    ! Note the [S&P] include the length of the edge in the H_0(curl) element
    ! and it (optionally) is included here. See notes in S_AND_T_MAKE_HIERARCHAL.
    ! 
    ! The numbering convention of nodes, edges and faces follows [S&P] in this
    ! sub-routine [Table II,S&P]. 
    !
    !*******************************************************************************
!!!
!!! Interface Variables
!!!
    INTEGER(I4B), INTENT(IN) :: element_num, element_order, local_face_num
    REAL(SP), DIMENSION(3), INTENT(IN) :: normal ! Outward unit normal on surface
    COMPLEX(SPC), INTENT(OUT), DIMENSION(8,8) :: Bs
    ! Dimensioned up to and including H_1(curl) (LT/QN) elements.
    LOGICAL(LGT), INTENT(IN) :: element_mixed_order
!!!
!!! Local Variables
!!!
    INTEGER(I4B) :: inode,jnode        ! Loop counters
    INTEGER(I4B) :: knode,lnode        !  "
    INTEGER(I4B), DIMENSION(3) :: temp_face_nodes
    REAL(SP), DIMENSION(4) :: a,b,c,d  ! gradient terms: See [eqn(6),S&P] 
    REAL(SP), DIMENSION(4,3) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4) :: b_s,c_s,d_s  
    ! gradient terms for surface expansion 
    ! following normal-cross product vector
    ! basis function.
    REAL(SP), DIMENSION(4,4) :: M      ! integration matrix M cf [eqn(15),S&P]
    REAL(SP), DIMENSION(3) :: ell      ! edge lengths
    REAL(SP), DIMENSION(4,4,4) :: N    ! integration matrix N [eqn(46),S&P]
    REAL(SP), DIMENSION(4,4,4,4) :: P  ! integration matrix P [eqn(47),S&P]
    REAL(SP), DIMENSION(4,4) ::   phi_s  ! See [eqn(8),S&P] - but for surface 
    ! expansion vbf's.  
    REAL(SP), DIMENSION(3):: tempvec1,tempvec2
    ! temporary vectors
    INTEGER(I4B) count                 ! counter

    IF ( ABS(SQRT(DOT_PRODUCT(normal,normal))-1.0_SP).GE.EPS) THEN
       STOP 'IE: B_MAKE_HIERARCHAL_ANALYTIC called with invalid unit normal.'
    END IF

    Bs = 0.0_SP ! Initialize

    CALL SIMPLEX_COEFFICIENTS(element_num,a,b,c,d)    

    IF (CUBATURE) THEN
       ! Calculate GRADIENT_LAMBDA for repeated later use:
       grad_lambda = GRADIENT_LAMBDA(element_num,.TRUE.) 
    END IF

    ! Find relevant SURFACE expansion terms. Note that 
    ! \vec{S}_i =  \hat{n} \times \vec{N}_i

    DO inode = 1,4 ! Node
       tempvec1 = (/ b(inode),c(inode),d(inode) /)
       tempvec2 = CROSS_PRODUCT(normal,tempvec1)
       b_s(inode) = tempvec2(1) ! x component
       c_s(inode) = tempvec2(2) ! y component
       d_s(inode) = tempvec2(3) ! z component
    END DO

    ! Compute dot product terms for SURFACE vector basis functions:  
    ! Extension of [eqns(7)&(8),S&P] 
    DO inode = 1,4 ! Node
       DO jnode = 1,4 ! Node
          phi_s(inode,jnode) = b_s(inode)*b_s(jnode) + c_s(inode)*c_s(jnode) + & 
               d_s(inode)*d_s(jnode)
       END DO
    END DO

    CALL INTEGRATION_MATRIX_2D_ORDER1

    ! Find the lengths of the element edges
    temp_face_nodes = GLOBAL_FACENODES(element_num,local_face_num)
    ell(1) = T_LENGTH(temp_face_nodes(1),temp_face_nodes(2))
    ell(2) = T_LENGTH(temp_face_nodes(1),temp_face_nodes(3))
    ell(3) = T_LENGTH(temp_face_nodes(2),temp_face_nodes(3))

    CALL B_EDGE_EDGE_ORDER1_ANALYTIC

    IF(element_order.EQ.2.AND.element_mixed_order) THEN   
       CALL INTEGRATION_MATRICES_2D_ORDER2
       CALL B_EDGE_EDGE_ORDER2_ANALYTIC
       CALL B_EDGE_FACE_ORDER2_ANALYTIC
       CALL B_FACE_FACE_ORDER2_ANALYTIC
    ELSE
       STOP 'IE: routine B_MAKE_HIERARCHAL_ANALYTIC only valid up to LT/QN elements'
    END IF

    ! CALL B_MAKE_DEBUGGING ! Only when debugging

    ! Geometrical scaling.
    Bs = Bs * FACE_AREA(element_num,local_face_num) 

  CONTAINS

    SUBROUTINE INTEGRATION_MATRIX_2D_ORDER1
      !*******************************************************************************
      ! Sets up integration matrices M  ([eqn. (15) S&P], with common denominator 
      ! modified for surface integral.
      !*******************************************************************************
      M = SPREAD(SPREAD(1.0_SP,1,4),2,4) ! Create a 4x4 matrix of ones
      M(1,1) = 2.0_SP ! Change diagonals
      M(2,2) = 2.0_SP
      M(3,3) = 2.0_SP
      M(4,4) = 2.0_SP
      M = M / 12.0_SP ! Normalize
    END SUBROUTINE INTEGRATION_MATRIX_2D_ORDER1


    SUBROUTINE INTEGRATION_MATRICES_2D_ORDER2
      !*******************************************************************************
      ! Sets up integration matrices N and P , with common denominator 
      ! modified for surface integral.
      !*******************************************************************************
      ! Note that N is such that the 2nd and 3rd nodes are different by definition.

      DO inode  = 1,4
         DO jnode = 1,4
            DO knode = 1,4
               IF (inode.EQ.jnode.OR.inode.EQ.knode) THEN ! One coincident node
                  N(inode,jnode,knode) = 2.0_SP
               ELSE ! Disjoint nodes.
                  N(inode,jnode,knode) = 1.0_SP
               END IF
            END DO
         END DO
      END DO
      N = N/60.0_SP    ! Scaling factor for surface integration

      ! Set up P matrix
      DO inode  = 1,4
         DO jnode = 1,4
            DO knode = 1,4
               DO lnode = 1,4
                  IF (inode.EQ.knode.AND.jnode.EQ.lnode) THEN ! Coincident edge       
                     P(inode,jnode,knode,lnode) = 4.0_SP
                  ELSE IF (inode.EQ.knode) THEN ! Coincident start nodes
                     P(inode,jnode,knode,lnode) = 2.0_SP
                  ELSE IF (jnode.EQ.lnode) THEN ! Coincident end nodes
                     P(inode,jnode,knode,lnode) = 2.0_SP
                  ELSE IF (jnode.EQ.knode) THEN ! Coincident end and start nodes
                     P(inode,jnode,knode,lnode) = 2.0_SP
                  ELSE IF (inode.EQ.lnode) THEN ! Coincident start and end nodes
                     P(inode,jnode,knode,lnode) = 2.0_SP
                  ELSE ! Completely disjoint nodes. This cannot occur here.
                     P(inode,jnode,knode,lnode) = 1.0_SP
                  END IF
               END DO
            END DO
         END DO
      END DO
      P = P/360.0_SP    ! Scaling factor
    END SUBROUTINE INTEGRATION_MATRICES_2D_ORDER2


    SUBROUTINE B_EDGE_EDGE_ORDER1_ANALYTIC
      IMPLICIT NONE
      !*******************************************************************************
      ! Fills edge-edge B submatrices for the lowest order elements and 
      ! stores in elemental B matrix.
      !*******************************************************************************
      COMPLEX(SPC), DIMENSION(3,3) :: Be1e1
      INTEGER(I4B) :: iedge,jedge,inode  ! Misc counters
      INTEGER(I4B) :: i1,i2,j1,j2
      INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
      INTEGER(I4B) quad_point            ! quadature/cubature point index.
      REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.
      INTEGER(I4B) :: num_qpoints        ! number of quadrature points to be used
      INTEGER(I4B), DIMENSION(3) :: tempfacenodes ! local nodes of the face

      ! Initialize (needed for cubature)
      Be1e1 = 0.0_SP
      num_qpoints = 3

      EDGELOOP1: DO iedge = 1,3
         EDGELOOP2: DO jedge = 1,iedge ! Fill lower half
            i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
            i1 = i_edge_nodes(1)
            i2 = i_edge_nodes(2)
            j_edge_nodes = LOCAL_TRI_EDGENODES(jedge,local_face_num)
            j1 = j_edge_nodes(1)
            j2 = j_edge_nodes(2)
            IF (.NOT.CUBATURE) THEN
               Be1e1(iedge,jedge) = ( phi_s(i2,j2)*M(i1,j1) - phi_s(i2,j1)*M(i1,j2)  & 
                    - phi_s(i1,j2)*M(i2,j1) + phi_s(i1,j1)*M(i2,j2) )     
            ELSE 
               ! Integrate numerically, using three-point 2nd order rule [Cowper]
               ! Note that the VBF is \hat{n} \times VBF. 
               tempfacenodes = LOCAL_FACENODES(local_face_num)
               DO quad_point = 1,num_qpoints

                  ! Associate the correct 3D simplex coordinates with the 2D integration:
                  lambda = 0.0
                  DO inode = 1,3
                     lambda(tempfacenodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
                  END DO

                  Be1e1(iedge,jedge) = Be1e1(iedge,jedge) +               &
                       quad_tri_rules(num_qpoints)%rule(quad_point,4) *      & 
                       DOT_PRODUCT(VBF_S(1,lambda,grad_lambda,normal,i1,i2), & 
                       VBF_S(1,lambda,grad_lambda,normal,j1,j2))
               END DO
            END IF
            ! Scale by edge length
            IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
               Be1e1(iedge,jedge) = ell(iedge)*ell(jedge)*Be1e1(iedge,jedge)
            END IF
            ! Symmetrize. Diagonal is overwritten  
            Be1e1(jedge,iedge) = Be1e1(iedge,jedge) 
         END DO EDGELOOP2
      END DO EDGELOOP1
      Bs(1:3,1:3) = Be1e1(1:3,1:3) 
    END SUBROUTINE B_EDGE_EDGE_ORDER1_ANALYTIC

    SUBROUTINE B_EDGE_EDGE_ORDER2_ANALYTIC
      IMPLICIT NONE
      !*******************************************************************************
      ! Fills edge-edge B submatrices for the lowest order elements and 
      ! stores in elemental B matrix. Those computed are Be1e2 and Be2e2; 
      ! Be2e1 is filled using symmetry.
      !*******************************************************************************
      COMPLEX(SPC), DIMENSION(3,3) :: Be1e2, Be2e1, Be2e2
      INTEGER(I4B) :: iedge,jedge,inode ! Misc counters
      INTEGER(I4B) :: i1,i2,j1,j2
      INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
      INTEGER(I4B) quad_point            ! quadature/cubature point index.
      REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.
      INTEGER(I4B) :: num_qpoints        ! number of quadrature points to be used
      INTEGER(I4B), DIMENSION(3) :: tempfacenodes ! local nodes of the face

      ! Initialize (needed for cubature)
      Be1e2 = 0.0_SP
      Be2e2 = 0.0_SP
      num_qpoints = 6

      EDGELOOP1: DO iedge = 1,3
         EDGELOOP2: DO jedge = 1,3
            i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
            i1 = i_edge_nodes(1)
            i2 = i_edge_nodes(2)
            j_edge_nodes = LOCAL_TRI_EDGENODES(jedge,local_face_num)
            j1 = j_edge_nodes(1)
            j2 = j_edge_nodes(2)
            IF (.NOT.CUBATURE) THEN
               Be1e2(iedge,jedge) = ( phi_s(i2,j2)*M(i1,j1) + phi_s(i2,j1)*M(i1,j2)  & 
                    - phi_s(i1,j2)*M(i2,j1) - phi_s(i1,j1)*M(i2,j2) )     
               Be2e2(iedge,jedge) = ( phi_s(i2,j2)*M(i1,j1) + phi_s(i2,j1)*M(i1,j2)  & 
                    + phi_s(i1,j2)*M(i2,j1) + phi_s(i1,j1)*M(i2,j2) )      
            ELSE 
               ! Integrate numerically, using six-point 4-th order rule [Cowper].
               ! (Some edge-based VBF's are 2nd order).
               ! Note that the VBF is \hat{n} \times VBF. 
               tempfacenodes = LOCAL_FACENODES(local_face_num)
               DO quad_point = 1,num_qpoints

                  ! Associate the correct 3D simplex coordinates with the 2D integration:
                  lambda = 0.0
                  DO inode = 1,3
                     lambda(tempfacenodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
                  END DO

                  ! Add the contribution of this integration point - with correct weight:
                  Be1e2(iedge,jedge) = Be1e2(iedge,jedge) +               &
                       quad_tri_rules(num_qpoints)%rule(quad_point,4) *      & 
                       DOT_PRODUCT(VBF_S(1,lambda,grad_lambda,normal,i1,i2), & 
                       VBF_S(2,lambda,grad_lambda,normal,j1,j2))
                  Be2e2(iedge,jedge) = Be2e2(iedge,jedge) +               &
                       quad_tri_rules(num_qpoints)%rule(quad_point,4) *      & 
                       DOT_PRODUCT(VBF_S(2,lambda,grad_lambda,normal,i1,i2), & 
                       VBF_S(2,lambda,grad_lambda,normal,j1,j2))
               END DO
            END IF

            IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
               Be1e2(iedge,jedge) = ell(iedge)*ell(jedge)*Be1e2(iedge,jedge)
               Be2e2(iedge,jedge) = ell(iedge)*ell(jedge)*Be2e2(iedge,jedge)      
            ENDIF

         END DO EDGELOOP2
      END DO EDGELOOP1
      Be2e1 = TRANSPOSE(Be1e2) 
      Bs(4:6,1:3)  = Be2e1(1:3,1:3)              
      Bs(1:3,4:6)  = Be1e2(1:3,1:3)       
      Bs(4:6,4:6) = Be2e2(1:3,1:3)      
    END SUBROUTINE B_EDGE_EDGE_ORDER2_ANALYTIC

    SUBROUTINE B_EDGE_FACE_ORDER2_ANALYTIC
      IMPLICIT NONE
      !*******************************************************************************
      ! Fills edge-face B submatrices for the lowest order elements and 
      ! stores in elemental B matrix. Those computed are Be1f1,Be2f1,Be1f2,Be2f2; 
      ! the others are filled using symmetry.
      !*******************************************************************************
      COMPLEX(SPC), DIMENSION(3,1) :: Be1f1, Be1f2, Be2f1, Be2f2
      COMPLEX(SPC), DIMENSION(1,3) :: Bf1e1, Bf1e2, Bf2e1, Bf2e2
      INTEGER(I4B) :: iedge,inode        ! Misc counters
      INTEGER(I4B) :: i1,i2,j1,j2,j3
      INTEGER(I4B), DIMENSION(2) :: i_edge_nodes
      INTEGER(I4B), DIMENSION(3) :: j_face_nodes
      INTEGER(I4B) quad_point            ! quadature/cubature point index.
      REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.
      INTEGER(I4B) :: num_qpoints        ! number of quadrature points to be used
      INTEGER(I4B), DIMENSION(3) :: tempfacenodes ! local nodes of the face

      ! Initialize (needed for cubature)
      Be1f1 = 0.0_SP
      Be2f1 = 0.0_SP
      Be1f2 = 0.0_SP
      Be2f2 = 0.0_SP
      num_qpoints = 6

      EDGELOOP1: DO iedge = 1,3
         i_edge_nodes = LOCAL_TRI_EDGENODES(iedge,local_face_num)
         i1 = i_edge_nodes(1)
         i2 = i_edge_nodes(2)
         j_face_nodes = LOCAL_FACENODES(local_face_num)
         j1 = j_face_nodes(1)
         j2 = j_face_nodes(2)
         j3 = j_face_nodes(3)
         IF (.NOT.CUBATURE) THEN
            Be1f1(iedge,1) =                                                   &
                 ( phi_s(i2,j3)*N(i1,j1,j2) - phi_s(i1,j3)*N(i2,j1,j2) & 
                 - phi_s(i2,j2)*N(i1,j1,j3) + phi_s(i1,j2)*N(i2,j1,j3) ) 
            Be2f1(iedge,1) =                                                   &
                 ( phi_s(i2,j3)*N(i1,j1,j2) + phi_s(i1,j3)*N(i2,j1,j2) & 
                 - phi_s(i2,j2)*N(i1,j1,j3) - phi_s(i1,j2)*N(i2,j1,j3) )   
            Be1f2(iedge,1) =                                                   &
                 ( phi_s(i2,j3)*N(i1,j1,j2) - phi_s(i1,j3)*N(i2,j1,j2) & 
                 - phi_s(i2,j1)*N(i1,j2,j3) + phi_s(i1,j1)*N(i2,j2,j3) ) 
            Be2f2(iedge,1) =                                                   &
                 ( phi_s(i2,j3)*N(i1,j1,j2) + phi_s(i1,j3)*N(i2,j1,j2) & 
                 - phi_s(i2,j1)*N(i1,j2,j3) - phi_s(i1,j1)*N(i2,j2,j3) ) 
         ELSE 
            ! Integrate numerically, using six-point 4-th order rule [Cowper].
            ! (Some edge-based VBF's are 2nd order).
            ! Note that the VBF is \hat{n} \times VBF. 
            tempfacenodes = LOCAL_FACENODES(local_face_num)
            DO quad_point = 1,num_qpoints
               ! Associate the correct 3D simplex coordinates with the 2D integration:
               lambda = 0.0
               DO inode = 1,3
                  lambda(tempfacenodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
               END DO

               ! Add the contribution of this quadrature point:
               Be1f1(iedge,1) = Be1f1(iedge,1) +                       &
                    quad_tri_rules(num_qpoints)%rule(quad_point,4) *      & 
                    DOT_PRODUCT(VBF_S(1,lambda,grad_lambda,normal,i1,i2), & 
                    VBF_S(3,lambda,grad_lambda,normal,j1,j2,j3))
               Be2f1(iedge,1) = Be2f1(iedge,1) +                       &
                    quad_tri_rules(num_qpoints)%rule(quad_point,4) *      & 
                    DOT_PRODUCT(VBF_S(2,lambda,grad_lambda,normal,i1,i2), & 
                    VBF_S(3,lambda,grad_lambda,normal,j1,j2,j3))
               Be1f2(iedge,1) = Be1f2(iedge,1) +                       &
                    quad_tri_rules(num_qpoints)%rule(quad_point,4) *      &
                    DOT_PRODUCT(VBF_S(1,lambda,grad_lambda,normal,i1,i2), & 
                    VBF_S(4,lambda,grad_lambda,normal,j1,j2,j3))
               Be2f2(iedge,1) = Be2f2(iedge,1) +                       &
                    quad_tri_rules(num_qpoints)%rule(quad_point,4) *      &
                    DOT_PRODUCT(VBF_S(2,lambda,grad_lambda,normal,i1,i2), & 
                    VBF_S(4,lambda,grad_lambda,normal,j1,j2,j3))

            END DO
         END IF

         IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
            Be1f1(iedge,1) = ell(iedge)*Be1f1(iedge,1)
            Be2f1(iedge,1) = ell(iedge)*Be2f1(iedge,1)
            Be1f2(iedge,1) = ell(iedge)*Be1f2(iedge,1)
            Be2f2(iedge,1) = ell(iedge)*Be2f2(iedge,1)
         ENDIF

      END DO EDGELOOP1
      Bf1e1 = TRANSPOSE(Be1f1)
      Bf1e2 = TRANSPOSE(Be2f1)
      Bf2e1 = TRANSPOSE(Be1f2)              
      Bf2e2 = TRANSPOSE(Be2f2)
      Bs(1:3,7)  = Be1f1(1:3,1)
      Bs(1:3,8)  = Be1f2(1:3,1)
      Bs(4:6,7) = Be2f1(1:3,1)
      Bs(4:6,8) = Be2f2(1:3,1)
      Bs(7,1:3)  = Bf1e1(1,1:3)
      Bs(8,1:3)  = Bf2e1(1,1:3)
      Bs(7,4:6) = Bf1e2(1,1:3)
      Bs(8,4:6) = Bf2e2(1,1:3)
    END SUBROUTINE B_EDGE_FACE_ORDER2_ANALYTIC

    SUBROUTINE B_FACE_FACE_ORDER2_ANALYTIC
      IMPLICIT NONE
      !*******************************************************************************
      ! Fills edge-face B submatrices for the lowest order elements and 
      ! stores in elemental B matrix. Those computed are Bf1f1,Bf1f2,Bf1f2; 
      ! Bf2f1 is filled using symmetry.
      !*******************************************************************************
      COMPLEX(SPC), DIMENSION(1,1) :: Bf1f1, Bf1f2, Bf2f1, Bf2f2
      INTEGER(I4B) :: i1,i2,i3,j1,j2,j3,inode
      INTEGER(I4B), DIMENSION(3) :: i_face_nodes,j_face_nodes
      INTEGER(I4B) quad_point            ! quadature/cubature point index.
      REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.
      INTEGER(I4B) :: num_qpoints        ! number of quadrature points to be used
      INTEGER(I4B), DIMENSION(3) :: tempfacenodes ! local nodes of the face

      ! Initialize (needed for cubature)
      Bf1f1 = 0.0_SP
      Bf1f2 = 0.0_SP  
      Bf2f2 = 0.0_SP
      num_qpoints = 6

      ! Fill face-face sub-matrices. 
      ! Note that the index quad (i,j,k,l) always involves at least one
      ! repeated index.
      i_face_nodes = LOCAL_FACENODES(local_face_num)
      i1 = i_face_nodes(1)
      i2 = i_face_nodes(2)
      i3 = i_face_nodes(3)
      j_face_nodes = LOCAL_FACENODES(local_face_num)
      j1 = i_face_nodes(1)
      j2 = i_face_nodes(2)
      j3 = i_face_nodes(3)
      IF (.NOT.CUBATURE) THEN
         Bf1f1(1,1) =                                                      &
              ( phi_s(i3,j3)*P(i1,i2,j1,j2) - phi_s(i2,j3)*P(i1,i3,j1,j2) & 
              - phi_s(i3,j2)*P(i1,i2,j1,j3) + phi_s(i2,j2)*P(i1,i3,j1,j3) ) 
         Bf1f2(1,1) =                                                      &
              ( phi_s(i3,j3)*P(i1,i2,j1,j2) - phi_s(i2,j3)*P(i1,i3,j1,j2) & 
              - phi_s(i3,j1)*P(i1,i2,j2,j3) + phi_s(i2,j1)*P(i1,i3,j2,j3) ) 
         Bf2f2(1,1) =                                                      &
              ( phi_s(i3,j3)*P(i1,i2,j1,j2) - phi_s(i1,j3)*P(i2,i3,j1,j2) & 
              - phi_s(i3,j1)*P(i1,i2,j2,j3) + phi_s(i1,j1)*P(i2,i3,j2,j3) )
      ELSE 
         ! Integrate numerically, using six-point 4-th order rule [Cowper].
         ! Note that the VBF is \hat{n} \times VBF. 
         tempfacenodes = LOCAL_FACENODES(local_face_num)
         DO quad_point = 1,num_qpoints

            ! Associate the correct 3D simplex coordinates with the 2D integration:
            lambda = 0.0
            DO inode = 1,3
               lambda(tempfacenodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
            END DO

            ! Add this quadrature point's contribution:   
            Bf1f1(1,1) = Bf1f1(1,1) +                                  &
                 quad_tri_rules(num_qpoints)%rule(quad_point,4) *         & 
                 DOT_PRODUCT(VBF_S(3,lambda,grad_lambda,normal,i1,i2,i3), & 
                 VBF_S(3,lambda,grad_lambda,normal,j1,j2,j3))
            Bf1f2(1,1) = Bf1f2(1,1) +                                  &
                 quad_tri_rules(num_qpoints)%rule(quad_point,4) *         & 
                 DOT_PRODUCT(VBF_S(3,lambda,grad_lambda,normal,i1,i2,i3), & 
                 VBF_S(4,lambda,grad_lambda,normal,j1,j2,j3))
            Bf2f2(1,1) = Bf2f2(1,1) +                                  &
                 quad_tri_rules(num_qpoints)%rule(quad_point,4) *         & 
                 DOT_PRODUCT(VBF_S(4,lambda,grad_lambda,normal,i1,i2,i3), & 
                 VBF_S(4,lambda,grad_lambda,normal,j1,j2,j3))
         END DO
      END IF

      Bf2f1 = TRANSPOSE(Bf1f2)
      Bs(7,7)  = Bf1f1(1,1)     
      Bs(7,8)  = Bf1f2(1,1)          
      Bs(8,7)  = Bf2f1(1,1)
      Bs(8,8)  = Bf2f2(1,1)     

    END SUBROUTINE B_FACE_FACE_ORDER2_ANALYTIC


    SUBROUTINE B_MAKE_DEBUGGING
      !*******************************************************************************
      ! Temporary debugging code. 
      !*******************************************************************************
      ! Testing only:
      REAL(SP), DIMENSION(4,3) :: grad_lambda
      REAL(SP), DIMENSION(3):: grad_lambda1,grad_lambda2,grad_lambda3
      REAL(SP) temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9

      ! Tests
      ! Only for one case, Xempty0, dof 7 (f1 type)
      grad_lambda = GRADIENT_LAMBDA(element_num,.true.)
      grad_lambda1 = grad_lambda(1,1:3)
      grad_lambda2 = grad_lambda(3,1:3)
      grad_lambda3 = grad_lambda(4,1:3)

      ! Bf1f1:
      grad_lambda1 = CROSS_PRODUCT(normal,grad_lambda1)
      grad_lambda2 = CROSS_PRODUCT(normal,grad_lambda2)
      grad_lambda3 = CROSS_PRODUCT(normal,grad_lambda3)

      temp1 = (4.0*dot_product(grad_lambda3,grad_lambda3) & 
           -2.0*dot_product(grad_lambda2,grad_lambda3) & 
           -2.0*dot_product(grad_lambda3,grad_lambda2) & 
           +4.0*dot_product(grad_lambda2,grad_lambda2))/360.0

      ! Bf1f2: 
      temp2 = (4.0*dot_product(grad_lambda3,grad_lambda3) & 
           -2.0*dot_product(grad_lambda2,grad_lambda3) & 
           -2.0*dot_product(grad_lambda3,grad_lambda1) & 
           +2.0*dot_product(grad_lambda2,grad_lambda1))/360.0


      ! Bf2f2: 
      temp3 = (4.0*dot_product(grad_lambda3,grad_lambda3) & 
           -2.0*dot_product(grad_lambda3,grad_lambda1) & 
           -2.0*dot_product(grad_lambda1,grad_lambda3) & 
           +4.0*dot_product(grad_lambda1,grad_lambda1))/360.0


      ! Be1f1(1,1)
      temp4  = (dot_product(grad_lambda2,grad_lambda3) & 
           - dot_product(grad_lambda1,grad_lambda3) & 
           - dot_product(grad_lambda2,grad_lambda2) & 
           + dot_product(grad_lambda1,grad_lambda2) ) /30.0
      ! Be2f1(1,1)
      temp5  = (dot_product(grad_lambda2,grad_lambda3) & 
           + dot_product(grad_lambda1,grad_lambda3) & 
           - dot_product(grad_lambda2,grad_lambda2) & 
           - dot_product(grad_lambda1,grad_lambda2) ) /30.0
      ! Be2f1(1,1)
      temp6  = (2*dot_product(grad_lambda2,grad_lambda3) & 
           - 2*dot_product(grad_lambda1,grad_lambda3) & 
           -   dot_product(grad_lambda2,grad_lambda1) & 
           + 2*dot_product(grad_lambda1,grad_lambda1) ) /60.0
      ! Be2f2(1,1)
      temp7  = (2*dot_product(grad_lambda2,grad_lambda3) & 
           + 2*dot_product(grad_lambda1,grad_lambda3) & 
           -   dot_product(grad_lambda2,grad_lambda1) & 
           - 2*dot_product(grad_lambda1,grad_lambda1) ) /60.0


      ! Be1e2(1,1)
      temp8  = (2*dot_product(grad_lambda2,grad_lambda2) & 
           -  dot_product(grad_lambda1,grad_lambda2) & 
           +  dot_product(grad_lambda2,grad_lambda1) & 
           -2*dot_product(grad_lambda1,grad_lambda1) ) /12.0

      ! Be2e2(1,1)
      temp9  = (2*dot_product(grad_lambda2,grad_lambda2) & 
           +  dot_product(grad_lambda1,grad_lambda2) & 
           +  dot_product(grad_lambda2,grad_lambda1) & 
           +2*dot_product(grad_lambda1,grad_lambda1) ) /12.0

      !if (element_num.EQ.1) THEN
      !print *,P(1,2,1,2)

      !print *,'n X nabla lambda'
      !print *,'Code:'
      !do inode=1,4
      !print *,b_s(inode),c_s(inode),d_s(inode)
      !end do

      !print *,'local nodes: ',i_face_nodes,j_face_nodes

      !print *,'Check:'
      !print *,grad_lambda1
      !print *,grad_lambda2
      !print *,grad_lambda3

      !print *,'Entry Bs(1,4) computed as per S&P:',Bs(1,4)
      !print *,'Entry Bs(1,4) checked:',temp8
      !print *,'Entry Bs(4,4) computed as per S&P:',Bs(4,4)
      !print *,'Entry Bs(4,4) checked:',temp9
      !print *,'Entry Bs(4,7) computed as per S&P:',Bs(4,7)
      !print *,'Entry Bs(4,7) checked:',temp5
      !print *,'Entry Bs(1,8) computed as per S&P:',Bs(1,8)
      !print *,'Entry Bs(1,8) checked:',temp6
      !print *,'Entry Bs(4,8) computed as per S&P:',Bs(4,8)
      !print *,'Entry Bs(4,8) checked:',temp7
      !print *,'Entry Bs(7,7) computed as per S&P:',Bs(7,7)
      !print *,'Entry Bs(7,7) checked:',temp1
      !print *,'Entry Bs(7,8) computed as per S&P:',Bs(7,8)
      !print *,'Entry Bs(7,8) checked:',temp2
      !print *,'Entry Bs(8,8) computed as per S&P:',Bs(8,8)
      !print *,'Entry Bs(8,8) checked:',temp3
      !end if 
    END SUBROUTINE B_MAKE_DEBUGGING
    !*******************************************************************************

  END SUBROUTINE B_MAKE_HIERARCHAL_ANALYTIC
  !*******************************************************************************

END MODULE B_matrix

!**********************************************************************
!* 
!* MODULE S_intg_funs
!* 
!* Input parameters: (module variables that need to be set by the user)
!*   
!*   V, iedge, jedge, i1, i2, j1, j2 (Up to LT/LN),  and
!*   iface, jface, j3                (Up to LT/QN) and
!*   i3                              (Up to QT/TN)
!*   need to be set per per element. (see S_EDGE_EDGE/S_EDGE_FACE
!*                                    and S_FACE_FACE)
!*   elinfo also needs to be set per element when curvilinear elements 
!*   used.
!*   
!* Output:
!*   
!*   The product of the curls of various basis function types. eg. 
!*   curl_e1e1 returns the dot product of the curl of two e1 type
!*   edge functions, as described by V, i1, i2, j1 and j2. 
!*   
!* Description:
!*   
!*   e1, e2, e3, f1, f2 and f3 refer to the edge and face based.              
!*   functions See [S&P] and [D&H].  e3 and f3 are  defined in 
!*   [Webb99], using the notation G_2^(e) and G_20^(f) respectively.
!*   These functions are in the gradient space (of order two), and 
!*   have zero curl.
!*   
!*   
!**********************************************************************
MODULE S_intg_funs
  USE nrtype
  USE geometry
  USE basis_function
  IMPLICIT NONE
  REAL(SP), DIMENSION(4,4,3) ::   V        ! See [eqn(7),S&P]
  TYPE(element_info_type) :: elinfo
  INTEGER(I4B) :: iedge,jedge              ! edge counters
  INTEGER(I4B) :: iface,jface              ! face counters
  INTEGER(I4B) :: i1,i2,i3,j1,j2,j3        ! node counters
  REAL(SP), DIMENSION(3,3)  :: sigma_tensor

CONTAINS

  REAL(SP) FUNCTION curl_e1e1(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    curl_e1e1 = DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(1,lambda,V,j1,j2))
  END FUNCTION curl_e1e1

  REAL(SP) FUNCTION curl_e1e1_curvi(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    ! Local coordinates on parent unitary tet for curvilinear elements.
    REAL(SP) :: uu,vv,ww   
    REAL(SP), DIMENSION(3,3) :: temp_mat 
    REAL(SP), DIMENSION(3) :: temp_vec
    uu=lambda(2)
    vv=lambda(3)
    ww=lambda(4)
    temp_mat= MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.),&
         TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    temp_vec= MATMUL(temp_mat,CURL_PARENT_VBF(1,jedge))
    curl_e1e1_curvi = DOT_PRODUCT(CURL_PARENT_VBF(1,iedge),temp_vec)&
         /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    ! The absolute value of the determinant is required since the parent 
    ! element is always "right hand" numbered, whereas the actual element 
    ! may be left or right hand numbered. 
  END FUNCTION curl_e1e1_curvi

  REAL(SP) FUNCTION curl_e1e1_tensor(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    curl_e1e1_tensor = DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2), &
         MATMUL(sigma_tensor,CURL_VBF(1,lambda,V,j1,j2)))
  END FUNCTION curl_e1e1_tensor

  REAL(SP) FUNCTION curl_e1e2(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    curl_e1e2 = DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(2,lambda,V,j1,j2))
  END FUNCTION curl_e1e2

  REAL(SP) FUNCTION curl_e2e2(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    curl_e2e2 = DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(2,lambda,V,j1,j2))
  END FUNCTION curl_e2e2

  REAL(sp) FUNCTION curl_e1f1(lambda)
    REAL(sp), DIMENSION(4), INTENT(in) :: lambda
    curl_e1f1 = DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(3,lambda,V,j1,j2,j3))
  END FUNCTION curl_e1f1

  REAL(sp) FUNCTION curl_e1f2(lambda)
    REAL(sp), DIMENSION(4), INTENT(in) :: lambda
    curl_e1f2 = DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(4,lambda,V,j1,j2,j3))
  END FUNCTION curl_e1f2

  REAL(SP) FUNCTION curl_e1f1_curvi(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    ! Local coordinates on parent unitary tet for curvilinear elements.
    REAL(SP) :: uu,vv,ww   
    REAL(SP), DIMENSION(3,3) :: temp_mat 
    REAL(SP), DIMENSION(3) :: temp_vec
    uu=lambda(2)
    vv=lambda(3)
    ww=lambda(4)
    temp_mat = MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
         TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    temp_vec = MATMUL(temp_mat,CURL_PARENT_VBF(3,jface,uu,vv,ww))
    curl_e1f1_curvi = DOT_PRODUCT(CURL_PARENT_VBF(1,iedge),temp_vec) &
         /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    ! The absolute value of the determinant is required since the parent 
    ! element is always "right hand" numbered, whereas the actual element 
    ! may be left or right hand numbered. 
  END FUNCTION curl_e1f1_curvi

  REAL(SP) FUNCTION curl_e1f2_curvi(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    ! Local coordinates on parent unitary tet for curvilinear elements.
    REAL(SP) :: uu,vv,ww   
    REAL(SP), DIMENSION(3,3) :: temp_mat 
    REAL(SP), DIMENSION(3) :: temp_vec
    uu=lambda(2)
    vv=lambda(3)
    ww=lambda(4)
    temp_mat = MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
         TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    temp_vec = MATMUL(temp_mat,CURL_PARENT_VBF(4,jface,uu,vv,ww))
    curl_e1f2_curvi = DOT_PRODUCT(CURL_PARENT_VBF(1,iedge),temp_vec) &
         /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    ! The absolute value of the determinant is required since the parent 
    ! element is always "right hand" numbered, whereas the actual element 
    ! may be left or right hand numbered. 
  END FUNCTION curl_e1f2_curvi

  REAL(sp) FUNCTION curl_e2f1(lambda)
    REAL(sp), DIMENSION(4), INTENT(in) :: lambda
    curl_e2f1 = DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(3,lambda,V,j1,j2,j3))
  END FUNCTION curl_e2f1

  REAL(sp) FUNCTION curl_e2f2(lambda)
    REAL(sp), DIMENSION(4), INTENT(in) :: lambda
    curl_e2f2 = DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(4,lambda,V,j1,j2,j3))
  END FUNCTION curl_e2f2

  REAL(SP) FUNCTION curl_e2f1_curvi(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    ! Local coordinates on parent unitary tet for curvilinear elements.
    REAL(SP) :: uu,vv,ww   
    REAL(SP), DIMENSION(3,3) :: temp_mat 
    REAL(SP), DIMENSION(3) :: temp_vec
    uu=lambda(2)
    vv=lambda(3)
    ww=lambda(4)
    temp_mat= MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
         TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    temp_vec= MATMUL(temp_mat,CURL_PARENT_VBF(3,jface,uu,vv,ww))
    curl_e2f1_curvi = DOT_PRODUCT(CURL_PARENT_VBF(2,iedge),temp_vec) &
         /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    ! The absolute value of the determinant is required since the parent 
    ! element is always "right hand" numbered, whereas the actual element 
    ! may be left or right hand numbered. 
  END FUNCTION curl_e2f1_curvi

  REAL(SP) FUNCTION curl_e2f2_curvi(lambda)
    REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
    ! Local coordinates on parent unitary tet for curvilinear elements.
    REAL(SP) :: uu,vv,ww   
    REAL(SP), DIMENSION(3,3) :: temp_mat 
    REAL(SP), DIMENSION(3) :: temp_vec
    uu=lambda(2)
    vv=lambda(3)
    ww=lambda(4)
    temp_mat= MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
         TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    temp_vec= MATMUL(temp_mat,CURL_PARENT_VBF(4,jface,uu,vv,ww))
    curl_e2f2_curvi = DOT_PRODUCT(CURL_PARENT_VBF(2,iedge),temp_vec) & 
         /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
    ! The absolute value of the determinant is required since the parent 
    ! element is always "right hand" numbered, whereas the actual element 
    ! may be left or right hand numbered. 
  END FUNCTION curl_e2f2_curvi


    REAL(sp) FUNCTION curl_f1f1(lambda)
      REAL(sp), DIMENSION(4), INTENT(in) :: lambda
      curl_f1f1 = DOT_PRODUCT(CURL_VBF(3,lambda,V,i1,i2,i3),CURL_VBF(3,lambda,V,j1,j2,j3))
    END FUNCTION curl_f1f1

    REAL(sp) FUNCTION curl_f1f2(lambda)
      REAL(sp), DIMENSION(4), INTENT(in) :: lambda
      curl_f1f2 = DOT_PRODUCT(CURL_VBF(3,lambda,V,i1,i2,i3),CURL_VBF(4,lambda,V,j1,j2,j3))
    END FUNCTION curl_f1f2

    REAL(sp) FUNCTION curl_f2f2(lambda)
      REAL(sp), DIMENSION(4), INTENT(in) :: lambda
      curl_f2f2 = DOT_PRODUCT(CURL_VBF(4,lambda,V,i1,i2,i3),CURL_VBF(4,lambda,V,j1,j2,j3))   
    END FUNCTION curl_f2f2

    REAL(SP) FUNCTION curl_f1f1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
           TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      temp_vec= MATMUL(temp_mat,CURL_PARENT_VBF(3,jface,uu,vv,ww))
      curl_f1f1_curvi = DOT_PRODUCT(CURL_PARENT_VBF(3,iface,uu,vv,ww),temp_vec) &
           /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION curl_f1f1_curvi

    REAL(SP) FUNCTION curl_f1f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat = MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
           TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      temp_vec = MATMUL(temp_mat,CURL_PARENT_VBF(4,jface,uu,vv,ww))
      curl_f1f2_curvi = DOT_PRODUCT(CURL_PARENT_VBF(3,iface,uu,vv,ww),temp_vec) &
           /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION curl_f1f2_curvi

    REAL(SP) FUNCTION curl_f2f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat = MATMUL(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.), &
           TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      temp_vec = MATMUL(temp_mat,CURL_PARENT_VBF(4,jface,uu,vv,ww))
      curl_f2f2_curvi = DOT_PRODUCT(CURL_PARENT_VBF(4,iface,uu,vv,ww),temp_vec) &
           /ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION curl_f2f2_curvi

  END MODULE S_intg_funs

  !**********************************************************************
  !* 
  !* MODULE S_intg_funs
  !* 
  !* Input parameters: (module variables that need to be set by the user)
  !*   
  !*   V, iedge, jedge, i1, i2, j1, j2 (Up to LT/LN),  and
  !*   iface, jface, j3                (Up to LT/QN) and
  !*   i3                              (Up to QT/TN)
  !*   need to be set per per element. (see S_EDGE_EDGE/S_EDGE_FACE
  !*                                    and S_FACE_FACE)
  !*   elinfo also needs to be set per element when curvilinear elements 
  !*   used.
  !*   
  !* Output:
  !*   
  !*   The product of the curls of various basis function types. eg. 
  !*   curl_e1e1 returns the dot product of the curl of two e1 type
  !*   edge functions, as described by V, i1, i2, j1 and j2. 
  !*   
  !* Description:
  !*   
  !*   e1, e2, e3, f1, f2 and f3 refer to the edge and face based.              
  !*   functions See [S&P] and [D&H].  e3 and f3 are  defined in 
  !*   [Webb99], using the notation G_2^(e) and G_20^(f) respectively.
  !*   These functions are in the gradient space (of order two), and 
  !*   have zero curl.
  !*   
  !*   
  !**********************************************************************
  MODULE T_intg_funs
    USE nrtype
    USE geometry
    USE basis_function
    IMPLICIT NONE
    REAL(SP), DIMENSION(4,3) :: grad_lambda ! Gradient of simplex coordinates
    TYPE(element_info_type) :: elinfo
    INTEGER(I4B) :: iedge,jedge              ! edge counters
    INTEGER(I4B) :: iface,jface              ! face counters
    INTEGER(I4B) :: i1,i2,i3,j1,j2,j3        ! node counters
    REAL(SP), DIMENSION(3,3)  :: sigma_tensor

  CONTAINS

    REAL(SP) FUNCTION e1e1(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1e1 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(1,lambda,grad_lambda,j1,j2))
    END FUNCTION e1e1

    REAL(SP) FUNCTION e1e1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(1,lambda,grad_lambda,j1,j2))
      e1e1_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1e1_curvi

    REAL(SP) FUNCTION e1e1_tensor(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1e1_tensor = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2), &
           MATMUL(sigma_tensor,VBF(1,lambda,grad_lambda,j1,j2)))
    END FUNCTION e1e1_tensor

    REAL(SP) FUNCTION e1e2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1e2 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(2,lambda,grad_lambda,j1,j2))
    END FUNCTION e1e2

    REAL(SP) FUNCTION e1e2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(2,lambda,grad_lambda,j1,j2))
      e1e2_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1e2_curvi

    REAL(SP) FUNCTION e2e2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e2e2 = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(2,lambda,grad_lambda,j1,j2))
    END FUNCTION e2e2

    REAL(SP) FUNCTION e2e2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(2,lambda,grad_lambda,j1,j2))
      e2e2_curvi = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e2e2_curvi

    REAL(SP) FUNCTION e1e3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1e3 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(5,lambda,grad_lambda,j1,j2))
    END FUNCTION e1e3

    REAL(SP) FUNCTION e2e3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e2e3 = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(5,lambda,grad_lambda,j1,j2))
    END FUNCTION e2e3

    REAL(SP) FUNCTION e3e3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e3e3 = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),VBF(5,lambda,grad_lambda,j1,j2))
    END FUNCTION e3e3

    REAL(SP) FUNCTION e1e3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(5,lambda,grad_lambda,j1,j2))      
      e1e3_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1e3_curvi

    REAL(SP) FUNCTION e2e3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(5,lambda,grad_lambda,j1,j2))      
      e2e3_curvi = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e2e3_curvi

    REAL(SP) FUNCTION e3e3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(5,lambda,grad_lambda,j1,j2))      
      e3e3_curvi = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e3e3_curvi

    REAL(SP) FUNCTION e1f1(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1f1 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(3,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e1f1

    REAL(SP) FUNCTION e2f1(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e2f1 = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(3,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e2f1

    REAL(SP) FUNCTION e1f2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1f2 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(4,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e1f2

    REAL(SP) FUNCTION e2f2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e2f2 = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(4,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e2f2

    REAL(SP) FUNCTION e1f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e1f3 = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e1f3

    REAL(SP) FUNCTION e2f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e2f3 = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e2f3

    REAL(SP) FUNCTION e3f1(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e3f1 = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),VBF(3,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e3f1

    REAL(SP) FUNCTION e3f2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e3f2 = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),VBF(4,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e3f2

    REAL(SP) FUNCTION e3f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      e3f3 = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION e3f3

    REAL(SP) FUNCTION e1f1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(3,lambda,grad_lambda,j1,j2,j3))
      e1f1_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1f1_curvi

    REAL(SP) FUNCTION e2f1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(3,lambda,grad_lambda,j1,j2,j3))
      e2f1_curvi = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e2f1_curvi

    REAL(SP) FUNCTION e1f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(4,lambda,grad_lambda,j1,j2,j3))
      e1f2_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1f2_curvi

    REAL(SP) FUNCTION e2f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(4,lambda,grad_lambda,j1,j2,j3))
      e2f2_curvi = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e2f2_curvi

    REAL(SP) FUNCTION e1f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      e1f3_curvi = DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))

      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e1f3_curvi

    REAL(SP) FUNCTION e2f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      e2f3_curvi = DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e2f3_curvi

    REAL(SP) FUNCTION e3f1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(3,lambda,grad_lambda,j1,j2,j3))      
      e3f1_curvi = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e3f1_curvi

    REAL(SP) FUNCTION e3f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(4,lambda,grad_lambda,j1,j2,j3))
      e3f2_curvi = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e3f2_curvi

    REAL(SP) FUNCTION e3f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      e3f3_curvi = DOT_PRODUCT(VBF(5,lambda,grad_lambda,i1,i2),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION e3f3_curvi

    REAL(SP) FUNCTION f1f1(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f1f1 = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),VBF(3,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f1f1

    REAL(SP) FUNCTION f1f2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f1f2 = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),VBF(4,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f1f2

    REAL(SP) FUNCTION f2f2(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f2f2 = DOT_PRODUCT(VBF(4,lambda,grad_lambda,i1,i2,i3),VBF(4,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f2f2

    REAL(SP) FUNCTION f1f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f1f3 = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f1f3

    REAL(SP) FUNCTION f2f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f2f3 = DOT_PRODUCT(VBF(4,lambda,grad_lambda,i1,i2,i3),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f2f3

    REAL(SP) FUNCTION f3f3(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      f3f3 = DOT_PRODUCT(VBF(6,lambda,grad_lambda,i1,i2,i3),VBF(6,lambda,grad_lambda,j1,j2,j3))
    END FUNCTION f3f3

    REAL(SP) FUNCTION f1f1_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(3,lambda,grad_lambda,j1,j2,j3))
      f1f1_curvi = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f1f1_curvi

    REAL(SP) FUNCTION f1f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(4,lambda,grad_lambda,j1,j2,j3))
      f1f2_curvi = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f1f2_curvi

    REAL(SP) FUNCTION f2f2_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(4,lambda,grad_lambda,j1,j2,j3))
      f2f2_curvi = DOT_PRODUCT(VBF(4,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f2f2_curvi

    REAL(SP) FUNCTION f1f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      f1f3_curvi = DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f1f3_curvi

    REAL(SP) FUNCTION f2f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      f2f3_curvi = DOT_PRODUCT(VBF(4,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f2f3_curvi

    REAL(SP) FUNCTION f3f3_curvi(lambda)
      REAL(SP), DIMENSION(4), INTENT(in) :: lambda   ! Simplex coordinates of point.
      ! Local coordinates on parent unitary tet for curvilinear elements.
      REAL(SP) :: uu,vv,ww   
      REAL(SP), DIMENSION(3,3) :: temp_mat 
      REAL(SP), DIMENSION(3) :: temp_vec
      uu=lambda(2)
      vv=lambda(3)
      ww=lambda(4)
      temp_mat= MATMUL(TRANSPOSE(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.)), &
           JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.TRUE.))
      temp_vec= MATMUL(temp_mat,VBF(6,lambda,grad_lambda,j1,j2,j3))
      f3f3_curvi = DOT_PRODUCT(VBF(6,lambda,grad_lambda,i1,i2,i3),temp_vec) &
           *ABS(DET_DIM3(JACOBIAN_POLY2(elinfo%num,uu,vv,ww,.FALSE.)))
      ! The absolute value of the determinant is required since the parent 
      ! element is always "right hand" numbered, whereas the actual element 
      ! may be left or right hand numbered. 
    END FUNCTION f3f3_curvi

  END MODULE T_intg_funs

MODULE S_and_T_matrix
  USE nrtype
  USE basis_function
  USE problem_info, ONLY: SCALE_BY_EDGE_LENGTH
  USE geometry
  USE quad_tables
  IMPLICIT NONE



CONTAINS

  SUBROUTINE S_AND_T_MAKE_HIERARCHAL(element_num,element_order,element_mixed_order,Se,Te,& 
       no_material_scaling,sigma_tensor,element_curvilinear)
    USE geometry
    USE problem_info
    USE quad_tables
    USE material_properties
    IMPLICIT NONE
    !*******************************************************************************
    ! SUBROUTINE DESCRIPTION
    !*******************************************************************************
    ! This subroutine constructs the S and T matrices for element i using 
    ! hierarchal elements. Currently implemented are H_0(curl) (CT/LN) and 
    ! H_1(curl) (LT/QN)  elements. Two versions of the latter are available. 
    ! The formulation used is that of:
    ! [S&P] Savage and Peterson, "Higher-order vector finite elements for tetrahedral 
    ! cells", IEEE Trans-MTT, June 96, pp. 874-879, 
    !
    ! with extensions: 
    !
    ! [D&H] DB Davidson, RH Hansmann, "Hierarchal 2D and 3D Vector Finite Elements 
    ! for Electromagnetic Wave Eigenvalue Problems",
    ! ACES 15th Annl. Review of Progress in Applied Computational Electromagnetics, 
    ! March 1999, pp.518--521
    !
    ! and 
    ! [A&V]: LS Andersen,JL Volakis, "Hierarhical tangential vector finite 
    ! elements for tetrahedra",
    ! IEEE M&GW Lttrs, March 1998, p.127--129.
    ! [Webb99]: J Webb, "Hierarchal vector basis functions of abritrary order for
    ! triangular and tetrehedral vector elements", IEEE AP-Trans, Aug 1999, p. 1244-1253.
    !
    ! Elemental contributions  are computed using quadrature. The [S&P] (and [D&H] extensions) 
    ! can optionally also be computed using closed form expressions. 
    !
    ! An option has now been added to include tensor material properties. This is useful
    ! for PML-type materials. 
    !
    !*******************************************************************************
    ! AUTHOR
    !*******************************************************************************
    ! D B Davidson
    !
    !*******************************************************************************
    ! Last revised
    !*******************************************************************************
    ! 24 Jan 1999 DBD. 
    ! 18 Feb 1999 DBD.
    ! 23 Jan 2000 DBD. Error handling streamlined.
    ! 15 Feb 2000 DBD. Error in double use of ipiv and work corrected.
    !  4 Mar 2000 DBD. Updated to use new data structures.
    ! 14 Mar 2000 DBD. References updated.
    ! 17 Mar 2000 DBD. Simplex coordinates computed in external routine.
    ! 18 May 2000 DBD. Minor documentation updates. 
    ! 31 Jan 2001 MMB. Error handling further streamlined.
    ! 12 Feb 2001 DBD. Cubature and A&V elements added. 
    ! 07 Mar 2002 DBD. Polynomial complete elements up to QT/QN order added. 
    ! 25 Mar 2003 DBD. Option added to NOT include eps and mu in S and T. 
    ! 11 Jun 2003 DBD. Some preliminary PML (diagonally anistropic material treatment 
    !                  added for CT/LN. 
    ! 13 Jul 2005 DBD. Support for 2nd order polynomial curvilinear elements added, to LT/QN order.
    !
    !*******************************************************************************
    ! INPUT 
    !*******************************************************************************
    ! Input data is the element number, order and element geometry.
    !
    !*******************************************************************************
    ! OUTPUT 
    !*******************************************************************************
    ! Output data is the elemental S and T matrices. The full system matrix
    ! assembly is done elsewhere, taking the BC's into account. 
    ! 
    ! The matrices are scaled by eps_r and mu_r unless the last (optional) flag is
    ! set. (The option to return unscaled values is primarily intended for time domain  
    ! analysis).
    !
    ! The storage convention used is the following:
    !
    ! | Se1e1 Se1e2 Se1f1 Se1f2 Se1e3 Se1f3|
    ! | Se2e1 Se2e2 Se2f1 Se2f2 Se2e3 Se2f3|
    ! | Sf1e1 Sf1e2 Sf1f1 Sf1f2 Sf1e3 Sf1f3|
    ! | Se1e1 Se1e2 Se1f1 Se1f2 Se1e3 Se1f3|
    ! | Se3e1 Se3e2 Se3f1 Se3f2 Se3e3 Se3f3|
    ! | Sf3e1 Sf3e2 Sf3f1 Sf3f2 Sf3e3 Sf3f3|
    !
    ! where e1, e2, e3, f1, f2 and f3 refer to the edge and face based functions.
    ! See [S&P] and [D&H]. Similarly for T. e3 and f3 are defined in [Webb99], using the 
    ! notation G_2^(e) and G_20^(f) respectively. These functions are in the gradient space (of 
    ! order two), and have zero curl. 
    !
    ! Note the- [S&P] include the length of the edge in the H_0(curl) element;
    ! [L&M] do not. An option is included to omit the length to permit 
    ! direct comparison between the matrices generated by the two. The 
    ! solution is not affected by this but post-processing is. 
    ! 
    ! The numbering convention of nodes, edges and faces follows [S&P] in this
    ! sub-routine [Table II,S&P]. The edge numbering scheme is as for Lee and
    ! Mittra, apart from the offset of one in the latter (i.e. e12 in S&P is
    ! e01 in L&M). 
    !
    ! The polynomial complete elements follow the work of J Webb:
    ! "Hierarchal vector basis functions of arbitrary order for 
    !  triangular and tetrahedral finite elements"
    ! IEEE T-AP, vol 47, No 8, Aug 1999 p. 1244ff. 
    !

    ! 
    !
    !*******************************************************************************   
!!!
!!! Interface Variables
!!!
    INTEGER(I4B), INTENT(IN) :: element_num, element_order
    LOGICAL(LGT), INTENT(IN) :: element_mixed_order
    LOGICAL(LGT), INTENT(IN), OPTIONAL :: no_material_scaling
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor
    LOGICAL(LGT), INTENT(IN), OPTIONAL :: element_curvilinear
!!!
!!! Output Variables
!!!  
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se,Te
!!!
!!! Internal Variables
!!!
    ! Dimensioned up to and including H_1(curl) (LT/QN)elements.
    INTEGER(I4B) :: knode,lnode        !  "
    INTEGER(I4B) count
    REAL(SP), DIMENSION(4,3) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,4) :: M      ! integration matrix M [eqn(15),S&P]
    REAL(SP), DIMENSION(4,4,4) :: N    ! integration matrix N [eqn(46),S&P]
    REAL(SP), DIMENSION(4,4,4,4) :: P  ! integration matrix P [eqn(47),S&P]
    REAL(SP), DIMENSION(6) :: ell      ! edge lengths
    REAL(SP), DIMENSION(4,4) ::   phi  ! See [eqn(8),S&P] 
    REAL(SP), DIMENSION(4,4,3) ::   V  ! See [eqn(7),S&P]
    TYPE(element_info_type) :: elinfo
    LOGICAL(LGT) no_material_scaling_flag
    LOGICAL(LGT) tensor_material_flag
    LOGICAL(LGT) curvilinear_flag
!!!
    tensor_material_flag = PRESENT(sigma_tensor)
    IF(PRESENT(element_curvilinear)) THEN 
       curvilinear_flag = element_curvilinear
    ELSE
       curvilinear_flag = .FALSE.
    END IF

    IF(tensor_material_flag.AND..NOT.CUBATURE.AND.ELEMENT_ORDER.GT.1.AND..NOT.MIXED_ORDER_FLAG.AND..NOT.curvilinear_flag) THEN
       STOP 'IE. PML FETD only implemented for cubature, CT/LN basis functions on rectilinear elements.'
    END IF

    IF (CUBATURE) THEN

       elinfo%num = element_num            ! Init elinfo structure to be passed
       elinfo%order = element_order        ! to lower level routines.
       elinfo%mixed_order = element_mixed_order

       ! Calculate GRADIENT_LAMBDA for repeated later use:
       IF(curvilinear_flag) THEN ! Added DBD 13 July 05
          ! These are computed directly from a parent unitary tet, with \lambda2=x,\lambda3=y & \lambda4=z 
          grad_lambda(1,:) = (/ -1.0_SP,-1.0_SP,-1.0_SP /)
          grad_lambda(2,:) = (/ 1.0_SP,0.0_SP,0.0_SP /)
          grad_lambda(3,:) = (/ 0.0_SP,1.0_SP,0.0_SP /)
          grad_lambda(4,:) = (/ 0.0_SP,0.0_SP,1.0_SP /) 
       ELSE 
          grad_lambda = GRADIENT_LAMBDA(element_num,.TRUE.)
       END IF ! End DBD 13 July 05.
    END IF

    ! Compute dot product terms:  See [eqns(7)&(8),S&P] 
    ! Needed for both closed form and cubature implementations. (Not used by the curvilinear formulation, 
    ! but can still be computed). 
    CALL S_and_T_make_V_phi(V, phi, element_num)
    ! Find the lengths of the element edges
    ell = T_LENGTHS(element_num) ! Needed even if not thus scaled.
    !Initialize
    Se = 0.0_SP 
    Te = 0.0_SP 

    IF (CUBATURE) THEN ! This will become the default in future releases.
       ! Note that these routines return only the upper half
       ! of the elemental matrices.
       CALL S_EDGE_EDGE(Se, elinfo, V, curvilinear_flag, sigma_tensor)
       CALL T_EDGE_EDGE(Te, elinfo, grad_lambda, curvilinear_flag, sigma_tensor)
       IF(element_order.GE.2) THEN 
          CALL S_EDGE_FACE(Se, elinfo, V, curvilinear_flag, sigma_tensor) 
          CALL S_FACE_FACE(Se, elinfo, V, curvilinear_flag, sigma_tensor) 
          CALL T_EDGE_FACE(Te, elinfo, grad_lambda, curvilinear_flag, sigma_tensor)
          CALL T_FACE_FACE(Te, elinfo, grad_lambda, curvilinear_flag, sigma_tensor)
       ENDIF
       CALL FILL_LOWER_ELEMENTAL_MATRICES(Se, ELEM_TET_MATRIX_SIZE)
       CALL FILL_LOWER_ELEMENTAL_MATRICES(Te, ELEM_TET_MATRIX_SIZE)
    ELSE
       ! This is depreciated code, may become redundant in future releases.
       ! These return the full elemental matrices.
       IF (element_mixed_order) THEN
          M = INTEGRATION_MATRIX_ORDER1()
          CALL S_EDGE_EDGE_ORDER1_ANALYTIC(Se, V, ell)
          CALL T_EDGE_EDGE_ORDER1_ANALYTIC(Te, grad_lambda, M, ell, phi)
          ! Compute additional (hierarchal) H_1(curl) (LT/QN) terms if required.
          IF(element_order.EQ.2) THEN 
             N = INTEGRATION_MATRICES_ORDER2_N() 
             P = INTEGRATION_MATRICES_ORDER2_P()     
             CALL S_EDGE_EDGE_ORDER2_ANALYTIC(Se, V, ell)
             CALL S_EDGE_FACE_ORDER2_ANALYTIC(Se, V, ell)
             CALL S_FACE_FACE_ORDER2_ANALYTIC(Se, V, M)
             CALL T_EDGE_EDGE_ORDER2_ANALYTIC(Te, grad_lambda, M, ell, phi)
             CALL T_EDGE_FACE_ORDER2_ANALYTIC(Te, grad_lambda, N, ell, phi)
             CALL T_FACE_FACE_ORDER2_ANALYTIC(Te, grad_lambda, P, phi)
          ELSE IF(element_order.GE.3) THEN 
             STOP 'IE: Analytical S&T evaluation only available up to 2nd order elements'
          END IF
       ELSE 
          STOP 'IE: Analytical S&T evaluation only available for mixed order elements'
       ENDIF
    ENDIF
    ! Geometrical and material scaling (if requested).
    ! DBD changes 25 March 2003
    IF(.NOT.PRESENT(no_material_scaling)) THEN ! set to false as default if not provided. 
       no_material_scaling_flag = .false.
    ELSE
       no_material_scaling_flag = no_material_scaling ! Otherwise just copy it.
    END IF

    IF (no_material_scaling_flag) THEN
       IF(curvilinear_flag) THEN
          ! Curvilinear elements are defined on a parent tet with unit sides, volume 1/6. 
          ! Actual volume is taken into account in Jacobian. 
          Se = Se /6.0_SP
          Te = Te /6.0_SP
       ELSE 
          Se = Se * VOLUME(element_num) 
          Te = Te * VOLUME(element_num) 
       END IF
    ELSE ! Default
       IF(curvilinear_flag) THEN
          ! As above.
          !         Se = Se * VOLUME(element_num) / mu_r(elements(element_num)%material) ! While only T is done .....
          Se = Se  / (6.0_SP *mu_r(elements(element_num)%material) )
          Te = Te * eps_r(elements(element_num)%material) / 6.0_SP
       ELSE
          Se = Se * VOLUME(element_num) / mu_r(elements(element_num)%material)
          Te = Te * VOLUME(element_num) * eps_r(elements(element_num)%material)
       END IF
    ENDIF
    ! End DBD changes 25 March 2003


    IF(DEBUG_SYSTEM_ELEMENTS) THEN
       CALL PRINT_ELEMENT_S_AND_T
    END IF

    !if (element_num.eq.1) THEN
    !print *,Te
    !print *,Te
    !end if
    !stop

  CONTAINS
    SUBROUTINE PRINT_ELEMENT_S_AND_T
      USE nrtype
      USE unit_numbers
      IMPLICIT NONE
      !*******************************************************************************
      ! This utility routine prints out the  S and T entries for each individual 
      ! element. For testing and debugging only. 
      ! It is modified version of the routine originally used when eigenvalue
      ! analysis was being debugged, extended to support complete elements to LT/QN order.
      ! DBD 18 July 2005
      !*******************************************************************************
      INTEGER(I4B) jj

      WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER


      WRITE(FILEOUT,'(//A,A)') '****************** S matrix elements ',&
           '******************'                                                     

      WRITE(FILEOUT,'(A,I8)') 'Element number: ',element_num
      SELECT CASE (elements(element_num)%order) 
      CASE(1)
         IF(element_mixed_order) THEN
            DO jj = 1,6
               WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
               WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                    Se(jj,1:jj)
            END DO
         ELSE
            DO jj = 1,12
               WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
               WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                    Se(jj,1:jj)
            END DO
         END IF
      CASE(2)          
         DO jj = 1,20
            WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
            WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                 REAL(Se(jj,1:jj)) ! HACK
         END DO
      CASE DEFAULT
         STOP 'IE: Error in routine PRINT_ELEMENT_S_AND_T - unimplemented hierarchal order'
      END SELECT

      WRITE(FILEOUT,'(//A,A)') '****************** T matrix elements ',&
           '******************'

      WRITE(FILEOUT,'(A,I8)') 'Element number: ',element_num
      SELECT CASE (elements(element_num)%order) 
      CASE(1)                    
         IF(element_mixed_order) THEN
            DO jj = 1,6
               WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
               WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                    Te(jj,1:jj)
            END DO
         ELSE      
            DO jj = 1,12
               WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
               WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                    Te(jj,1:jj)
            END DO
         END IF
      CASE(2)          
         DO jj = 1,20
            WRITE(FILEOUT,'(A,I8)') 'Row = ',jj
            WRITE(FILEOUT,'(3(G10.4,1X,G10.4,3X),/,3X,3(G10.4,1X,G10.4,3X))' ) & 
                 REAL(Te(jj,1:jj)) ! HACK
         END DO
      CASE DEFAULT
         STOP 'IE: Error in routine PRINT_ELEMENT_S_AND_T - unimplemented hierarchal order'
      END SELECT

    END SUBROUTINE PRINT_ELEMENT_S_AND_T

  END SUBROUTINE S_AND_T_MAKE_HIERARCHAL

  SUBROUTINE S_and_T_make_V_phi(V, phi, element_num)
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    REAL(SP), DIMENSION(4,4,3), INTENT(out) ::  V  ! See [eqn(7),S&P]    
    REAL(SP), DIMENSION(4,4), INTENT(out) ::   phi  ! See [eqn(8),S&P] 
    INTEGER(i4b), INTENT(in) :: element_num
!!!
!!! Local Variables
!!!
    INTEGER(I4B) :: inode,jnode        ! Loop counters
    REAL(SP), DIMENSION(4) :: a,b,c,d  ! Gradient of basis function terms: 
    ! See [eqn(6),S&P] 

    CALL SIMPLEX_COEFFICIENTS(element_num,a,b,c,d)    

    ! Compute dot product terms:  See [eqns(7)&(8),S&P] 
    ! Needed for both closed form and cubature implementations. 
    DO inode = 1,4 ! Node
       DO jnode = 1,4 ! Node
          phi(inode,jnode) = b(inode)*b(jnode) + c(inode)*c(jnode) + & 
               d(inode)*d(jnode)
          V(inode,jnode,1) = c(inode)*d(jnode) - c(jnode)*d(inode) ! V_x 
          V(inode,jnode,2) = b(jnode)*d(inode) - b(inode)*d(jnode) ! V_y
          V(inode,jnode,3) = b(inode)*c(jnode) - b(jnode)*c(inode) ! V_z
       END DO
    END DO
  END SUBROUTINE S_and_T_make_V_phi

  SUBROUTINE FILL_LOWER_ELEMENTAL_MATRICES(Se, matrix_size)
    !*******************************************************************************
    ! Fills in (uncomputed) lower half of S and T matrices. Note that some elements are
    ! overwritten - this is for coding simplicity. 
    !*******************************************************************************
    IMPLICIT NONE
    COMPLEX(SPC), INTENT(OUT), DIMENSION(:,:) :: Se
    INTEGER(i4b), INTENT(in) :: matrix_size
    INTEGER(i4b) ii,jj
    DO ii=2,matrix_size
       DO jj = 1,ii-1
          Se(ii,jj) = Se(jj,ii)
       END DO
    END DO
  END SUBROUTINE FILL_LOWER_ELEMENTAL_MATRICES

  SUBROUTINE S_EDGE_EDGE(Se, elinfo_in, V_in, curvilinear_flag, sigma_tensor_in)

    !*******************************************************************************
    ! Computes edge-edge S submatrices using quadrature 
    ! and stores in the elemental S matrix.
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    ! Changed 11 June 2003 DBD: tensor material property included to CT/LN order 
    ! Extended 20 July 2005 DBD to support 2nd order curvilinear elements.
    !*******************************************************************************
    USE integrate
    USE S_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,4,3), INTENT(in) ::   V_in  ! See [eqn(7),S&P]
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Se1e1
    COMPLEX(SPC), DIMENSION(6,6) :: Se1e2, Se2e2
    COMPLEX(SPC), DIMENSION(6,6) :: Se1e3, Se2e3, Se3e3
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) :: intg_order  ! Desired integration order
    LOGICAL(lgt)  :: tensor_material_flag  
    ! Jacobian (or inverse) for curvlinear elements; temporary matrix using these.

    tensor_material_flag = PRESENT(sigma_tensor_in)
!!!
!!! Assign input parameters to module S_intg_funs
!!!
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    V = V_in
    elinfo = elinfo_in
!!!
    IF (SCALE_BY_EDGE_LENGTH) CALL ERROR_FEMFEKO(1,4106)

    IF(.NOT.curvilinear_flag) THEN
       intg_order = 2              ! 2nd order rule [Jinyun]. Function    
       ! is at most quadratic - for e2e2 elements.   
    ELSE 
       intg_order = 4              
    END IF

    Se1e1(1:6,1:6) = 0.0_SP

    Se1e2(1:6,1:6) = 0.0_SP 
    Se2e2(1:6,1:6) = 0.0_SP 

    Se1e3(1:6,1:6) = 0.0_SP
    Se2e3(1:6,1:6) = 0.0_SP
    Se3e3(1:6,1:6) = 0.0_SP

    EDGE_LOOP1: DO iedge = 1,6
       EDGE_LOOP2: DO jedge = 1,6
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)
          IF(.NOT.tensor_material_flag) THEN
             IF(.NOT.curvilinear_flag) THEN
                Se1e1(iedge,jedge) = integrate_tet(curl_e1e1, intg_order)
             ELSE 
                Se1e1(iedge,jedge) = integrate_tet(curl_e1e1_curvi, intg_order)
             END IF
          ELSE 
             Se1e1(iedge,jedge) = integrate_tet(curl_e1e1_tensor, intg_order)
          END IF
          IF (elinfo%order.GE.2.AND.elinfo%mixed_order.OR.&
               elinfo%order.GE.1.AND..NOT.elinfo%mixed_order) THEN
             IF(.NOT.curvilinear_flag) THEN
                Se1e2(iedge,jedge) = integrate_tet(curl_e1e2, intg_order)
                Se2e2(iedge,jedge) = integrate_tet(curl_e2e2, intg_order)
             ELSE 
                IF (element_type.EQ.3) THEN
                   Se1e2 = 0.0 ! Known to be zero.
                   Se2e2 = 0.0 ! Known to be zero.
                ELSE
                   STOP 'IE: Only Webb LT/LN curvilinear elements implemented in S_EDGE_EDGE.'
                END IF
             END IF
          END IF
       END DO EDGE_LOOP2
    END DO EDGE_LOOP1
    IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
         elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
       IF (element_type.EQ.3) THEN
          Se1e3 = 0.0 ! Known to be zero.
          Se2e3 = 0.0 ! Known to be zero.
          Se3e3 = 0.0 ! Known to be zero.
       ELSE
          STOP 'IE: Only Webb QT/QN elements implemented in S_EDGE_EDGE'
       END IF
    END IF
    Se(1:6,1:6) = Se1e1(1:6,1:6)
    ! Fill elemental matrix. Note that for the lowest order elements,
    ! the full matrix is filled, but only the first 6x6 entries are non-zero.
    Se(1:6,7:12)    = Se1e2(1:6,1:6)
    Se(7:12,7:12)   = Se2e2(1:6,1:6)

    Se(1:6,21:26)   = Se1e3(1:6,1:6)
    Se(7:12,21:26)  = Se2e3(1:6,1:6)
    Se(21:26,21:26) = Se3e3(1:6,1:6)

  END SUBROUTINE S_EDGE_EDGE


  SUBROUTINE S_EDGE_FACE(Se, elinfo_in, V_in, curvilinear_flag, sigma_tensor_in)
    !*******************************************************************************
    ! Computes edge-face S submatrices and stores in elemental S matrix. 
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    ! Extended 20 July 2005 DBD to support 2nd order curvilinear elements.
    !*******************************************************************************
    USE integrate
    USE S_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,4,3) ::   V_in  ! See [eqn(7),S&P]
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,4) :: Se1f1, Se1f2, Se2f1, Se2f2
    COMPLEX(SPC), DIMENSION(6,4) :: Se1f3, Se2f3
    COMPLEX(SPC), DIMENSION(6,4) :: Se3f1, Se3f2, Se3f3
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes 
    INTEGER(I4B), DIMENSION(3) :: j_face_nodes
    REAL(SP), DIMENSION(3) :: temp2    ! temp. storage                            
    INTEGER(I4B) :: intg_order  ! Desired integration order
    LOGICAL(lgt)  :: tensor_material_flag  

!!!
!!! Assign input parameters to module S_intg_funs
!!!
    tensor_material_flag = PRESENT(sigma_tensor_in)
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    V = V_in
    elinfo = elinfo_in
!!!
    IF (SCALE_BY_EDGE_LENGTH) CALL ERROR_FEMFEKO(1,4106)

    IF(.NOT.curvilinear_flag) THEN
       intg_order = 4           ! 4-th order rule.
    ELSE 
       intg_order = 6           ! 6-th order rule
    END IF

    ! Initialize (needed for cubature)
    Se1f1(1:6,1:4) = 0.0_SP
    Se1f2(1:6,1:4) = 0.0_SP
    Se2f1(1:6,1:4) = 0.0_SP
    Se2f2(1:6,1:4) = 0.0_SP

    Se1f3(1:6,1:4) = 0.0_SP
    Se2f3(1:6,1:4) = 0.0_SP

    Se3f1(1:6,1:4) = 0.0_SP
    Se3f2(1:6,1:4) = 0.0_SP
    Se3f3(1:6,1:4) = 0.0_SP

    ! Fill edge-face submatrices
    EDGELOOP: DO iedge = 1,6
       FACELOOP: DO jface = 1,4 
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)

          IF(.NOT.curvilinear_flag) THEN
             Se1f1(iedge,jface) = integrate_tet(curl_e1f1, intg_order)
             Se1f2(iedge,jface) = integrate_tet(curl_e1f2, intg_order)
          ELSE 
             Se1f1(iedge,jface) = integrate_tet(curl_e1f1_curvi, intg_order)
             Se1f2(iedge,jface) = integrate_tet(curl_e1f2_curvi, intg_order)
          END IF

          IF(.NOT.curvilinear_flag) THEN
             Se2f1(iedge,jface) = integrate_tet(curl_e2f1, intg_order)
             Se2f2(iedge,jface) = integrate_tet(curl_e2f2, intg_order)
          ELSE 
             Se2f1(iedge,jface) = integrate_tet(curl_e2f1_curvi, intg_order)
             Se2f2(iedge,jface) = integrate_tet(curl_e2f2_curvi, intg_order)
          END IF
       END DO FACELOOP
    END DO EDGELOOP
    IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
         elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
       IF (element_type.EQ.3) THEN
          Se1f3 = 0.0 ! Known to be zero.
          Se2f3 = 0.0 !   Ditto
          Se3f1 = 0.0 ! 
          Se3f2 = 0.0 ! 
          Se3f3 = 0.0 ! 
       ELSE
          STOP 'IE: Only Webb QT/QN elements implemented in S_EDGE_EDGE'
       END IF
    END IF
    Se(1:6,13:16)   = Se1f1(1:6,1:4)
    Se(1:6,17:20)   = Se1f2(1:6,1:4)     
    Se(7:12,13:16)  = Se2f1(1:6,1:4)
    Se(7:12,17:20)  = Se2f2(1:6,1:4)

    Se(1:6,27:30)   = Se1f3(1:6,1:4)
    Se(7:12,27:30)  = Se2f3(1:6,1:4)
    Se(21:26,13:16) = Se3f1(1:6,1:4)
    Se(21:26,17:20) = Se3f2(1:6,1:4)
    Se(21:26,27:30) = Se3f3(1:6,1:4)

  END SUBROUTINE S_EDGE_FACE

  SUBROUTINE S_FACE_FACE(Se, elinfo_in, V_in, curvilinear_flag, sigma_tensor_in)
    !*******************************************************************************
    ! Computes face-face S submatrices and stores in elemental S matrix. 
    ! others are filled using symmetry.
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    ! Extended 20 July 2005 DBD to support 2nd order curvilinear elements.
    !*******************************************************************************
    USE integrate
    USE S_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,4,3) ::   V_in  ! See [eqn(7),S&P]
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(4,4) :: Sf1f1, Sf1f2, Sf2f2
    COMPLEX(SPC), DIMENSION(4,4) :: Sf1f3, Sf2f3, Sf3f3
    INTEGER(I4B), DIMENSION(3) :: i_face_nodes, j_face_nodes
    INTEGER(I4B) intg_order            ! Integration order to use
    LOGICAL(lgt)  :: tensor_material_flag  

    tensor_material_flag = PRESENT(sigma_tensor_in)
!!!
!!! Assign input parameters to module S_intg_funs
!!!
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    V = V_in
    elinfo = elinfo_in
!!!
    IF(.NOT.curvilinear_flag) THEN
       intg_order = 4 ! 4-th order rule.
    ELSE 
       intg_order = 6 ! 6-th order rule
    END IF

    ! Initialize (needed for cubature)
    Sf1f1(1:4,1:4) = 0.0_SP
    Sf1f2(1:4,1:4) = 0.0_SP
    Sf2f2(1:4,1:4) = 0.0_SP

    Sf1f3(1:4,1:4) = 0.0_SP
    Sf2f3(1:4,1:4) = 0.0_SP
    Sf3f3(1:4,1:4) = 0.0_SP

    FACELOOP1: DO iface = 1,4
       FACELOOP2: DO jface = 1,4 
          i_face_nodes = LOCAL_FACENODES(iface)
          i1 = i_face_nodes(1)
          i2 = i_face_nodes(2)
          i3 = i_face_nodes(3)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)
          IF(.NOT.curvilinear_flag) THEN
             Sf1f1(iface,jface) = integrate_tet(curl_f1f1, intg_order)
             Sf1f2(iface,jface) = integrate_tet(curl_f1f2, intg_order)
             Sf2f2(iface,jface) = integrate_tet(curl_f2f2, intg_order)
          ELSE 
             Sf1f1(iface,jface) = integrate_tet(curl_f1f1_curvi, intg_order)
             Sf1f2(iface,jface) = integrate_tet(curl_f1f2_curvi, intg_order)
             Sf2f2(iface,jface) = integrate_tet(curl_f2f2_curvi, intg_order)
          END IF
       END DO FACELOOP2
    END DO FACELOOP1
    IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
         elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
       IF (element_type.EQ.3) THEN
          Sf1f3 = 0.0 ! Known to be zero.
          Sf2f3 = 0.0 !   Ditto
          Sf3f3 = 0.0 ! 
       ELSE
          STOP 'IE: Only Webb QT/QN elements implemented in S_EDGE_EDGE'
       END IF
    END IF
    Se(13:16,13:16)  = Sf1f1(1:4,1:4)    
    Se(13:16,17:20)  = Sf1f2(1:4,1:4)     
    Se(17:20,17:20)  = Sf2f2(1:4,1:4) 

    Se(13:16,27:30)  = Sf1f3(1:4,1:4)    
    Se(17:20,27:30)  = Sf2f3(1:4,1:4)    
    Se(27:30,27:30)  = Sf3f3(1:4,1:4)    

  END SUBROUTINE S_FACE_FACE

  SUBROUTINE T_EDGE_EDGE(Te, elinfo_in, grad_lambda_in, curvilinear_flag, sigma_tensor_in)
    !*******************************************************************************
    ! Computes edge-edge T submatrices using quadrature 
    ! and stores in the elemental T matrix.
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    ! Changed 11 June 2003 DBD: tensor material property included to CT/LN order 
    !*******************************************************************************
    USE integrate
    USE T_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,3), INTENT(in) :: grad_lambda_in ! Gradient of simplex coordinates
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Te1e1
    COMPLEX(SPC), DIMENSION(6,6) :: Te1e2, Te2e2
    COMPLEX(SPC), DIMENSION(6,6) :: Te1e3, Te2e3, Te3e3
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) intg_order            ! Requested integration order
    LOGICAL(lgt)  :: tensor_material_flag  
!!!
!!! Assign input parameters to module T_intg_funs
!!!
    tensor_material_flag = PRESENT(sigma_tensor_in)
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    grad_lambda = grad_lambda_in
    elinfo = elinfo_in
!!!    
    IF (SCALE_BY_EDGE_LENGTH) CALL ERROR_FEMFEKO(1,4106)

    IF(.NOT.curvilinear_flag) THEN
       intg_order = 4 ! 4-th order rule.
    ELSE 
       intg_order = 6 ! 6-th order rule
    END IF

    ! Initialize (needed for cubature)
    Te1e1(1:6,1:6) = 0.0_SP

    Te1e2(1:6,1:6) = 0.0_SP 
    Te2e2(1:6,1:6) = 0.0_SP 

    Te1e3(1:6,1:6) = 0.0_SP
    Te2e3(1:6,1:6) = 0.0_SP
    Te3e3(1:6,1:6) = 0.0_SP


    EDGE_LOOP1: DO iedge = 1,6
       EDGE_LOOP2: DO jedge = 1,6
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)
          IF(.NOT.tensor_material_flag) THEN
             IF(.NOT.curvilinear_flag) THEN
                ! Standard rectilinear elements
                Te1e1(iedge,jedge) = integrate_tet(e1e1, intg_order)
             ELSE 
                ! Curvilinear elements
                Te1e1(iedge,jedge) = integrate_tet(e1e1_curvi, intg_order)
             END IF
          ELSE 
             ! Rectilinear elements with material tensor properties
             Te1e1(iedge,jedge) = integrate_tet(e1e1_tensor, intg_order)
          END IF

          IF (elinfo%order.GE.2.AND.elinfo%mixed_order.OR.&
               elinfo%order.GE.1.AND..NOT.elinfo%mixed_order) THEN
             IF(.NOT.curvilinear_flag) THEN
                Te1e2(iedge,jedge)  = integrate_tet(e1e2, intg_order)
                Te2e2(iedge,jedge)  = integrate_tet(e2e2, intg_order)
             ELSE 
                Te1e2(iedge,jedge) = integrate_tet(e1e2_curvi, intg_order)
                Te2e2(iedge,jedge) = integrate_tet(e2e2_curvi, intg_order)
             END IF

          END IF
          IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
               elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
             IF (element_type.EQ.3) THEN
                IF(.NOT.curvilinear_flag) THEN
                   Te1e3(iedge,jedge) = integrate_tet(e1e3, intg_order)
                   Te2e3(iedge,jedge) = integrate_tet(e2e3, intg_order)
                   Te3e3(iedge,jedge) = integrate_tet(e3e3, intg_order)
                ELSE 
                   Te1e3(iedge,jedge) = integrate_tet(e1e3_curvi, intg_order)
                   Te2e3(iedge,jedge) = integrate_tet(e2e3_curvi, intg_order)
                   Te3e3(iedge,jedge) = integrate_tet(e3e3_curvi, intg_order)
                END IF

             ELSE
                STOP 'IE: Only Webb QT/QN elements implemented in T_EDGE_EDGE'
             END IF
          END IF
       END DO EDGE_LOOP2
    END DO EDGE_LOOP1
    Te(1:6,1:6) = Te1e1(1:6,1:6)
    ! Fill elemental matrix. Note that for the lowest order elements,
    ! the full matrix is filled, but only the first 6x6 entries are non-zero.
    Te(1:6,7:12)    = Te1e2(1:6,1:6)
    Te(7:12,7:12)   = Te2e2(1:6,1:6)

    Te(1:6,21:26)   = Te1e3(1:6,1:6)
    Te(7:12,21:26)  = Te2e3(1:6,1:6)
    Te(21:26,21:26) = Te3e3(1:6,1:6)

  END SUBROUTINE T_EDGE_EDGE


  SUBROUTINE T_EDGE_FACE(Te, elinfo_in, grad_lambda_in, curvilinear_flag, sigma_tensor_in)
    !*******************************************************************************
    ! Fills edge-face T submatrices.
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    ! Extended 16 July 2005 DBD to support 2nd order curvilinear elements.
    !*******************************************************************************
    USE integrate
    USE T_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,3), INTENT(in)  :: grad_lambda_in ! Gradient of simplex coordinates
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,4) :: Te1f1, Te1f2, Te2f1, Te2f2
    COMPLEX(SPC), DIMENSION(6,4) :: Te1f3, Te2f3
    COMPLEX(SPC), DIMENSION(6,4) :: Te3f1, Te3f2, Te3f3

    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes
    INTEGER(I4B), DIMENSION(3) :: j_face_nodes
    INTEGER(I4B) intg_order            ! Requested integration order
    LOGICAL(lgt)  :: tensor_material_flag  
!!!
!!! Assign input parameters to module T_intg_funs
!!!
    tensor_material_flag = PRESENT(sigma_tensor_in)
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    grad_lambda = grad_lambda_in
    elinfo = elinfo_in
!!!
    IF (SCALE_BY_EDGE_LENGTH) CALL ERROR_FEMFEKO(1,4106)

    IF(.NOT.curvilinear_flag) THEN
       intg_order = 4 ! 4-th order rule.
    ELSE 
       intg_order = 6 ! 6-th order rule
    END IF

    ! Initialize (needed for cubature)
    Te1f1(1:6,1:4) = 0.0_SP
    Te1f2(1:6,1:4) = 0.0_SP
    Te2f1(1:6,1:4) = 0.0_SP
    Te2f2(1:6,1:4) = 0.0_SP

    Te1f3(1:6,1:4) = 0.0_SP
    Te2f3(1:6,1:4) = 0.0_SP

    Te3f1(1:6,1:4) = 0.0_SP
    Te3f2(1:6,1:4) = 0.0_SP
    Te3f3(1:6,1:4) = 0.0_SP

    EDGE_LOOP: DO iedge = 1,6
       FACE_LOOP: DO jface = 1,4 
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)               
          IF(.NOT.curvilinear_flag) THEN
             Te1f1(iedge,jface) = integrate_tet(e1f1, intg_order)
             Te2f1(iedge,jface) = integrate_tet(e2f1, intg_order)
             Te1f2(iedge,jface) = integrate_tet(e1f2, intg_order)
             Te2f2(iedge,jface) = integrate_tet(e2f2, intg_order)
          ELSE 
             Te1f1(iedge,jface) = integrate_tet(e1f1_curvi, intg_order)
             Te2f1(iedge,jface) = integrate_tet(e2f1_curvi, intg_order)
             Te1f2(iedge,jface) = integrate_tet(e1f2_curvi, intg_order)
             Te2f2(iedge,jface) = integrate_tet(e2f2_curvi, intg_order)
          END IF
          IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
               elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
             IF (element_type.EQ.3) THEN
                IF(.NOT.curvilinear_flag) THEN
                   Te1f3(iedge,jface) =  integrate_tet(e1f3, intg_order)
                   Te2f3(iedge,jface) =  integrate_tet(e2f3, intg_order)
                   Te3f1(iedge,jface) =  integrate_tet(e3f1, intg_order)
                   Te3f2(iedge,jface) =  integrate_tet(e3f2, intg_order)
                   Te3f3(iedge,jface) =  integrate_tet(e3f3, intg_order)
                ELSE 
                   Te1f3(iedge,jface) = integrate_tet(e1f3_curvi, intg_order)
                   Te2f3(iedge,jface) = integrate_tet(e2f3_curvi, intg_order)
                   Te3f1(iedge,jface) = integrate_tet(e3f1_curvi, intg_order)
                   Te3f2(iedge,jface) = integrate_tet(e3f2_curvi, intg_order)
                   Te3f3(iedge,jface) = integrate_tet(e3f3_curvi, intg_order)
                END IF
             ELSE
                STOP 'IE: Only Webb QT/QN elements implemented in T_EDGE_FACE'
             END IF
          END IF
       END DO FACE_LOOP
    END DO EDGE_LOOP

    Te(1:6,13:16)   = Te1f1(1:6,1:4)
    Te(1:6,17:20)   = Te1f2(1:6,1:4)     
    Te(7:12,13:16)  = Te2f1(1:6,1:4)
    Te(7:12,17:20)  = Te2f2(1:6,1:4)

    Te(1:6,27:30)   = Te1f3(1:6,1:4)
    Te(7:12,27:30)  = Te2f3(1:6,1:4)
    Te(21:26,13:16) = Te3f1(1:6,1:4)
    Te(21:26,17:20) = Te3f2(1:6,1:4)
    Te(21:26,27:30) = Te3f3(1:6,1:4)

  END SUBROUTINE T_EDGE_FACE

  SUBROUTINE T_FACE_FACE(Te, elinfo_in, grad_lambda_in, curvilinear_flag, sigma_tensor_in)
    !*******************************************************************************
    ! Fills face-face T submatrices.
    ! Presently implemented for up 2nd order QT/QN elements. 
    ! Returns upper half of matrix. 
    !*******************************************************************************
    USE integrate
    USE T_intg_funs
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    TYPE(element_info_type), INTENT(in) :: elinfo_in
    REAL(SP), DIMENSION(4,3) :: grad_lambda_in ! Gradient of simplex coordinates
    LOGICAL(lgt), INTENT(in) :: curvilinear_flag
    REAL(SP), DIMENSION(3,3), INTENT(IN), OPTIONAL  :: sigma_tensor_in
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(4,4) :: Tf1f1, Tf1f2, Tf2f2
    COMPLEX(SPC), DIMENSION(4,4) :: Tf1f3, Tf2f3, Tf3f3
    INTEGER(I4B), DIMENSION(3) :: i_face_nodes, j_face_nodes
    INTEGER(I4B) intg_order            ! Requested integration order
    LOGICAL(lgt)  :: tensor_material_flag  
!!!
!!! Assign input parameters to module T_intg_funs
!!!
    tensor_material_flag = PRESENT(sigma_tensor_in)
    IF (tensor_material_flag) sigma_tensor = sigma_tensor_in
    grad_lambda = grad_lambda_in
    elinfo = elinfo_in
!!!
    IF(.NOT.curvilinear_flag) THEN
       intg_order = 4 ! 4-th order rule.
    ELSE 
       intg_order = 6 ! 6-th order rule
    END IF

    ! Initialize (needed for cubature)
    Tf1f1(1:4,1:4) = 0.0_SP
    Tf1f2(1:4,1:4) = 0.0_SP
    Tf2f2(1:4,1:4) = 0.0_SP

    Tf1f3(1:4,1:4) = 0.0_SP
    Tf2f3(1:4,1:4) = 0.0_SP
    Tf3f3(1:4,1:4) = 0.0_SP

    FACE_LOOP1: DO iface = 1,4
       FACE_LOOP2: DO jface = 1,4 
          i_face_nodes = LOCAL_FACENODES(iface)
          i1 = i_face_nodes(1)
          i2 = i_face_nodes(2)
          i3 = i_face_nodes(3)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)
          IF(.NOT.curvilinear_flag) THEN
             Tf1f1(iface,jface) = integrate_tet(f1f1, intg_order)
             Tf1f2(iface,jface) = integrate_tet(f1f2, intg_order) 
             Tf2f2(iface,jface) = integrate_tet(f2f2, intg_order) 
          ELSE
             Tf1f1(iface,jface) = integrate_tet(f1f1_curvi, intg_order) 
             Tf1f2(iface,jface) = integrate_tet(f1f2_curvi, intg_order) 
             Tf2f2(iface,jface) = integrate_tet(f2f2_curvi, intg_order) 
          END IF
          IF (elinfo%order.GE.3.AND.elinfo%mixed_order.OR.&
               elinfo%order.GE.2.AND..NOT.elinfo%mixed_order) THEN 
             IF (element_type.EQ.3) THEN
                IF(.NOT.curvilinear_flag) THEN
                   Tf1f3(iface,jface) = integrate_tet(f1f3, intg_order)
                   Tf2f3(iface,jface) = integrate_tet(f2f3, intg_order)
                   Tf3f3(iface,jface) = integrate_tet(f3f3, intg_order)
                ELSE 
                   Tf1f3(iface,jface) = integrate_tet(f1f3_curvi, intg_order)
                   Tf2f3(iface,jface) = integrate_tet(f2f3_curvi, intg_order)
                   Tf3f3(iface,jface) = integrate_tet(f3f3_curvi, intg_order)
                END IF
             ELSE
                STOP 'IE: Only Webb QT/QN elements implemented in T_FACE_FACE'
             END IF
          END IF
       END DO FACE_LOOP2
    END DO FACE_LOOP1
    Te(13:16,13:16)  = Tf1f1(1:4,1:4)    
    Te(13:16,17:20)  = Tf1f2(1:4,1:4)     
    Te(17:20,17:20)  = Tf2f2(1:4,1:4) 

    Te(13:16,27:30)  = Tf1f3(1:4,1:4)    
    Te(17:20,27:30)  = Tf2f3(1:4,1:4)    
    Te(27:30,27:30)  = Tf3f3(1:4,1:4)    

  END SUBROUTINE T_FACE_FACE

  ! Following group of routines will probably eventually be removed. 
  ! Retained for speed tests to compare quadrature to analytical 
  ! evaluation. 

  FUNCTION INTEGRATION_MATRIX_ORDER1()
    !*******************************************************************************
    ! Sets up integration matrices M  ([eqn. (15) S&P]).  
    !*******************************************************************************
    REAL(SP), DIMENSION(4,4) :: INTEGRATION_MATRIX_ORDER1    ! integration matrix M [eqn(15),S&P]
    INTEGRATION_MATRIX_ORDER1 = SPREAD(SPREAD(1.0_SP,1,4),2,4) ! Create a 4x4 matrix of ones
    INTEGRATION_MATRIX_ORDER1(1,1) = 2.0_SP ! Change diagonals
    INTEGRATION_MATRIX_ORDER1(2,2) = 2.0_SP
    INTEGRATION_MATRIX_ORDER1(3,3) = 2.0_SP
    INTEGRATION_MATRIX_ORDER1(4,4) = 2.0_SP
    INTEGRATION_MATRIX_ORDER1= INTEGRATION_MATRIX_ORDER1 / 20.0_SP ! Normalize
  END FUNCTION INTEGRATION_MATRIX_ORDER1


  FUNCTION INTEGRATION_MATRICES_ORDER2_N()
    !*******************************************************************************
    ! Sets up integration matrix N. (See S&P).  
    !*******************************************************************************
    ! Set up N matrix. Note that N is such that the 2nd and 3rd nodes
    ! are different by definition.
    IMPLICIT NONE
    REAL(SP), DIMENSION(4,4,4) :: INTEGRATION_MATRICES_ORDER2_N ! integration matrix N [eqn(46),S&P]
    INTEGER(i4b) :: inode, jnode, knode
    DO inode  = 1,4
       DO jnode = 1,4
          DO knode = 1,4
             IF (inode.EQ.jnode.OR.inode.EQ.knode) THEN ! One coincident node
                INTEGRATION_MATRICES_ORDER2_N(inode,jnode,knode) = 2.0_SP
             ELSE ! Disjoint nodes.
                INTEGRATION_MATRICES_ORDER2_N(inode,jnode,knode) = 1.0_SP
             END IF
          END DO
       END DO
    END DO
    INTEGRATION_MATRICES_ORDER2_N = INTEGRATION_MATRICES_ORDER2_N/120.0_SP    ! Scaling factor
  END FUNCTION INTEGRATION_MATRICES_ORDER2_N

  FUNCTION INTEGRATION_MATRICES_ORDER2_P()
    !*******************************************************************************
    ! Sets up integration matrix P. (See S&P).  
    !*******************************************************************************
    IMPLICIT NONE
    REAL(SP), DIMENSION(4,4,4,4) :: INTEGRATION_MATRICES_ORDER2_P !integration matrix P [eqn(47),S&P]
    INTEGER(i4b) :: inode, jnode, knode, lnode
    ! Set up P matrix
    DO inode  = 1,4
       DO jnode = 1,4
          DO knode = 1,4
             DO lnode = 1,4
                IF (inode.EQ.knode.AND.jnode.EQ.lnode) THEN ! Coincident edge     
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 4.0_SP
                ELSE IF (inode.EQ.knode) THEN ! Coincident start nodes
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 2.0_SP
                ELSE IF (jnode.EQ.lnode) THEN ! Coincident end nodes
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 2.0_SP
                ELSE IF (jnode.EQ.knode) THEN ! Coincident end and start nodes
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 2.0_SP
                ELSE IF (inode.EQ.lnode) THEN ! Coincident start and end nodes
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 2.0_SP
                ELSE ! Completely disjoint nodes.
                   INTEGRATION_MATRICES_ORDER2_P(inode,jnode,knode,lnode) = 1.0_SP
                END IF
             END DO
          END DO
       END DO
    END DO
    INTEGRATION_MATRICES_ORDER2_P = INTEGRATION_MATRICES_ORDER2_P/840.0_SP    ! Scaling factor
  END FUNCTION INTEGRATION_MATRICES_ORDER2_P

  SUBROUTINE S_EDGE_EDGE_ORDER1_ANALYTIC(Se, V, ell)
    !*******************************************************************************
    ! Fills edge-edge S submatrices for the lowest order elements and 
    ! stores in elemental S matrix (E in [eqn.(10), S&P]).
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
    REAL(SP), DIMENSION(4,4,3), INTENT(in) ::   V  ! See [eqn(7),S&P]
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Se1e1
    INTEGER(I4B) :: iedge,jedge        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 4 ! 2nd order rule [Jinyun]. Note that this is over-kill for 
    ! most lowest-order vbf's; the gradient of these is usually a constant.
    ! Initialize (needed for cubature)
    Se1e1 = 0.0_SP
    EDGE_LOOP1: DO iedge = 1,6
       EDGE_LOOP2: DO jedge = 1,iedge ! Fill lower half
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)

          IF(.NOT.CUBATURE) THEN
             Se1e1(iedge,jedge) = 4* DOT_PRODUCT(V(i1,i2,1:3),V(j1,j2,1:3))
          ELSE 
             ! Integrate numerically, 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Se1e1(iedge,jedge) = Se1e1(iedge,jedge) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * & 
                     DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(1,lambda,V,j1,j2))
             END DO
          END IF
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             Se1e1(iedge,jedge) = ell(iedge)*ell(jedge)*Se1e1(iedge,jedge)
          ENDIF
          Se1e1(jedge,iedge) = Se1e1(iedge,jedge) ! Symmetrize. Diagonal is overwritten.
       END DO EDGE_LOOP2
    END DO EDGE_LOOP1
    Se(1:6,1:6) = Se1e1(1:6,1:6)     
  END SUBROUTINE S_EDGE_EDGE_ORDER1_ANALYTIC

  SUBROUTINE S_EDGE_EDGE_ORDER2_ANALYTIC(Se, V, ell)
    !*******************************************************************************
    ! Fills additional edge-edge S submatrices for the LT/QN elements and 
    ! stores in elemental S matrix. Those computed are Se1e2 and Se2e2; Se2e1 is 
    ! filled using symmetry.
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
    REAL(SP), DIMENSION(4,4,3), INTENT(in) ::   V  ! See [eqn(7),S&P]
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Se1e2, Se2e1, Se2e2
    INTEGER(I4B) :: iedge,jedge        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 4 ! 2nd order rule [Jinyun]. Function 
    ! is at most quadratic - for e2e2 elements. 
    ! Initialize (needed for cubature)
    Se1e2(1:6,1:6) = 0.0_SP 
    Se2e2(1:6,1:6) = 0.0_SP 

    EDGELOOP1: DO iedge = 1,6
       EDGELOOP2: DO jedge = 1,6
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)
          SELECT CASE(ELEMENT_TYPE) 
          CASE(1,3) ! S&P and Webb99
             Se1e2(1:6,1:6) = 0.0_SP ! Remains zero
             Se2e2(1:6,1:6) = 0.0_SP 
          CASE(2) ! A&V
             ! Analytical expression not available. 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Se1e2(iedge,jedge) = Se1e2(iedge,jedge) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(2,lambda,V,j1,j2))
                Se2e2(iedge,jedge) = Se2e2(iedge,jedge) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(2,lambda,V,j1,j2))
             END DO
          CASE DEFAULT
             STOP 'IE: Invalid element type in S_AND_T_MAKE_HIERARCHAL'
          END SELECT
          ! Note - the A&V elements are NOT scaled by edge-length.
       END DO EDGELOOP2
    END DO EDGELOOP1
    ! Fill other sub-matrices by symmetry.
    Se2e1 = TRANSPOSE(Se1e2)
    Se(1:6,7:12)  = Se1e2(1:6,1:6)
    Se(7:12,1:6)  = Se2e1(1:6,1:6)  
    Se(7:12,7:12) = Se2e2(1:6,1:6)   
  END SUBROUTINE S_EDGE_EDGE_ORDER2_ANALYTIC

  SUBROUTINE S_EDGE_FACE_ORDER2_ANALYTIC(Se, V, ell)
    !*******************************************************************************
    ! Fills edge-face S submatrices for the LT/QN elements and 
    ! stores in elemental S matrix. These computed are Se1f1,Se1f2,Se2f1,Se2f2; the 
    ! others are filled using symmetry. 
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
    REAL(SP), DIMENSION(4,4,3), INTENT(in) ::   V  ! See [eqn(7),S&P]
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,4) :: Se1f1, Se1f2, Se2f1, Se2f2
    COMPLEX(SPC), DIMENSION(4,6) :: Sf1e1, Sf1e2, Sf2e1, Sf2e2
    INTEGER(I4B) :: iedge,jface        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2,j3
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes 
    INTEGER(I4B), DIMENSION(3) :: j_face_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.
    REAL(SP), DIMENSION(3) :: temp2    ! temp. storage                            

    num_cpts = 11 ! Use 4-th order rule.

    ! Initialize (needed for cubature)
    Se1f1(1:6,1:4) = 0.0_SP
    Se1f2(1:6,1:4) = 0.0_SP
    Se2f1(1:6,1:4) = 0.0_SP
    Se2f2(1:6,1:4) = 0.0_SP

    ! Fill edge-face submatrices
    EDGELOOP: DO iedge = 1,6
       FACELOOP: DO jface = 1,4 
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)
          IF (.NOT.CUBATURE) THEN
             IF (ELEMENT_TYPE.NE.1) THEN
                STOP 'IE: Analytical integration only implemented for S&P LT/QN elements'
             END IF
             temp2(1:3) = 2.0_SP*V(j2,j3,1:3) + V(j1,j3,1:3) - V(j1,j2,1:3)               
             Se1f1(iedge,jface) = 0.5_SP*DOT_PRODUCT(V(i1,i2,1:3),temp2(1:3))  
             temp2(1:3) = V(j2,j3,1:3) + 2.0_SP*V(j1,j3,1:3) + V(j1,j2,1:3)               
             Se1f2(iedge,jface) = 0.5_SP*DOT_PRODUCT(V(i1,i2,1:3),temp2(1:3))
          ELSE
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Se1f1(iedge,jface)  = Se1f1(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) *& 
                     DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(3,lambda,V,j1,j2,j3))
                Se1f2(iedge,jface)  = Se1f2(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) *& 
                     DOT_PRODUCT(CURL_VBF(1,lambda,V,i1,i2),CURL_VBF(4,lambda,V,j1,j2,j3))
             END DO
          END IF
          SELECT CASE(ELEMENT_TYPE) 
          CASE(1,3) ! S&P edge-based type 2 element contributions are zero by definition.
             ! As are Webb99 (e2 elements identical).
             Se2f1(iedge,jface) = 0.0_SP
             Se2f2(iedge,jface) = 0.0_SP
          CASE(2) ! Cubature only
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Se2f1(iedge,jface)  = Se2f1(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) *& 
                     DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(3,lambda,V,j1,j2,j3))
                Se2f2(iedge,jface)  = Se2f2(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) *& 
                     DOT_PRODUCT(CURL_VBF(2,lambda,V,i1,i2),CURL_VBF(4,lambda,V,j1,j2,j3))
             END DO
          CASE DEFAULT
             STOP 'IE: Invalid element type in S_AND_T_MAKE_HIERARCHAL'
          END SELECT
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             ! Scale all - even if zero.
             Se1f1(iedge,jface) = ell(iedge)*Se1f1(iedge,jface)
             Se1f2(iedge,jface) = ell(iedge)*Se1f2(iedge,jface)
             Se2f1(iedge,jface) = ell(iedge)*Se2f1(iedge,jface)
             Se2f2(iedge,jface) = ell(iedge)*Se2f2(iedge,jface)
          ENDIF
       END DO FACELOOP
    END DO EDGELOOP
    Sf1e1 = TRANSPOSE(Se1f1)
    Sf1e2 = TRANSPOSE(Se2f1)              
    Sf2e1 = TRANSPOSE(Se1f2)
    Sf2e2 = TRANSPOSE(Se2f2)          
    Se(1:6,13:16)  = Se1f1(1:6,1:4)
    Se(1:6,17:20)  = Se1f2(1:6,1:4)     
    Se(7:12,13:16) = Se2f1(1:6,1:4)
    Se(7:12,17:20) = Se2f2(1:6,1:4)
    Se(13:16,1:6)  = Sf1e1(1:4,1:6)
    Se(17:20,1:6)  = Sf2e1(1:4,1:6)     
    Se(13:16,7:12) = Sf1e2(1:4,1:6)
    Se(17:20,7:12) = Sf2e2(1:4,1:6)
  END SUBROUTINE S_EDGE_FACE_ORDER2_ANALYTIC


  SUBROUTINE S_FACE_FACE_ORDER2_ANALYTIC(Se, V, M)
    !*******************************************************************************
    ! Fills face-face S submatrices for the LT/QN elements and 
    ! stores in elemental S matrix. These computed are Sf1f1,Sf1f2,Sf2f2;
    ! Sf2f1 is filled using symmetry. 
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se
    REAL(SP), DIMENSION(4,4,3), INTENT(in) ::   V  ! See [eqn(7),S&P]
    REAL(SP), DIMENSION(4,4), INTENT(in) :: M      ! integration matrix M [eqn(15),S&P]!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(4,4) :: Sf1f1, Sf1f2, Sf2f1, Sf2f2
    INTEGER(I4B) :: iface,jface        ! Misc counters
    INTEGER(I4B) :: i1,i2,i3,j1,j2,j3
    INTEGER(I4B), DIMENSION(3) :: i_face_nodes, j_face_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 11 ! Use 4-th order rule.
    ! Initialize (needed for cubature)
    Sf1f1(1:4,1:4) = 0.0_SP
    Sf1f2(1:4,1:4) = 0.0_SP
    Sf2f2(1:4,1:4) = 0.0_SP

    FACELOOP1: DO iface = 1,4
       FACELOOP2: DO jface = 1,4 
          i_face_nodes = LOCAL_FACENODES(iface)
          i1 = i_face_nodes(1)
          i2 = i_face_nodes(2)
          i3 = i_face_nodes(3)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)
          IF (.NOT.CUBATURE) THEN
             Sf1f1(iface,jface) = & 
                  4.0_SP*M(i1,j1)* DOT_PRODUCT(V(i2,i3,:),V(j2,j3,:)) &
                  +2.0_SP*M(i1,j2)* DOT_PRODUCT(V(i2,i3,:),V(j1,j3,:)) & 
                  -2.0_SP*M(i1,j3)* DOT_PRODUCT(V(i2,i3,:),V(j1,j2,:)) & 
                  +2.0_SP*M(i2,j1)* DOT_PRODUCT(V(i1,i3,:),V(j2,j3,:)) & 
                  +       M(i2,j2)* DOT_PRODUCT(V(i1,i3,:),V(j1,j3,:)) & 
                  -       M(i2,j3)* DOT_PRODUCT(V(i1,i3,:),V(j1,j2,:)) & 
                  -2.0_SP*M(i3,j1)* DOT_PRODUCT(V(i1,i2,:),V(j2,j3,:)) & 
                  -       M(i3,j2)* DOT_PRODUCT(V(i1,i2,:),V(j1,j3,:)) &   
                  +       M(i3,j3)* DOT_PRODUCT(V(i1,i2,:),V(j1,j2,:))        
             Sf1f2(iface,jface) = & 
                  2.0_SP*M(i1,j1)* DOT_PRODUCT(V(i2,i3,:),V(j2,j3,:)) &
                  +4.0_SP*M(i1,j2)* DOT_PRODUCT(V(i2,i3,:),V(j1,j3,:)) & 
                  +2.0_SP*M(i1,j3)* DOT_PRODUCT(V(i2,i3,:),V(j1,j2,:)) & 
                  +       M(i2,j1)* DOT_PRODUCT(V(i1,i3,:),V(j2,j3,:)) & 
                  +2.0_SP*M(i2,j2)* DOT_PRODUCT(V(i1,i3,:),V(j1,j3,:)) & 
                  +       M(i2,j3)* DOT_PRODUCT(V(i1,i3,:),V(j1,j2,:)) & 
                  -       M(i3,j1)* DOT_PRODUCT(V(i1,i2,:),V(j2,j3,:)) & 
                  -2.0_SP*M(i3,j2)* DOT_PRODUCT(V(i1,i2,:),V(j1,j3,:)) &   
                  -       M(i3,j3)* DOT_PRODUCT(V(i1,i2,:),V(j1,j2,:))        
             Sf2f2(iface,jface) = & 
                  M(i1,j1)* DOT_PRODUCT(V(i2,i3,:),V(j2,j3,:)) &
                  +2.0_SP*M(i1,j2)* DOT_PRODUCT(V(i2,i3,:),V(j1,j3,:)) & 
                  +       M(i1,j3)* DOT_PRODUCT(V(i2,i3,:),V(j1,j2,:)) & 
                  +2.0_SP*M(i2,j1)* DOT_PRODUCT(V(i1,i3,:),V(j2,j3,:)) & 
                  +4.0_SP*M(i2,j2)* DOT_PRODUCT(V(i1,i3,:),V(j1,j3,:)) & 
                  +2.0_SP*M(i2,j3)* DOT_PRODUCT(V(i1,i3,:),V(j1,j2,:)) & 
                  +       M(i3,j1)* DOT_PRODUCT(V(i1,i2,:),V(j2,j3,:)) & 
                  +2.0_SP*M(i3,j2)* DOT_PRODUCT(V(i1,i2,:),V(j1,j3,:)) &   
                  +       M(i3,j3)* DOT_PRODUCT(V(i1,i2,:),V(j1,j2,:))        
          ELSE 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Sf1f1(iface,jface)  = Sf1f1(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * & 
                     DOT_PRODUCT(CURL_VBF(3,lambda,V,i1,i2,i3),CURL_VBF(3,lambda,V,j1,j2,j3))
                Sf1f2(iface,jface)  = Sf1f2(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * & 
                     DOT_PRODUCT(CURL_VBF(3,lambda,V,i1,i2,i3),CURL_VBF(4,lambda,V,j1,j2,j3))
                Sf2f2(iface,jface)  = Sf2f2(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * & 
                     DOT_PRODUCT(CURL_VBF(4,lambda,V,i1,i2,i3),CURL_VBF(4,lambda,V,j1,j2,j3))   
             END DO
          END IF
       END DO FACELOOP2
    END DO FACELOOP1
    Sf2f1 = TRANSPOSE(Sf1f2)
    Se(13:16,13:16)  = Sf1f1(1:4,1:4)    
    Se(13:16,17:20)  = Sf1f2(1:4,1:4)     
    Se(17:20,13:16)  = Sf2f1(1:4,1:4)
    Se(17:20,17:20)  = Sf2f2(1:4,1:4)     
  END SUBROUTINE S_FACE_FACE_ORDER2_ANALYTIC


  SUBROUTINE T_EDGE_EDGE_ORDER1_ANALYTIC(Te, grad_lambda, M, ell, phi)
    !*******************************************************************************
    ! Fills edge-edge T submatrices for the lowest order elements and 
    ! stores in elemental T matrix (F in [eqn(16), S&P]).  
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!  
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
    REAL(SP), DIMENSION(4,3), INTENT(in) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,4), INTENT(in) :: M      ! integration matrix M [eqn(15),S&P]  
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
    REAL(SP), DIMENSION(4,4), INTENT(in) ::   phi  ! See [eqn(8),S&P] 
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Te1e1
    INTEGER(I4B) :: iedge,jedge        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 4 ! 2nd order rule
    ! Initialize (needed for cubature)
    Te1e1 = 0.0_SP

    EDGE_LOOP1: DO iedge = 1,6
       EDGE_LOOP2: DO jedge = 1,iedge ! Fill lower half
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)
          IF(.NOT.CUBATURE) THEN
             Te1e1(iedge,jedge) = ( phi(i2,j2)*M(i1,j1) - phi(i2,j1)*M(i1,j2)  & 
                  - phi(i1,j2)*M(i2,j1) + phi(i1,j1)*M(i2,j2) )     
          ELSE
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Te1e1(iedge,jedge) = Te1e1(iedge,jedge) + &
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(1,lambda,grad_lambda,j1,j2))
             END DO
          END IF
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             Te1e1(iedge,jedge) = ell(iedge)*ell(jedge)*Te1e1(iedge,jedge)
          ENDIF
          Te1e1(jedge,iedge) = Te1e1(iedge,jedge) ! Symmetrize. Diagonal is overwritten  
       END DO EDGE_LOOP2
    END DO EDGE_LOOP1
    Te(1:6,1:6) = Te1e1(1:6,1:6) 
  END SUBROUTINE T_EDGE_EDGE_ORDER1_ANALYTIC

  SUBROUTINE T_EDGE_EDGE_ORDER2_ANALYTIC(Te, grad_lambda, M, ell, phi)
    !*******************************************************************************
    ! Fills additional edge-edge T submatrices for the LT/QN elements and 
    ! stores in elemental T matrix. Those computed are Te1e2 and Te2e2; Te2e1 is 
    ! filled using symmetry
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!  
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
    REAL(SP), DIMENSION(4,3), INTENT(in) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,4), INTENT(in) :: M      ! integration matrix M [eqn(15),S&P]  
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
    REAL(SP), DIMENSION(4,4), INTENT(in) ::   phi  ! See [eqn(8),S&P] 
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,6) :: Te1e2, Te2e1, Te2e2
    INTEGER(I4B) :: iedge,jedge        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes, j_edge_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 11 ! 4th order rule. A&V elements are quadratic, product thus 4-th order. 
    ! Initialize (needed for cubature)
    Te1e2 = 0.0_SP ! Intialize
    Te2e2 = 0.0_SP ! Intialize

    EDGE_LOOP1: DO iedge = 1,6
       EDGE_LOOP2: DO jedge = 1,6
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_edge_nodes = LOCAL_EDGENODES(jedge)
          j1 = j_edge_nodes(1)
          j2 = j_edge_nodes(2)
          IF (.NOT.CUBATURE) THEN
             Te1e2(iedge,jedge) = ( phi(i2,j2)*M(i1,j1) + phi(i2,j1)*M(i1,j2)  & 
                  - phi(i1,j2)*M(i2,j1) - phi(i1,j1)*M(i2,j2) )     
             Te2e2(iedge,jedge) = ( phi(i2,j2)*M(i1,j1) + phi(i2,j1)*M(i1,j2)  & 
                  + phi(i1,j2)*M(i2,j1) + phi(i1,j1)*M(i2,j2) )
          ELSE 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Te1e2(iedge,jedge)  = Te1e2(iedge,jedge) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(2,lambda,grad_lambda,j1,j2))
                Te2e2(iedge,jedge)  = Te2e2(iedge,jedge) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(2,lambda,grad_lambda,j1,j2))
             END DO
          END IF
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             Te1e2(iedge,jedge) = ell(iedge)*ell(jedge)*Te1e2(iedge,jedge)
             Te2e2(iedge,jedge) = ell(iedge)*ell(jedge)*Te2e2(iedge,jedge)      
          ENDIF
       END DO EDGE_LOOP2
    END DO EDGE_LOOP1
    Te2e1 = TRANSPOSE(Te1e2) 
    Te(7:12,1:6)  = Te2e1(1:6,1:6)              
    Te(1:6,7:12)  = Te1e2(1:6,1:6)       
    Te(7:12,7:12) = Te2e2(1:6,1:6)      
  END SUBROUTINE T_EDGE_EDGE_ORDER2_ANALYTIC


  SUBROUTINE T_EDGE_FACE_ORDER2_ANALYTIC(Te, grad_lambda, N, ell, phi)
    !*******************************************************************************
    ! Fills edge-face S submatrices for the LT/QN elements and 
    ! stores in elemental S matrix. These computed are Te1f1,Te1f2,Te2f1,Te2f2; the 
    ! others are filled using symmetry. 
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!  
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
    REAL(SP), DIMENSION(4,3), INTENT(in) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,4,4), INTENT(in) :: N    ! integration matrix N [eqn(46),S&P]
    REAL(SP), DIMENSION(6), INTENT(in) :: ell      ! edge lengths
    REAL(SP), DIMENSION(4,4), INTENT(in) ::   phi  ! See [eqn(8),S&P] 
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(6,4) :: Te1f1, Te1f2, Te2f1, Te2f2
    COMPLEX(SPC), DIMENSION(4,6) :: Tf1e1, Tf1e2, Tf2e1, Tf2e2
    INTEGER(I4B) :: iedge,jface        ! Misc counters
    INTEGER(I4B) :: i1,i2,j1,j2,j3
    INTEGER(I4B), DIMENSION(2) :: i_edge_nodes
    INTEGER(I4B), DIMENSION(3) :: j_face_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 11 ! 4th order rule. A&V elements are quadratic, product thus 4-th order. 
    ! Initialize (needed for cubature)
    Te1f1(1:6,1:4) = 0.0_SP
    Te1f2(1:6,1:4) = 0.0_SP
    Te2f1(1:6,1:4) = 0.0_SP
    Te2f2(1:6,1:4) = 0.0_SP

    EDGE_LOOP: DO iedge = 1,6
       FACE_LOOP: DO jface = 1,4 
          i_edge_nodes = LOCAL_EDGENODES(iedge)
          i1 = i_edge_nodes(1)
          i2 = i_edge_nodes(2)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)               
          IF (.NOT.CUBATURE) THEN
             Te1f1(iedge,jface) =                                           &
                  ( phi(i2,j3)*N(i1,j1,j2) - phi(i1,j3)*N(i2,j1,j2) & 
                  - phi(i2,j2)*N(i1,j1,j3) + phi(i1,j2)*N(i2,j1,j3) ) 
             Te2f1(iedge,jface) =                                           &
                  ( phi(i2,j3)*N(i1,j1,j2) + phi(i1,j3)*N(i2,j1,j2) & 
                  - phi(i2,j2)*N(i1,j1,j3) - phi(i1,j2)*N(i2,j1,j3) )   
             Te1f2(iedge,jface) =                                           &
                  ( phi(i2,j3)*N(i1,j1,j2) - phi(i1,j3)*N(i2,j1,j2) & 
                  - phi(i2,j1)*N(i1,j2,j3) + phi(i1,j1)*N(i2,j2,j3) ) 
             Te2f2(iedge,jface) =                                           &
                  ( phi(i2,j3)*N(i1,j1,j2) + phi(i1,j3)*N(i2,j1,j2) & 
                  - phi(i2,j1)*N(i1,j2,j3) - phi(i1,j1)*N(i2,j2,j3) ) 
          ELSE 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Te1f1(iedge,jface)  = Te1f1(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(3,lambda,grad_lambda,j1,j2,j3))
                Te2f1(iedge,jface)  = Te2f1(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(3,lambda,grad_lambda,j1,j2,j3))
                Te1f2(iedge,jface)  = Te1f2(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(1,lambda,grad_lambda,i1,i2),VBF(4,lambda,grad_lambda,j1,j2,j3))
                Te2f2(iedge,jface)  = Te2f2(iedge,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(2,lambda,grad_lambda,i1,i2),VBF(4,lambda,grad_lambda,j1,j2,j3))
             END DO
          END IF
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             Te1f1(iedge,jface) = ell(iedge)*Te1f1(iedge,jface)
             Te2f1(iedge,jface) = ell(iedge)*Te2f1(iedge,jface)
             Te1f2(iedge,jface) = ell(iedge)*Te1f2(iedge,jface)
             Te2f2(iedge,jface) = ell(iedge)*Te2f2(iedge,jface)
          ENDIF
       END DO FACE_LOOP
    END DO EDGE_LOOP
    Tf1e1 = TRANSPOSE(Te1f1)        
    Tf1e2 = TRANSPOSE(Te2f1)        
    Tf2e1 = TRANSPOSE(Te1f2)              
    Tf2e2 = TRANSPOSE(Te2f2)           
    Te(1:6,13:16)  = Te1f1(1:6,1:4)
    Te(1:6,17:20)  = Te1f2(1:6,1:4)
    Te(7:12,13:16) = Te2f1(1:6,1:4)
    Te(7:12,17:20) = Te2f2(1:6,1:4)
    Te(13:16,1:6)  = Tf1e1(1:4,1:6)
    Te(17:20,1:6)  = Tf2e1(1:4,1:6)
    Te(13:16,7:12) = Tf1e2(1:4,1:6)
    Te(17:20,7:12) = Tf2e2(1:4,1:6)
  END SUBROUTINE T_EDGE_FACE_ORDER2_ANALYTIC

  SUBROUTINE T_FACE_FACE_ORDER2_ANALYTIC(Te, grad_lambda, P, phi)
    !*******************************************************************************
    ! Fills face-face T submatrices for the LT/QN elements and 
    ! stores in elemental T matrix. These computed are Tf1f1,Tf1f2,Tf2f2;
    ! Tf2f1 is filled using symmetry. 
    !*******************************************************************************
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!  
    COMPLEX(SPC), INTENT(OUT), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Te
    REAL(SP), DIMENSION(4,3), INTENT(in) :: grad_lambda ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,4,4,4), INTENT(in) :: P  ! integration matrix P [eqn(47),S&P]
    REAL(SP), DIMENSION(4,4), INTENT(in) ::   phi  ! See [eqn(8),S&P] 
!!!
!!! Local Variables
!!!
    COMPLEX(SPC), DIMENSION(4,4) :: Tf1f1, Tf1f2, Tf2f1, Tf2f2
    INTEGER(I4B) :: iface,jface        ! Misc counters
    INTEGER(I4B) :: i1,i2,i3,j1,j2,j3
    INTEGER(I4B), DIMENSION(3) :: i_face_nodes, j_face_nodes
    INTEGER(I4B) num_cpts, cube_point  ! number and index of cubature points
    REAL(SP), DIMENSION(4) :: lambda   ! Simplex coordinates of point.

    num_cpts = 11 ! 4th order rule. A&V elements are quadratic, product thus 4-th order. 
    ! Initialize (needed for cubature)
    Tf1f1(1:4,1:4) = 0.0_SP
    Tf1f2(1:4,1:4) = 0.0_SP
    Tf2f2(1:4,1:4) = 0.0_SP

    FACE_LOOP1: DO iface = 1,4
       FACE_LOOP2: DO jface = 1,4 
          i_face_nodes = LOCAL_FACENODES(iface)
          i1 = i_face_nodes(1)
          i2 = i_face_nodes(2)
          i3 = i_face_nodes(3)
          j_face_nodes = LOCAL_FACENODES(jface)
          j1 = j_face_nodes(1)
          j2 = j_face_nodes(2)
          j3 = j_face_nodes(3)
          IF (.NOT.CUBATURE) THEN
             Tf1f1(iface,jface) =                                              &
                  ( phi(i3,j3)*P(i1,i2,j1,j2) - phi(i2,j3)*P(i1,i3,j1,j2) & 
                  - phi(i3,j2)*P(i1,i2,j1,j3) + phi(i2,j2)*P(i1,i3,j1,j3) ) 
             Tf1f2(iface,jface) =                                              &
                  ( phi(i3,j3)*P(i1,i2,j1,j2) - phi(i2,j3)*P(i1,i3,j1,j2) & 
                  - phi(i3,j1)*P(i1,i2,j2,j3) + phi(i2,j1)*P(i1,i3,j2,j3) ) 
             Tf2f2(iface,jface) =                                              &
                  ( phi(i3,j3)*P(i1,i2,j1,j2) - phi(i1,j3)*P(i2,i3,j1,j2) & 
                  - phi(i3,j1)*P(i1,i2,j2,j3) + phi(i1,j1)*P(i2,i3,j2,j3) )
          ELSE 
             DO cube_point=1,num_cpts
                lambda = cube_tet_rules(num_cpts)%rule(cube_point,1:4)
                Tf1f1(iface,jface)  = Tf1f1(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),VBF(3,lambda,grad_lambda,j1,j2,j3))
                Tf1f2(iface,jface)  = Tf1f2(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(3,lambda,grad_lambda,i1,i2,i3),VBF(4,lambda,grad_lambda,j1,j2,j3))
                Tf2f2(iface,jface)  = Tf2f2(iface,jface) + & 
                     cube_tet_rules(num_cpts)%rule(cube_point,5) * &
                     DOT_PRODUCT(VBF(4,lambda,grad_lambda,i1,i2,i3),VBF(4,lambda,grad_lambda,j1,j2,j3))
             END DO
          END IF
       END DO FACE_LOOP2
    END DO FACE_LOOP1
    Tf2f1 = TRANSPOSE(Tf1f2)
    Te(13:16,13:16)  = Tf1f1(1:4,1:4)     
    Te(13:16,17:20)  = Tf1f2(1:4,1:4)          
    Te(17:20,13:16)  = Tf2f1(1:4,1:4)
    Te(17:20,17:20)  = Tf2f2(1:4,1:4)     
  END SUBROUTINE T_FACE_FACE_ORDER2_ANALYTIC

  !*******************************************************************************
END MODULE S_and_T_matrix
