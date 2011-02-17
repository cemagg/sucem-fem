! Last changed:
! 20 Jul 04: Changes to B_MAKE_HIERARCHAL. DBD.
! 2002-05-09: Changes to EVALUATE_ELEMENTAL_FUNCTIONS. MMB.
! Last updated DBD 05 Mar 2002. 
! NABLA_VBF renamed CURL_VBF (more accurate name). 
! NB: new arguments added to B_MAKE_HIERARCHAL_ANAYTIC and S_AND_T_MAKE_HIERARCHAL. 
! Routines still needing updating:
! POINT_EVALUATE_FACE_FUNCTIONS (add Webb99 face elements) [Matthys!]

MODULE inv_curvi
  USE nrtype
  USE debugvar
  USE unit_numbers
  IMPLICIT NONE
!!$  PRIVATE
!!$  PUBLIC :: INVERSE_CURVILINEAR
  !*******************************************************************************
  REAL(SP) :: tol_xx,tol_func       ! Tolerances on vector xx and function func.
  REAL(SP), DIMENSION(10) :: xnod, ynod, znod ! JPS added
  REAL(SP), DIMENSION(3) :: FUNC, delta_xx
  REAL(SP), DIMENSION(3,3) :: JACOB_FUNC
  INTEGER(I4B) :: iter
  INTEGER(I4B), PARAMETER :: NDIM=3 ! Dimensionality of problem.	
  INTEGER(I4B), PARAMETER :: NTRIAL=20000 ! Number of Newton Raphson iterations.
  INTEGER(I4B) :: info               ! Info from LAPACK factorization.
  INTEGER(I4B) :: lwork              ! LAPACK dimensioning data
  PARAMETER(lwork=3)
  REAL(SP), DIMENSION(lwork) :: temp_work ! LAPACK temporary workspace.
  INTEGER(I4B), DIMENSION(3) :: temp_ipiv ! Pivot indices from factorization. 

  ! JPS added declarations
  INTEGER, PARAMETER :: mpn=3
  ! INTEGER, PARAMETER :: mpldfjac=3
  INTEGER :: mpinfo
  INTEGER, PARAMETER :: mplwa=(mpn*(3*mpn+13))/2
  DOUBLE PRECISION, PARAMETER :: mptol = 1.0e-6
  DOUBLE PRECISION :: xmp(mpn),mpfvec(mpn),mpwa(mplwa)
  !EXTERNAL MAPPING_FUNC
  ! end jps added declarations
  REAL(SP), DIMENSION(3) :: r_vec_m ! Coordinates of point in real coordinate system..

CONTAINS

  ! JPS moved subroutine INVERSE_CURVILINEAR here from outside subroutine evalueate_elemental_functions
  SUBROUTINE INVERSE_CURVILINEAR(elem, r_vec,u,v,w,point_found)
    USE geometry, ONLY: avg_edge_length, MID_EDGE_NODES
    
    !*******************************************************************************
    ! Using a Newton-Raphson method, attempt to invert the 2nd order mapping and find the point u,v,w corresponding to 
    ! r_vec (i.e. x,y,z) within this element. 
    ! As a first guess, the point corresponding to the 
    ! rectilinear mapping is used - this is the value of u,v,w on entry. 
    ! Note that it is possible, especially for a point lying close to the element faces in the rectilinear element, 
    ! that under a curvilinear mapping, the point may now lie within another element. 
    ! In this case, the routine terminates after a specified number of iterations, returning 
    ! a failure flag to the calling routine.
    !
    ! The Newton-Raphson algorithm is described in "Numerical Recipes in Fortran", 
    ! Press et al, 2nd Edn 1992. This is, however, an independent implementation. 
    !
    ! CAUTION!! This algorithm has not been properly tested to date for elements with all sides curved.
    !
    ! DBD 03 August 2005.
    !
    ! Modified by JPS on 8 September 2005 
    ! to USE the minpack globally convergent nonlinear solver
    ! this solver implements a line search algorithm to overcome the convergence
    ! problems experienced by the standard Newton's method taking the full
    ! newton step.
    !!*******************************************************************************
    REAL(SP), DIMENSION(3), INTENT(IN) :: r_vec ! Coordinates of point in real coordinate system..
    INTEGER(I4B), INTENT(IN) :: elem
    REAL(SP), INTENT(out) :: u,v,w ! Coordinates of point in reference coordinate system (if found, 0 otherwise)
    LOGICAL(LGT), OPTIONAL, INTENT(out) :: point_found

    tol_xx = 0.001    ! u,v,w vary from 0 to 1, so this is effectively normalized.
    tol_func = tol_xx*avg_edge_length ! Normalized by average edge length in mesh.
    r_vec_m = r_vec             ! Initialise the module global variable
    point_found = .FALSE.

    ! added by JPS to output element nodes in debug mode
    CALL MID_EDGE_NODES(elem,xnod,ynod,znod)

    IF (DEBUG_NEWTON_RAPHSON) THEN
       WRITE(FILEOUT,'(/,1X,A)') 'Newton Raphson scheme'
       WRITE(FILEOUT, *) 'elem = ', elem
       WRITE(FILEOUT,*) 'nodal x-coords = ', xnod
       WRITE(FILEOUT,*) 'nodal y-coords = ', ynod
       WRITE(FILEOUT,*) 'nodal z-coords = ', znod
       WRITE(FILEOUT,*) '(x,y,z) = ',  r_vec(1),r_vec(2),r_vec(3)
       !  WRITE(FILEOUT,'(1X,A)') 'iteration u    v     w          function values'
    END IF
    !xx = (/u,v,w/)! Initialize vector. 

    ! JPS addition
    ! initialize xx
    xmp = (/ 0.0, 0.0, 0.0 /)
    ! call minpack solver
    CALL hybrd1(MAPPING_FUNC,mpn,xmp,mpfvec,mptol,mpinfo,mpwa,mplwa)
    u = xmp(1)
    v = xmp(2)
    w = xmp(3)
    IF (mpinfo == 1) THEN
       point_found = .TRUE.
       IF (DEBUG_NEWTON_RAPHSON) THEN
          WRITE(FILEOUT,'(1X,A,3(F6.4,2X))') 'Newton Raphson successful. u v w = ',xmp
       END IF
    ELSE
       WRITE(FILEOUT,*) 'Nonlinear system solve using minpack was unsuccessful in routine INVERSE_CURVILINEAR'
       WRITE(FILEOUT,*) 'u v w = ', xmp
       WRITE(*,*) 'JPS: minpack routine was unsuccessful with (u,v,w) = ', xmp
       CONTINUE ! Routine exits with u,v,w unchanged.
    END IF


  END SUBROUTINE INVERSE_CURVILINEAR

  ! nonlinear mapping function to be passed to minpack nonlinear solver
  SUBROUTINE MAPPING_FUNC(nn, xx, fvec, iflag)
    ! subroutine to specify test nonlinear equation system
    ! without it's jacobian
    ! written by Julian P. Swartz
    ! last modified on 6 September 2005

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nn, iflag
    DOUBLE PRECISION, INTENT(INOUT) :: xx(nn)
    DOUBLE PRECISION, INTENT(OUT) :: fvec(nn)
    DOUBLE PRECISION :: uu, vv, ww
    DOUBLE PRECISION :: N_e(10)
    INTEGER :: i

    ! initialize fvec to zero
    fvec = 0

    uu = xx(1)
    vv = xx(2)
    ww = xx(3)

    ! geometrical transformation functions
    N_e = (/(2.0*(1.0-uu-vv-ww)-1.0)*(1.0-uu-vv-ww), &
         (2.0*uu-1.0)*uu, &
         (2.0*vv-1.0)*vv, &
         (2.0*ww-1.0)*ww, &
         4.0*(1.0-uu-vv-ww)*uu, &
         4.0*(1.0-uu-vv-ww)*vv, &
         4.0*(1.0-uu-vv-ww)*ww, &
         4.0*uu*vv, &
         4.0*uu*ww, &
         4.0*vv*ww/)


    DO i = 1,10
       fvec(1) = fvec(1) + N_e(i)*xnod(i)
       fvec(2) = fvec(2) + N_e(i)*ynod(i)
       fvec(3) = fvec(3) + N_e(i)*znod(i)
    END DO

    fvec(1) = fvec(1) - r_vec_m(1)
    fvec(2) = fvec(2) - r_vec_m(2)

    fvec(3) = fvec(3) - r_vec_m(3)

  END SUBROUTINE MAPPING_FUNC

  ! end jps move
END MODULE inv_curvi


MODULE basis_function
  USE geometry
  USE math_tools, ONLY: CROSS_PRODUCT,VECTOR_LENGTH
  USE nrtype
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!!$  PRIVATE
!!$  PUBLIC :: CURL_CURL_VBF, &
!!$       CURL_VBF, &
!!$       DIV_VBF, &
!!$       ELEMENTAL_DOF_FLAGS, &
!!$       EVALUATE_ELEMENTAL_FUNCTIONS, &
!!$       REAL_GAUSS_QUAD_LINE_VBF, &
!!$       REAL_INTERPOLATE_FUNCTION, &
!!$       COMPLEX_INTERPOLATE_FUNCTION, &
!!$       COMPLEX_GAUSS_QUAD_LINE_VBF, &
!!$       POINT_EVALUATE_FACE_FUNCTIONS, &
!!$       VBF, &
!!$       VBF_S, &
!!$       CURL_PARENT_VBF


  !*******************************************************************************
  ! This module contains routines that explicitly state the basis functions.
  !
  ! Explicit basis functions definitions are contained in the following routines:
  !
  ! CURL_VBF
  ! POINT_EVALUATE_FACE_FUNCTIONS
  ! VBF
  ! VBF_S
  !
  ! Created MMB 2001 Sept 26 from pre-existing routines. 
  ! Extensive changes in S_AND_T_MAKE_HIERARCHAL, with some new 
  ! 2002-05-07: Added ELEMENTAL_DOF_FLAGS. MMB.
  ! internal subprograms. 
  ! Last changed 04 March 2002. DBD:
  !   Documentation (comments) updated. 
  !   Work started on polynomial complete basis functions;  
  !   presently implemented up to LT/QN, with dummy code for QT/QN elements 
  !   (presently returns zero for the higher-order terms.)
  !
  !*******************************************************************************
  SAVE
  !*******************************************************************************
  ! Internal routines:
  ! (A routine is either listed as an interface, or in the private statement.)

  !*******************************************************************************
CONTAINS
  !*******************************************************************************



  ! AddedMMB
  FUNCTION CURL_CURL_VBF(vbf_type,lambda,grad_lambda,node1,node2,node3)
    USE math_tools, ONLY: CROSS_PRODUCT
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! Analytical expression for the curl-curl of the vector basis functions.
    ! Note that scaling by the edge length (if appropriate) is done elsewhere. 
    ! Future extensions will include:
    ! 5 -> e3, 6 -> f3 (Up to and including QT/QN)
    ! "Dummy" code at present to return zeros.
    !
    ! 2002-03-06: Created. MMB
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::  vbf_type              ! Type of vector basis function:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda       ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda ! Gradient of simplex coordinates
    INTEGER(I4B), INTENT(IN) ::  node1,node2
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3       ! Only required for face-based functions.
    REAL(SP), DIMENSION (3) :: CURL_CURL_VBF

    REAL(SP), DIMENSION (3) :: temp_vec

    IF (vbf_type.GE.3.AND..NOT.PRESENT(node3)) THEN
       STOP 'IE: Internal error in FUNCTION CURL_CURL_VBF, third node required'
    END IF

    CURL_CURL_VBF = 0.0 ! initialise

    SELECT CASE(vbf_type) 

    CASE (1) ! e1
       CURL_CURL_VBF(1:3) = 0.0_SP ! For all element types.

    CASE (2) ! e2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,3) ! Savage98 and Webb99 types. 
          CURL_CURL_VBF(1:3) = 0.0_SP
       CASE (2)   ! Andersen & Volakis types. 
          temp_vec = CROSS_PRODUCT(grad_lambda(node1,1:3),grad_lambda(node2,1:3))
	  CURL_CURL_VBF(1:3) =   3.0*CROSS_PRODUCT(grad_lambda(node1,1:3),temp_vec) &
               - 3.0*CROSS_PRODUCT(grad_lambda(node2,1:3),temp_vec)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION NABLA_VBF'
       END SELECT

    CASE (3)! f1
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2) ! Savage98 and A&V types. 
          temp_vec = CROSS_PRODUCT(grad_lambda(node2,1:3),grad_lambda(node3,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               + 2.0*CROSS_PRODUCT(grad_lambda(node1,1:3),temp_vec)
          temp_vec = CROSS_PRODUCT(grad_lambda(node1,1:3),grad_lambda(node3,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               + CROSS_PRODUCT(grad_lambda(node2,1:3),temp_vec)
          temp_vec = CROSS_PRODUCT(grad_lambda(node1,1:3),grad_lambda(node2,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               - CROSS_PRODUCT(grad_lambda(node3,1:3),temp_vec)
       CASE (3) ! Webb
          STOP 'IE: not implemented in CURL_CURL_VBF.'
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION NABLA_VBF'
       END SELECT

    CASE (4) !f2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2) ! Savage98 and A&V types. 
          temp_vec = CROSS_PRODUCT(grad_lambda(node2,1:3),grad_lambda(node3,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               + CROSS_PRODUCT(grad_lambda(node1,1:3),temp_vec)
          temp_vec = CROSS_PRODUCT(grad_lambda(node1,1:3),grad_lambda(node3,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               + 2.0*CROSS_PRODUCT(grad_lambda(node2,1:3),temp_vec)
          temp_vec = CROSS_PRODUCT(grad_lambda(node2,1:3),grad_lambda(node1,1:3))
          CURL_CURL_VBF(1:3) = CURL_CURL_VBF(1:3) &
               - CROSS_PRODUCT(grad_lambda(node3,1:3),temp_vec)
       CASE (3) ! Webb
          STOP 'IE: not implemented in CURL_CURL_VBF.'
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION NABLA_VBF'
       END SELECT

    CASE (5) ! e3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements.
          CURL_CURL_VBF(1:3) = 0.0_SP
       CASE DEFAULT
          STOP 'IE: not implemented in CURL_CURL_VBF.'
       END SELECT

    CASE (6) ! f3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements.
          CURL_CURL_VBF(1:3) = 0.0_SP
       CASE DEFAULT
          STOP 'IE: not implemented in CURL_CURL_VBF.'
       END SELECT

    CASE DEFAULT
       STOP 'IE: Invalid vector basis function type'
    END SELECT

  END FUNCTION CURL_CURL_VBF
  !*******************************************************************************
  ! AddedMMB


  FUNCTION CURL_VBF(vbf_type,lambda,V,node1,node2,node3)
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! Analytical expression for the curl of the vector basis functions
    ! Note that scaling by the edge length (if appropriate) is done elsewhere. 
    ! Extended 24 Feb 02 to include Webb99 LT/QN elements.
    ! Ditto    05 Mar 02 to include Webb99 QT/QN elements.
    !          Routine re-named CURL_VBF (was originally NABLA_VBF, name was misleading)
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::  vbf_type 
    ! Type of vector basis function:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2, 5 -> e3, 6 -> f3
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
    ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4,4,3), INTENT(IN) ::   V  
    ! See [eqn(7),S&P]
    INTEGER(I4B), INTENT(IN) ::  node1,node2
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 ! Only required for face-based
    ! functions.
    REAL(SP), DIMENSION (3) :: CURL_VBF

    IF (vbf_type.GE.3.AND..NOT.PRESENT(node3)) THEN
       STOP 'IE: Internal error in FUNCTION CURL_VBF, third node required'
    END IF

    SELECT CASE(vbf_type) 
    CASE (1) ! e1
       CURL_VBF(1:3) = 2.0_SP*V(node1,node2,1:3) ! For all element types.
    CASE (2) ! e2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,3) ! Savage98 and Webb99 types. 
          CURL_VBF(1:3) = 0.0_SP
       CASE (2)   ! Andersen & Volakis types. 
          CURL_VBF(1:3) = & 
               - lambda(node1) * V(node2,node2,1:3) & 
               + (3.0_SP*lambda(node1)-3.0_SP*lambda(node2)) * V(node1,node2,1:3) & 
               - lambda(node2) * V(node1,node1,1:3) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION CURL_VBF'
       END SELECT
    CASE (3)! f1
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2) ! Savage98 and A&V types. 
          CURL_VBF(1:3) = 2.0_SP*lambda(node1) * V(node2,node3,1:3) & 
               + lambda(node2) * V(node1,node3,1:3) & 
               - lambda(node3) * V(node1,node2,1:3) 
       CASE (3) ! Webb
          CURL_VBF(1:3) = -3.0_SP*lambda(node2) * V(node1,node3,1:3) &
               -3.0_SP*lambda(node1) * V(node2,node3,1:3)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION CURL_VBF'
       END SELECT
    CASE (4) !f2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2) ! Savage98 and A&V types. 
          CURL_VBF(1:3) =        lambda(node1) * V(node2,node3,1:3) & 
               + 2.0_SP*lambda(node2) * V(node1,node3,1:3) & 
               + lambda(node3) * V(node1,node2,1:3) 
       CASE (3) ! Webb
          CURL_VBF(1:3) = 3.0_SP*lambda(node3) * V(node1,node2,1:3) &
               +3.0_SP*lambda(node2) * V(node1,node3,1:3)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION CURL_VBF'
       END SELECT
    CASE (5) ! e3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only. These elements are in the gradient space
          ! and have no curl by definition.
          CURL_VBF(1:3) = 0.0_SP
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION CURL_VBF'
       END SELECT
    CASE (6) ! f3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only. Ditto e3 elements. 
          CURL_VBF(1:3) = 0.0_SP
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION CURL_VBF'
       END SELECT
    CASE DEFAULT
       STOP 'IE: Invalid vector basis function type'
    END SELECT
  END FUNCTION CURL_VBF
  !*******************************************************************************


  FUNCTION CURL_PARENT_VBF(vbf_type,ii,u,v,w)
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! Analytical expression for the curl of the vector basis functions
    ! on a parent "unitary" tetrahedron. 
    ! Edge length scaling is not supported.
    ! DBD, 20 July 05. Implemented only for Webb99 basis functions. 
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::  vbf_type 
    ! Type of vector basis function:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2, 5 -> e3, 6 -> f3
    INTEGER(I4B), INTENT(IN) :: ii ! edge or face number 
    REAL(SP),INTENT(IN), OPTIONAL :: u,v,w  ! Local coordinates on parent triangle (associated with u1,u2,u3 or x,y,z or xi,eta,zeta etc)  
    REAL(SP), DIMENSION (3) :: CURL_PARENT_VBF

    IF (vbf_type.GE.3.AND..NOT.PRESENT(u).AND..NOT.PRESENT(v).AND..NOT.PRESENT(w)) THEN
       STOP 'IE: Internal error in FUNCTION CURL_PARENT_VBF, local coordinates must be provided for this type of basis function'
    END IF


    SELECT CASE(vbf_type) 
    CASE (1) ! e1
       IF (SCALE_BY_EDGE_LENGTH) THEN
          STOP 'IE: FUNCTION CURL_PARENT_VBF does not support this element type with edge length scaling'
       END IF
       SELECT CASE(ii) ! edge number
       CASE (1) 
          CURL_PARENT_VBF = (/  0.0_SP, -2.0_SP,  2.0_SP/)
       CASE (2) 
          CURL_PARENT_VBF = (/  2.0_SP,  0.0_SP, -2.0_SP/)
       CASE (3) 
          CURL_PARENT_VBF = (/ -2.0_SP,  2.0_SP,  0.0_SP/)
       CASE (4) 
          CURL_PARENT_VBF = (/  0.0_SP,  0.0_SP,  2.0_SP/)
       CASE (5) 
          CURL_PARENT_VBF = (/  0.0_SP, -2.0_SP,  0.0_SP/)
       CASE (6) 
          CURL_PARENT_VBF = (/  2.0_SP,  0.0_SP,  0.0_SP/)
       CASE DEFAULT
          STOP 'IE: Invalid edge number in FUNCTION CURL_PARENT_VBF'
       END SELECT
    CASE (2) ! e2
       CURL_PARENT_VBF = (/  0.0_SP, 0.0_SP,  0.0_SP/)
    CASE (3)! f1
       IF(ELEMENT_TYPE.NE.3) THEN
          STOP 'IE: Only Webb99-type basis functions supported in FUNCTION CURL_PARENT_VBF'
       END IF
       SELECT CASE(ii)
       CASE (1) 
          CURL_PARENT_VBF = (/                          -3.0_SP*u,                             0.0_SP, &
               -3.0_SP+6.0_SP*u+3.0_SP*v+3.0_SP*w/)
       CASE (2) 
          CURL_PARENT_VBF = (/                           3.0_SP*u,  3.0_SP-6.0_SP*u-3.0_SP*v-3.0_SP*w, &
               0.0_SP/)
       CASE (3) 
          CURL_PARENT_VBF = (/ -3.0_SP+3.0_SP*u+6.0_SP*v+3.0_SP*w,                          -3.0_SP*v, &
               0.0_SP/)
       CASE (4) 
          CURL_PARENT_VBF = (/                          -3.0_SP*u,                          +3.0_SP*v, &
               0.0_SP/)
       CASE DEFAULT
          STOP 'IE: Invalid face number in FUNCTION CURL_PARENT_VBF'
       END SELECT
    CASE (4)! f2
       IF(ELEMENT_TYPE.NE.3) THEN
          STOP 'IE: Only Webb99-type basis functions supported in FUNCTION CURL_PARENT_VBF'
       END IF
       SELECT CASE(ii)
       CASE (1) 
          CURL_PARENT_VBF = (/                           3.0_SP*u,                          -3.0_SP*v, &
               -3.0_SP*u+3.0_SP*v/)
       CASE (2) 
          CURL_PARENT_VBF = (/                          -3.0_SP*u,                  3.0_SP*u-3.0_SP*w, &
               3.0_SP*w/)
       CASE (3) 
          CURL_PARENT_VBF = (/                 -3.0_SP*v+3.0_SP*w,                          +3.0_SP*v, &
               -3.0_SP*w/)
       CASE (4) 
          CURL_PARENT_VBF = (/                             0.0_SP,                          -3.0_SP*v, &
               3.0_SP*w/)
       CASE DEFAULT
          STOP 'IE: Invalid face number in FUNCTION CURL_PARENT_VBF'
       END SELECT
    CASE (5) ! e3
       CURL_PARENT_VBF = (/  0.0_SP, 0.0_SP,  0.0_SP/)
    CASE (6) ! f3
       CURL_PARENT_VBF = (/  0.0_SP, 0.0_SP,  0.0_SP/)
    CASE DEFAULT
       STOP 'IE: Invalid or unimplemented vector basis function type in FUNCTION CURL_PARENT_VBF'
    END SELECT
  END FUNCTION CURL_PARENT_VBF
  !*******************************************************************************



  FUNCTION DIV_VBF(vbf_type,lambda,grad_lambda,node1,node2,node3)
    USE nrtype
    USE problem_info
    IMPLICIT NONE

    !*******************************************************************************
    ! Analytical expression for the divergence of the vector basis functions
    ! Note that scaling by the edge length (if appropriate) is done elsewhere. 
    ! Written DBD 29 Sep 04, based on CURL_VBF. Only supports Webb99 elements. 
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::  vbf_type 
    ! Type of vector basis function:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2, 5 -> e3, 6 -> f3
    ! (presently only supported up to LT/QN, ie 5 and 6 unimplented).
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
    ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda 
    ! Gradient of simplex coordinates
    INTEGER(I4B), INTENT(IN) ::  node1,node2
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 ! Only required for face-based
    ! functions.
    REAL(SP) DIV_VBF


    IF (vbf_type.GE.3.AND..NOT.PRESENT(node3)) THEN
       STOP 'IE: Internal error in FUNCTION CURL_VBF, third node required'
    END IF

    SELECT CASE(vbf_type) 
    CASE (1) ! e1
       DIV_VBF = 0.0 ! For all element types. Identity for Whitney elements.
    CASE (2) ! e2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3) ! Webb99 types. 
          DIV_VBF = 2.0_SP*DOT_PRODUCT(grad_lambda(node1,:),grad_lambda(node2,:))
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION DIV_VBF'
       END SELECT
    CASE (3)! f1
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3) ! Webb
          DIV_VBF = - DOT_PRODUCT(grad_lambda(node2,:),VBF((1),lambda,grad_lambda,node1,node3)) & 
               - DOT_PRODUCT(grad_lambda(node1,:),VBF((1),lambda,grad_lambda,node2,node3)) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION DIV_VBF'
       END SELECT
    CASE (4) !f2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3) ! Webb
          DIV_VBF = - DOT_PRODUCT(grad_lambda(node3,:),VBF((1),lambda,grad_lambda,node2,node1)) & 
               - DOT_PRODUCT(grad_lambda(node2,:),VBF((1),lambda,grad_lambda,node3,node1)) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION DIV_VBF'
       END SELECT
    CASE (5) ! e3
       STOP 'IE: Unimplemented element type in FUNCTION DIV_VBF'
    CASE (6) ! f3
       STOP 'IE: Unimplemented element type in FUNCTION DIV_VBF'
    CASE DEFAULT
       STOP 'IE: Invalid vector basis function type'
    END SELECT
  END FUNCTION DIV_VBF
  !*******************************************************************************


  SUBROUTINE ELEMENTAL_DOF_FLAGS(elem,dof_flags)
    USE geometry
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! This routine returns the array <dof_flags> with 0 at positions where the local
    ! basis function of element <elem> does not correspond to a d.o.f. and 1 otherwise.
    ! The length of the vector is the maximum number of basis functions that the code 
    ! has to offer. Thus the routine indicate which basis functions would correspond to
    ! dof's if the element was of the highest possible order.  This routine assumes
    ! that the <faces%free> and <edges%free> arrays are already assigned.
    !
    ! 2002-05-07: Created. MMB.
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: elem
    INTEGER(I4B), DIMENSION(ELEM_TET_MATRIX_SIZE), INTENT(OUT) :: dof_flags

    INTEGER(I4B) :: icount, global_num

    dof_flags = 0 ! initialise

    ! Cycle through the edges:
    DO icount = 1,6
       global_num = elements(elem)%edges(icount)
       IF (edges(global_num)%free) THEN
          dof_flags(     icount) = 1 ! e1 basis functions
          dof_flags(6 +  icount) = 1 ! e2 basis functions
          dof_flags(20 + icount) = 1 ! e3 basis functions
       END IF
    END DO

    ! Cycle through the faces:
    DO icount = 1,4
       global_num = elements(elem)%faces(icount)
       IF (faces(global_num)%free) THEN
          dof_flags(12 + icount) = 1 ! f1 basis functions
          dof_flags(16 + icount) = 1 ! f2 basis functions
          dof_flags(26 + icount) = 1 ! f3 basis functions
       END IF
    END DO

  END SUBROUTINE ELEMENTAL_DOF_FLAGS
  !*******************************************************************************


  SUBROUTINE EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec,curl_power,func_values,order_request,mixed_request)
    USE geometry
    USE math_tools, ONLY: CROSS_PRODUCT
    USE problem_info
    use inv_curvi
    IMPLICIT NONE
    !*******************************************************************************
    ! Returns the values of the elemental (= vector basis) functions 
    ! in element <elem>, at the xyz position <r_vec>.
    !
    ! MMB - 19 Dec 2000
    ! Re-written DBD 16 Feb 2001, to call function VBF rather than explicitly
    ! define the VBF's in this routine. 
    ! Extended DBD 04 Mar 2002 for up to QT/QN elements.
    ! 2002-03-07: Added the argument <curl_power> = 0/1/2, to indicate whether 
    !             (curl)^0, (curl)^1 or (curl)^2 of the basis functions must be 
    !             returned. MMB.
    ! 2002-05-09: Added the capability to explicitly request the order of the
    !             result. MMB.
    ! 2004-09-29: Added the option to compute the divergence of an element, by setting
    !             <curl_power> = -1. In this case, what is returned is the scalar value 
    !             of the divergence in func_values(1), with the other two components zero. DBD.
    ! 2005-08-03: Added limited capability to compute the elemental functions in 
    !             curvilinear elements. DBD.
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: elem
    REAL(SP), DIMENSION(3), INTENT(IN) :: r_vec
    INTEGER(I4B), INTENT(IN) :: curl_power
    REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,3), INTENT(OUT) :: func_values
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: order_request
    LOGICAL(LGT), OPTIONAL, INTENT(IN) :: mixed_request

    INTEGER(I4B) :: count,vcount1,vcount2         ! counters
    REAL(SP), DIMENSION(4,4) :: vert_mat,sim_mat  ! Matrix in eqn.(6) of S+P 1996
    REAL(SP), DIMENSION(4) :: simp_c,cart_c       ! Store coordinates
    INTEGER(I4B), DIMENSION(2) :: tempnodes
    INTEGER(I4B), DIMENSION(3) :: tempfacenodes
    REAL(SP), DIMENSION(6) :: edge_lengths
    REAL(SP), DIMENSION(4,3) :: grad_lam
    REAL(SP), DIMENSION(4,4,3) :: V_matrix        ! See [eqn(7),S&P]
    INTEGER(I4B) :: order_value                   ! Values which will be used to
    LOGICAL(LGT) :: mixed_value                   ! decide the function to be evaluated
    ! Added DBD 3 Aug 05
    INTEGER(I4B) :: iface                         ! more counters
    LOGICAL(LGT) :: curvilinear_flag              ! Flag determining whether element is curvilinear or not.
    REAL(SP) :: uu, vv, ww                        ! Coordinates of point in reference coordinate system
    LOGICAL(LGT) :: point_found                   ! Flag indicating whether above point was found.
    ! End added DBD 3 Aug 05


    IF (ELEMENT_TYPE.NE.1.AND.SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: EVALUATE_ELEMENTAL_FUNCTIONS does not support this element type with edge length scaling'
    END IF

    ! Esteblish the order & mixed values to which the functions must be evaluated:
    IF (PRESENT(order_request)) THEN
       order_value = order_request
    ELSE
       order_value = MAX_ORDER(elem)
    END IF
    IF (PRESENT(mixed_request)) THEN
       mixed_value = mixed_request
    ELSE
       mixed_value = MIXED_ORDER(elem)
    END IF

    ! Determine if the element is curvilinear or rectilinear (or if test mode with all curvilinear has been requested).
    curvilinear_flag=.FALSE. ! Default
    DO iface = 1,4
       IF (faces(elements(elem)%faces(iface))%curvilinear.OR.ALL_CURVILINEAR) THEN
          WRITE(*,*) "found a curvilinear face in EVALUATE_ELEMENTAL_FUNCTIONS"
          curvilinear_flag =.TRUE. 
          EXIT ! Found a face which is curvilinear, hence element is curvilinear. Exit loop.
       END IF
    END DO


    IF(curl_power.NE.0.AND.curvilinear_flag) THEN!
       STOP 'IE. ROUTINE EVALUATE_ELEMENTAL_FUNCTIONS CAN ONLY COMPUTE BASIS FUNCTIONS AT PRESENT, NOT CURL OR DIV THEREOF.'
    END IF

    ! Set all values to zero (This is a precaution to ensure the lower order elements
    ! do not return non-zero values for non-existent higher order terms.)
    func_values = 0.0_SP

    ! Calculate the simplex co-ordinates:
    CALL SIMPLEX_COEFFICIENTS(elem,sim_mat(1:4,4),sim_mat(1:4,1), &
         sim_mat(1:4,2),sim_mat(1:4,3),vert_mat)
    cart_c(1:3) = r_vec
    cart_c(4) = 1.0
    simp_c = MATMUL(sim_mat,cart_c)

    ! If curvilinear, try to find point uvw corresponding to r_vec
    IF (curvilinear_flag) THEN
       uu = simp_c(2)
       vv = simp_c(3)
       ww = simp_c(4)
       !       WRITE(*,*) "Calling subroutine INVERSE_CURVILINEAR"
       CALL INVERSE_CURVILINEAR(elem, r_vec,uu,vv,ww,point_found)
    END IF

    ! Calculate the 6 edge lengths for repeated later use:
    edge_lengths(1:6) = T_LENGTHS(elem)

    ! Calculate GRADIENT_LAMBDA for repeated later use:
    grad_lam = GRADIENT_LAMBDA(elem,.true.) 
    
    
    ! Check that <curl_power> is in the valid range, and calculate the
    ! neccesary V-matrix argument of CURL_VBF if <curl_power=1>:
    SELECT CASE(curl_power)
    CASE (-1:0)
       CONTINUE
    CASE (1)
       DO vcount1 = 1,4
          DO vcount2 = 1,4
             V_matrix(vcount1,vcount2,1:3) = &
		  CROSS_PRODUCT(grad_lam(vcount1,1:3),grad_lam(vcount2,1:3))
          END DO
       END DO
    CASE (2)
       CONTINUE
    CASE DEFAULT
       STOP 'IE: <curl_power> invalid in EVALUATE_ELEMENTAL_FUNCTIONS.'
    END SELECT

    ! CT/LN and higher order terms: 
    ! Evaluate the E1 functions.
    DO count = 1,6
       tempnodes = LOCAL_EDGENODES(count)
       SELECT CASE(curl_power)
       CASE (-1)
          func_values(count,1) = DIV_VBF(1,simp_c,grad_lam,tempnodes(1),tempnodes(2))   
       CASE (0)
	  IF (curvilinear_flag) THEN
             func_values(count,1:3) = PARENT_VBF(1,count,uu,vv,ww)
          ELSE
             func_values(count,1:3) = VBF(1,simp_c,grad_lam,tempnodes(1),tempnodes(2))   
          END IF
       CASE (1)
          func_values(count,1:3) = CURL_VBF(1,simp_c,V_matrix,tempnodes(1),tempnodes(2))
       CASE (2)
          func_values(count,1:3) = CURL_CURL_VBF(1,simp_c,grad_lam,tempnodes(1),tempnodes(2))
       END SELECT
       IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
          func_values(count,1:3) = edge_lengths(count)*func_values(count,1:3)
       END IF
    END DO

    ! LT/LN and higher order terms.
    IF (order_value.GE.2.OR.&
         .NOT.mixed_value.AND.order_value.GE.1) THEN
       ! Evaluate the E2 functions:
       DO count = 1,6
          tempnodes = LOCAL_EDGENODES(count)
          SELECT CASE(curl_power)
          CASE (-1)
             func_values(count+6,1) = DIV_VBF(2,simp_c,grad_lam,tempnodes(1),tempnodes(2))
          CASE (0)
             IF (curvilinear_flag) THEN
                func_values(count+6,1:3) = PARENT_VBF(2,count,uu,vv,ww)
             ELSE 
                func_values(count+6,1:3) = VBF(2,simp_c,grad_lam,tempnodes(1),tempnodes(2))
             END IF
          CASE (1)
             func_values(count+6,1:3) = CURL_VBF(2,simp_c,V_matrix,tempnodes(1),tempnodes(2))
          CASE (2)
             func_values(count+6,1:3) = CURL_CURL_VBF(2,simp_c,grad_lam,tempnodes(1),tempnodes(2))
          END SELECT
          IF (ELEMENT_TYPE.EQ.1.AND.SCALE_BY_EDGE_LENGTH) THEN
             func_values(count+6,1:3) = edge_lengths(count)*func_values(count+6,1:3)
          END IF
       END DO
    END IF

    ! LT/QN and higher order terms.
    IF (order_value.GE.2) THEN
       ! Evaluate the F1 functions:
       DO count = 1,4
          tempfacenodes = LOCAL_FACENODES(count)
          SELECT CASE(curl_power)
          CASE (-1)
             func_values(count+12,1) = DIV_VBF(3,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (0)
             IF (curvilinear_flag) THEN
                func_values(count+12,1:3) = PARENT_VBF(3,count,uu,vv,ww)
             ELSE 
                func_values(count+12,1:3) = VBF(3,simp_c,grad_lam,tempfacenodes(1), & 
                     tempfacenodes(2),tempfacenodes(3))
             END IF
          CASE (1)
             func_values(count+12,1:3) = CURL_VBF(3,simp_c,V_matrix,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (2)
             func_values(count+12,1:3) = CURL_CURL_VBF(3,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          END SELECT
       END DO
       ! Evaluate the F2 functions:
       DO count = 1,4
          tempfacenodes = LOCAL_FACENODES(count)
          SELECT CASE(curl_power)
          CASE (-1)
             func_values(count+16,1) = DIV_VBF(4,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (0)
             IF (curvilinear_flag) THEN
                func_values(count+16,1:3) = PARENT_VBF(4,count,uu,vv,ww)
             ELSE 
                func_values(count+16,1:3) = VBF(4,simp_c,grad_lam,tempfacenodes(1), & 
                     tempfacenodes(2),tempfacenodes(3))
             END IF
          CASE (1)
             func_values(count+16,1:3) = CURL_VBF(4,simp_c,V_matrix,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (2)
             func_values(count+16,1:3) = CURL_CURL_VBF(4,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          END SELECT
       END DO
    END IF

    ! QT/QN and higher order terms.
    IF (order_value.GE.3.OR.&
         .NOT.mixed_value.AND.order_value.GE.2) THEN
       ! Evaluate the three E3 functions:
       IF (curvilinear_flag) THEN
	  STOP 'IE: Unimplemented hierarchal order for curvilinear elements in EVALUATE_ELEMENTAL_FUNCTIONS.'
       END IF
       DO count = 1,6
          tempnodes = LOCAL_EDGENODES(count)
          SELECT CASE(curl_power)
          CASE (-1)
             STOP 'Unimplemented code in routine EVALUATE_ELEMENTAL_FUNCTIONS'
          CASE (0)
             func_values(count+20,1:3) = VBF(5,simp_c,grad_lam,tempnodes(1),tempnodes(2))   
          CASE (1)
             func_values(count+20,1:3) = CURL_VBF(5,simp_c,V_matrix,tempnodes(1),tempnodes(2))
          CASE (2)
             func_values(count+20,1:3) = CURL_CURL_VBF(5,simp_c,grad_lam,tempnodes(1),tempnodes(2))
          END SELECT
       END DO
       ! Evaluate the F2 functions:
       DO count = 1,4
          tempfacenodes = LOCAL_FACENODES(count)
          SELECT CASE(curl_power)
          CASE (-1)
             STOP 'Unimplemented code in routine EVALUATE_ELEMENTAL_FUNCTIONS'
          CASE (0)
             func_values(count+26,1:3) = VBF(6,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (1)
             func_values(count+26,1:3) = CURL_VBF(6,simp_c,V_matrix,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          CASE (2)
             func_values(count+26,1:3) = CURL_CURL_VBF(6,simp_c,grad_lam,tempfacenodes(1), & 
                  tempfacenodes(2),tempfacenodes(3))
          END SELECT
       END DO
    END IF

    IF (order_value.GE.4.OR..NOT.mixed_value.AND.order_value.GE.3) THEN
       STOP 'IE: Unimplemented hierarchal order in EVALUATE_ELEMENTAL_FUNCTIONS.'
    END IF

    IF (curvilinear_flag) THEN ! Transform functions from (u,v,w) to (x,y,z).
       DO count = 1,ELEM_TET_MATRIX_SIZE
          func_values(count,1:3) = MATMUL(JACOBIAN_POLY2(elem,uu,vv,ww,.TRUE.),func_values(count,1:3))
       END DO
    END IF


  CONTAINS

    FUNCTION POLY2(xx)
      !*******************************************************************************
      ! This function returns x,y,z (in vector POLY2) given u,v,w (in vector xx) for a 2nd order polynomial curvilinear mapping.
      !
      ! DBD 03 August 2005.
      !*******************************************************************************

      REAL(SP), DIMENSION(3) :: POLY2 
      REAL(SP), DIMENSION(3), INTENT(IN) :: xx ! Vector of coordinates of point in reference coordinate system 
      REAL(SP)  :: u,v,w                      ! Vector of coordinates of point in reference coordinate system 
      REAL(SP), DIMENSION(10) :: N_e
      INTEGER(I4B) inode
      REAL(SP), DIMENSION(10) ::  x, y, z        ! Nodal cofficients.

      u = xx(1);
      v = xx(2); 
      w = xx(3);

      CALL MID_EDGE_NODES(elem,x,y,z)

      N_e = (/ (2.0_SP*(1.0_SP-u-v-w)-1.0_SP)*(1.0_SP-u-v-w), & 
           (2.0_SP*u-1.0_SP)*u, & 
           (2.0_SP*v-1.0_SP)*v, & 
           (2.0_SP*w-1.0_SP)*w, & 
           4.0_SP*(1.0_SP-u-v-w)*u, & 
           4.0_SP*(1.0_SP-u-v-w)*v, & 
           4.0_SP*(1.0_SP-u-v-w)*w, & 
           4.0_SP*u*v, & 
           4.0_SP*u*w, & 
           4.0_SP*v*w  /)
      POLY2 = 0.0 ! Initialize on each call.
      DO inode = 1,10
         POLY2(1)= POLY2(1) +N_e(inode)*x(inode)
         POLY2(2)= POLY2(2) +N_e(inode)*y(inode)
         POLY2(3)= POLY2(3) +N_e(inode)*z(inode)
      END DO

    END FUNCTION POLY2

  END SUBROUTINE EVALUATE_ELEMENTAL_FUNCTIONS





  FUNCTION REAL_GAUSS_QUAD_LINE_VBF(element_num,local_edge_num,&
       VBF_type,time_derivative_order,num_qpoints,ell)

    USE boundary_conditions
    USE geometry
    USE problem_info
    USE quad_tables
    USE scattering_analysis_data
    USE TD_source, ONLY: TD_INC_FIELD,TD_D_BY_DT_INC_FIELD,TD_D2_BY_DT2_INC_FIELD
    IMPLICIT NONE
    !*******************************************************************************
    ! Real-valued Gaussian quadrature along an edge for 
    ! the product of a vector basis function of type VBF_type and the incident field 
    ! of time derivative order as specified. 
    ! This is the incident field in the 4th term in (12.80) (Jin 2nd), for scattered field analysis.
    ! 
    ! Note that this function supports ONLY the Webb99 type elements, since some analytical preprocessing
    ! used here relies on this. The routine supports coefficient matching as described in JP Webb,
    ! "Matching a Given Field using hierarchal basis functions", Electromagnetics, 2003 (to appear). 
    !
    ! 1, 2, 3 and 4 point rules are available. 
    !
    ! 04 June 2003: First version. DBD.
    ! 
    ! Note that this function is grouped under basis_function, otherwise circular
    ! referencing results. 
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: element_num,local_edge_num,VBF_type,time_derivative_order,num_qpoints
    REAL(SP), INTENT(IN) :: ell        ! Edge length

    REAL(SP) :: REAL_GAUSS_QUAD_LINE_VBF
    REAL(SP), DIMENSION(3) :: tempcoord
    REAL(SP) :: x,y,z                             ! 3D coordinates. 
    REAL(SP), DIMENSION(3) :: temp_vec1,temp_vec2 ! Temporary vectors. 
    REAL(SP) :: N2eg,N1eg                         ! See eqn. (24) Webb03.  
    INTEGER(I4B) :: ii,inode                      ! Misc. counter.
    REAL(SP), DIMENSION(4) :: lambda              ! Simplex coordinates of point
    REAL(SP), DIMENSION(4,3) :: nabla_lambda      ! Gradients of simplex
    ! coordinates. Constant
    ! within element.
    INTEGER(I4B), DIMENSION(2) :: tempedgenodes   ! local nodes of the edge

    ! Test for correct data:
    IF (SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: Edge length scaling not supported in GAUSS_QUAD_LINE_VBF.'
    END IF
    IF (ELEMENT_TYPE.NE.3) THEN 
       STOP 'IE: Only Webb99 type elements supported in GAUSS_QUAD_LINE_VBF.'
    END IF
    IF (VBF_type.NE.2) THEN 
       STOP 'IE: GAUSS_QUAD_LINE_VBF only supports edge functions of grad type.'
    END IF

    ! Get the (normalized) gradients of the simplex coordinates.
    nabla_lambda =  GRADIENT_LAMBDA(element_num,.true.)

    REAL_GAUSS_QUAD_LINE_VBF = 0.0 ! Initialize
    tempedgenodes = LOCAL_EDGENODES(local_edge_num)

    ! Perform quadrature:
    DO ii = 1,num_qpoints

       ! Associate the correct 3D simplex coordinates with the 1D integration:
       lambda = 0.0
       DO inode = 1,2
          lambda(tempedgenodes(inode)) = quad_line_rules(num_qpoints)%rule(ii,inode)
       END DO

       ! Find values for x,y and z corresponding to sample points (all coordinates used). 
       tempcoord(1:3) = XYZ_COORDINATES(element_num,lambda)
       x = tempcoord(1)
       y = tempcoord(2)
       z = tempcoord(3)

       IF (TD_ANALYSIS) THEN
	  IF(.NOT.SCAT_FIELD) THEN
             STOP 'IE: Unsupported feature in GAUSS_QUAD_LINE_VBF. Should only be called for SCAT_FIELD option'
	  END IF

          temp_vec1 = EDGE_UNIT_VECTOR(tempedgenodes(1),tempedgenodes(2),element_num)
          SELECT CASE (time_derivative_order)
          CASE(0) 
             temp_vec2 = TD_INC_FIELD(x,y,z,'E')
          CASE(1) 
             temp_vec2 = TD_D_BY_DT_INC_FIELD(x,y,z,'E')
          CASE(2)
             temp_vec2 = TD_D2_BY_DT2_INC_FIELD(x,y,z,'E')
          CASE DEFAULT
             STOP 'IE: GAUSS_QUAD_LINE_VBF called with incorrect <derivative_order>'
	  END SELECT

          N2eg = lambda(tempedgenodes(1))
          N1eg = lambda(tempedgenodes(2))

          ! Add weighted function value.  
          REAL_GAUSS_QUAD_LINE_VBF = REAL_GAUSS_QUAD_LINE_VBF +                           &
               quad_line_rules(num_qpoints)%rule(ii,3) * & 
               (N2eg-N1eg)* DOT_PRODUCT(temp_vec1,temp_vec2)
       ELSE 
	  STOP 'IE: Unsupported feature in GAUSS_QUAD_LINE_VBF. Only implemented for TD analysis.'
       END IF

    END DO ! Quadrature loop

    ! Above computes essentially b^{eg}_j, j=1, in Webb03, (24). Scale by edge length when finished:
    REAL_GAUSS_QUAD_LINE_VBF = REAL_GAUSS_QUAD_LINE_VBF * ell

  END FUNCTION REAL_GAUSS_QUAD_LINE_VBF
  !*******************************************************************************




  SUBROUTINE REAL_INTERPOLATE_FUNCTION(time_derivative_order,elem,coefficients,order_request,mixed_request)
    USE geometry
    USE math_tools, ONLY: CROSS_PRODUCT,DET_DIM2
    USE problem_info
    USE TD_source
    IMPLICIT NONE
    !*******************************************************************************
    ! Returns the values of the coefficients of the basis functions that interpolate 
    ! the real-valued time-domain field of derivative order <derivative_order> 
    ! within element <elem>. 
    !
    ! Note that this is only uniquely defined for CT/LN elements. The theory used is 
    ! based on JP Webb, "Matching a given field using hierarchal vector  basis functions",
    ! Electromagnetics, to appear. 
    !
    ! Written DBD 21 May 2003.
    ! Extended DBD 04 June 2003 to include elements up to LT/QN. 
    ! Correction in LT/QN in function BJ terms 31 Aug 2004 DBD.
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: time_derivative_order
    INTEGER(I4B), INTENT(IN) :: elem
    REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE), INTENT(OUT) :: coefficients
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: order_request
    LOGICAL(LGT), OPTIONAL, INTENT(IN) :: mixed_request

    INTEGER(I4B) :: edgenum,facenum,vcount1,vcount2 ! counters
    INTEGER(I4B), DIMENSION(2) :: tempnodes,tempglobalnodes
    INTEGER(I4B), DIMENSION(3) :: tempfacenodes
    REAL(SP), DIMENSION(6) :: edge_lengths
    INTEGER(I4B) :: order_value                   ! Values which will be used to
    LOGICAL(LGT) :: mixed_value                   ! decide the function to be evaluated
    REAL(SP), DIMENSION(3) :: temp1,temp2,avg_field
    REAL(SP) :: x_1,y_1,z_1,x_2,y_2,z_2
    REAL(SP) ::  Meg_11, beg_1                    ! eqn (23) & (24) , [Webb03].
    REAL(SP), DIMENSION(2,2) :: Mfr,A1,A2         ! eqn (29) [Webb03]
    REAL(SP), DIMENSION(2) :: bfr                 ! eqn (31) [Webb03]
    REAL(SP), DIMENSION(4,4,3) :: V_matrix        ! See [eqn(7),S&P 96]
    REAL(SP), DIMENSION(4,3) :: grad_lam
    REAL(SP), DIMENSION(3) :: face_normal
    INTEGER(I4B) :: num_qpoints 

    IF (SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: REAL_INTERPOLATE_FUNCTION does not support edge length scaling.'
    END IF

    ! Establish the order & mixed values to which the functions must be evaluated:
    IF (PRESENT(order_request)) THEN
       order_value = order_request
    ELSE
       order_value = MAX_ORDER(elem)
    END IF
    IF (PRESENT(mixed_request)) THEN
       mixed_value = mixed_request
    ELSE
       mixed_value = MIXED_ORDER(elem)
    END IF

    ! Set all values to zero (This is a precaution to ensure the lower order elements
    ! do not return non-zero values for non-existent higher order terms.)
    coefficients = 0.0_SP

    ! Calculate the 6 edge lengths for repeated later use:
    edge_lengths(1:6) = T_LENGTHS(elem)

    ! Calculate GRADIENT_LAMBDA for repeated later use:
    ! grad_lam = GRADIENT_LAMBDA(elem,.true.) 

    ! CT/LN and higher order terms: 
    ! Evaluate the E1 functions.
    DO edgenum = 1,6
       tempnodes = LOCAL_EDGENODES(edgenum)
       tempglobalnodes   = GLOBAL_EDGENODES(elem,edgenum)
       ! Find field at start node
       x_1 = vertices(tempglobalnodes(1))%coord(1) 
       y_1 = vertices(tempglobalnodes(1))%coord(2) 
       z_1 = vertices(tempglobalnodes(1))%coord(3)
       ! Find field at end node
       x_2 = vertices(tempglobalnodes(2))%coord(1) 
       y_2 = vertices(tempglobalnodes(2))%coord(2) 
       z_2 = vertices(tempglobalnodes(2))%coord(3)
       SELECT CASE (time_derivative_order)
       CASE(0) 
          temp1 = TD_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = TD_INC_FIELD(x_2,y_2,z_2,'E')
       CASE(1) 
          temp1 = TD_D_BY_DT_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = TD_D_BY_DT_INC_FIELD(x_2,y_2,z_2,'E')
       CASE(2)
          temp1 = TD_D2_BY_DT2_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = TD_D2_BY_DT2_INC_FIELD(x_2,y_2,z_2,'E')
       CASE DEFAULT
	  STOP 'IE: REAL_INTERPOLATE_FUNCTION called with incorrect <derivative_order>'
       END SELECT
       ! Find average field
       avg_field = (temp1+temp2)/2
       ! Match coefficient - eqn(12) in Webb. 
       coefficients(edgenum) = DOT_PRODUCT(avg_field,EDGE_UNIT_VECTOR(tempnodes(1),tempnodes(2),elem)) & 
            * edge_lengths(edgenum)
    END DO

    ! LT/LN terms.
    IF (order_value.EQ.2.OR.&
         .NOT.mixed_value.AND.order_value.EQ.1) THEN
       print *,'Warning: this LT/LN implementation of REAL_INTERPOLATE_FUNCTION has not yet been validated.'

       ! Match edge gradient functions. Note: following is special case, valid only up to LT/QN order.
       ! This is eqn. 22-24, [Webb03]
       DO edgenum = 1,6
          beg_1 = REAL_GAUSS_QUAD_LINE_VBF(elem,edgenum,(2),time_derivative_order,(2),edge_lengths(edgenum))
          Meg_11 = 1.0_SP/3.0_SP ! Computed analytically beforehand.
          coefficients(6+edgenum) = beg_1 / Meg_11
       END DO
    END IF

    ! LT/QN terms.
    IF (order_value.EQ.2) THEN
       print *,'Warning: this LT/QN implementation of REAL_INTERPOLATE_FUNCTION has not yet been validated.'
       ! Calculate GRADIENT_LAMBDA and [V] for repeated later use:
       grad_lam = GRADIENT_LAMBDA(elem,.true.) 
       DO vcount1 = 1,4
          DO vcount2 = 1,4
             V_matrix(vcount1,vcount2,1:3) = &
                  CROSS_PRODUCT(grad_lam(vcount1,1:3),grad_lam(vcount2,1:3))
          END DO
       END DO
       ! Match face rot functions. Eqns. 29-31 [Webb03] 
       DO facenum = 1,4
	  face_normal = ELEMENT_NORMAL_VECTOR(elem,facenum)
          tempfacenodes = LOCAL_FACENODES(facenum)
	  num_qpoints = 3 ! degree of precision 2, sufficient for two first order functions. 
          Mfr(1,1) =  MJI(3,3,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(1,2) =  MJI(3,4,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(2,1) =  MJI(4,3,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(2,2) =  MJI(4,4,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          ! Find bfr:
          bfr(1) = BJ(3,time_derivative_order,elem,V_matrix,face_normal,tempfacenodes,num_qpoints)
          bfr(2) = BJ(4,time_derivative_order,elem,V_matrix,face_normal,tempfacenodes,num_qpoints)

	  ! Solve using Cramer's rule for speed:
	  A1 = Mfr
	  A2 = Mfr
	  A1(1,1) = bfr(1)
	  A1(2,1) = bfr(2)
	  A2(1,2) = bfr(1)
	  A2(2,2) = bfr(2)
          coefficients(12+facenum)  = DET_DIM2(A1)/DET_DIM2(Mfr)
          coefficients(16+facenum)  = DET_DIM2(A2)/DET_DIM2(Mfr)
       END DO
    END IF

    ! QT/QN and higher order terms.
    IF (order_value.GE.3.OR.&
         .NOT.mixed_value.AND.order_value.GE.2) THEN
       STOP 'Only up to LT/QN elements implemented in INTERPOLATE_FUNCTION'
    END IF


  CONTAINS 



    FUNCTION BJ(VBF_j,time_derivative_order,elem,V_matrix,normal,local_face_nodes,num_qpoints)
      USE nrtype
      USE material_properties
      IMPLICIT NONE
      !*******************************************************************************
      ! This function computes the equivalent of eqn(31) using quadrature. 
      ! The result is returned WITHOUT the area, since this will also be 
      ! omitted when evaluating the equivalent of (30). 
      !*******************************************************************************
      REAL(SP) :: BJ
      INTEGER(I4B), INTENT(IN) ::  VBF_j,time_derivative_order,elem,num_qpoints
      REAL(SP), DIMENSION(4,4,3),INTENT(IN) :: V_matrix        ! See [eqn(7),S&P 96]
      INTEGER(I4B), DIMENSION(3), INTENT(IN)  :: local_face_nodes
      REAL(SP), DIMENSION(3), INTENT(IN) :: normal
      REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
      REAL(SP), DIMENSION(3) :: curl_Njfr, curl_E
      REAL(SP) :: N_dot_curl_Njfr, N_dot_curl_E
      REAL(SP), DIMENSION(3) :: tempcoord
      INTEGER(I4B) :: quad_point,inode,i1,i2,i3
      REAL(SP) :: Y_char, x, y, z
      BJ = 0.0_SP
      DO quad_point = 1,num_qpoints
         ! Associate the correct 3D simplex coordinates with the 2D integration:
         lambda = 0.0
         DO inode = 1,3
            lambda(local_face_nodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
         END DO

         i1 = local_face_nodes(1)
         i2 = local_face_nodes(2)
         i3 = local_face_nodes(3)
         curl_Njfr = CURL_VBF(VBF_j,lambda,V_matrix,i1,i2,i3)
         N_dot_curl_Njfr = DOT_PRODUCT(normal,curl_Njfr )

         ! Find \curl E; this is computed from the time derivative of H. (\curl e = - dB/dt). Hence the TD source must 
         ! be called with _derivative_order+1
         Y_char = REAL(eps_r(HOMOG_MEDIUM)/mu_r(HOMOG_MEDIUM))/Z_zero ! Real-valued characteristic admittance of 
         ! homogeneous background medium. 

         ! Find values for x,y and z corresponding to sample points (all coordinates used). 
         tempcoord(1:3) = XYZ_COORDINATES(elem,lambda)
         x = tempcoord(1)
         y = tempcoord(2)
         z = tempcoord(3)

         SELECT CASE (time_derivative_order)
         CASE(0) 
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* TD_D_BY_DT_INC_FIELD(x,y,z,'H',Y_char)
         CASE(1)
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* TD_D2_BY_DT2_INC_FIELD(x,y,z,'H',Y_char)
         CASE(2)
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* TD_D3_BY_DT3_INC_FIELD(x,y,z,'H',Y_char)
         CASE DEFAULT
            STOP 'IE: BJ called with incorrect <derivative_order>'
         END SELECT
         N_dot_curl_E = DOT_PRODUCT(normal,curl_E )
         ! Add the contribution of this quadrature point:
         BJ = BJ + quad_tri_rules(num_qpoints)%rule(quad_point,4) * N_dot_curl_Njfr*N_dot_curl_E
      END DO

      !Scaling  by area not needed. 

    END FUNCTION BJ

  END SUBROUTINE REAL_INTERPOLATE_FUNCTION
  !*******************************************************************************


  FUNCTION MJI(VBF_j,VBF_i,V_matrix,normal,local_face_nodes,num_qpoints)
    USE nrtype
    IMPLICIT NONE
    !*******************************************************************************
    ! This function computes the equivalent of eqn(30) using quadrature. 
    ! The result is returned WITHOUT the are, since this will also be 
    ! omitted when evaluating the equivalent of (31). Since this is based only
    ! the basis functions, the real-valued function works for both
    ! real and complex-valued fields. 
    !*******************************************************************************
    REAL(SP) :: MJI
    INTEGER(I4B), INTENT(IN) ::  VBF_j,VBF_i,num_qpoints
    REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
    REAL(SP), DIMENSION(4,4,3),INTENT(IN) :: V_matrix        ! See [eqn(7),S&P 96]
    REAL(SP), DIMENSION(3), INTENT(IN) :: normal
    INTEGER(I4B), DIMENSION(3), INTENT(IN)  :: local_face_nodes
    REAL(SP), DIMENSION(3) :: curl_Njfr, curl_Nifr 
    REAL(SP)  :: N_dot_curl_Njfr, N_dot_curl_Nifr 
    INTEGER(I4B) :: quad_point,inode,i1,i2,i3

    MJI = 0.0_SP
    DO quad_point = 1,num_qpoints
       ! Associate the correct 3D simplex coordinates with the 2D integration:
       lambda = 0.0
       DO inode = 1,3
          lambda(local_face_nodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
       END DO

       i1 = local_face_nodes(1)
       i2 = local_face_nodes(2)
       i3 = local_face_nodes(3)
       curl_Njfr = CURL_VBF(VBF_j,lambda,V_matrix,i1,i2,i3)
       N_dot_curl_Njfr = DOT_PRODUCT(normal,curl_Njfr )
       curl_Nifr = CURL_VBF(VBF_i,lambda,V_matrix,i1,i2,i3)
       N_dot_curl_Nifr = DOT_PRODUCT(normal,curl_Nifr )
       ! Add the contribution of this quadrature point:
       MJI = MJI + quad_tri_rules(num_qpoints)%rule(quad_point,4) * N_dot_curl_Njfr * N_dot_curl_Nifr
    END DO

    !Scaling  by area not needed. 

  END FUNCTION MJI


  SUBROUTINE COMPLEX_INTERPOLATE_FUNCTION(time_derivative_order,elem,coefficients,order_request,mixed_request)
    USE frequency_data
    USE geometry
    USE math_tools, ONLY: CROSS_PRODUCT,DET_DIM2
    USE problem_info
    USE FD_scat_source
    IMPLICIT NONE
    !*******************************************************************************
    ! Returns the values of the coefficients of the basis functions that interpolate 
    ! the complex-valued frequency-domain field of derivative order <derivative_order> 
    ! within element <elem>. 
    !
    ! Note that this is only uniquely defined for CT/LN elements. The theory used is 
    ! based on JP Webb, "Matching a given field using hierarchal vector  basis functions",
    ! Electromagnetics, to appear. 
    !
    ! Based on REAL_INTERPOLATE_FUNCTION, and written so that arguments are compatible.
    ! Written DBD 07 Aug 2004. 
    ! Extended up to LT/QN 31 Aug 04 DBD.
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: time_derivative_order
    INTEGER(I4B), INTENT(IN) :: elem
    COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), INTENT(OUT) :: coefficients
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: order_request
    LOGICAL(LGT), OPTIONAL, INTENT(IN) :: mixed_request

    INTEGER(I4B) :: edgenum,facenum,vcount1,vcount2 ! counters
    INTEGER(I4B), DIMENSION(2) :: tempnodes,tempglobalnodes
    INTEGER(I4B), DIMENSION(3) :: tempfacenodes
    REAL(SP), DIMENSION(6) :: edge_lengths
    INTEGER(I4B) :: order_value                   ! Values which will be used to
    LOGICAL(LGT) :: mixed_value                   ! decide the function to be evaluated
    COMPLEX(SPC), DIMENSION(3) :: temp1,temp2,avg_field
    REAL(SP) :: x_1,y_1,z_1,x_2,y_2,z_2
    REAL(SP) ::  Meg_11                       ! eqn (23) & (24) , [Webb03].
    COMPLEX(SPC) ::  beg_1                        ! eqn (23) & (24) , [Webb03].
    REAL(SP), DIMENSION(2,2) :: Mfr               ! eqn (29) [Webb03]
    COMPLEX(SPC), DIMENSION(2,2) :: A1,A2         ! eqn (29) [Webb03]
    COMPLEX(SPC), DIMENSION(2) :: bfr             ! eqn (31) [Webb03]
    REAL(SP), DIMENSION(4,4,3) :: V_matrix        ! See [eqn(7),S&P 96]
    REAL(SP), DIMENSION(4,3) :: grad_lam
    REAL(SP), DIMENSION(3) :: face_normal
    INTEGER(I4B) :: num_qpoints 

    IF (SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: COMPLEX_INTERPOLATE_FUNCTION does not support edge length scaling.'
    END IF

    ! Establish the order & mixed values to which the functions must be evaluated:
    IF (PRESENT(order_request)) THEN
       order_value = order_request
    ELSE
       order_value = MAX_ORDER(elem)
    END IF
    IF (PRESENT(mixed_request)) THEN
       mixed_value = mixed_request
    ELSE
       mixed_value = MIXED_ORDER(elem)
    END IF

    ! Set all values to zero (This is a precaution to ensure the lower order elements
    ! do not return non-zero values for non-existent higher order terms.)
    coefficients = 0.0_SPC

    ! Calculate the 6 edge lengths for repeated later use:
    edge_lengths(1:6) = T_LENGTHS(elem)

    ! Calculate GRADIENT_LAMBDA for repeated later use:
    ! grad_lam = GRADIENT_LAMBDA(elem,.true.) 

    ! CT/LN and higher order terms: 
    ! Evaluate the E1 functions.
    DO edgenum = 1,6
       tempnodes = LOCAL_EDGENODES(edgenum)
       tempglobalnodes   = GLOBAL_EDGENODES(elem,edgenum)
       ! Find field at start node
       x_1 = vertices(tempglobalnodes(1))%coord(1) 
       y_1 = vertices(tempglobalnodes(1))%coord(2) 
       z_1 = vertices(tempglobalnodes(1))%coord(3)
       ! Find field at end node
       x_2 = vertices(tempglobalnodes(2))%coord(1) 
       y_2 = vertices(tempglobalnodes(2))%coord(2) 
       z_2 = vertices(tempglobalnodes(2))%coord(3)
       SELECT CASE (time_derivative_order)
       CASE(0) 
          temp1 = FD_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = FD_INC_FIELD(x_2,y_2,z_2,'E')
       CASE(1) 
          temp1 = j*2.0_SP*PI*frequency*FD_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = j*2.0_SP*PI*frequency*FD_INC_FIELD(x_2,y_2,z_2,'E')
       CASE(2)
          temp1 = -(2.0_SP*PI*frequency)**2*FD_INC_FIELD(x_1,y_1,z_1,'E')
          temp2 = -(2.0_SP*PI*frequency)**2*FD_INC_FIELD(x_2,y_2,z_2,'E')
       CASE DEFAULT
	  STOP 'IE: COMPLEX_INTERPOLATE_FUNCTION called with incorrect <derivative_order>'
       END SELECT
       ! Find average field
       avg_field = (temp1+temp2)/2
       ! Match coefficient - eqn(12) in Webb. Note that DOT_PRODUCT conjugates the first argument, 
       ! which would be incorrect in this context, hence the additional conjugation. 
       coefficients(edgenum) = DOT_PRODUCT(CONJG(avg_field),EDGE_UNIT_VECTOR(tempnodes(1),tempnodes(2),elem)) & 
            * edge_lengths(edgenum)
    END DO

    ! LT/LN terms.
    IF (order_value.EQ.2.OR.&
         .NOT.mixed_value.AND.order_value.EQ.1) THEN
       ! Match edge gradient functions. Note: following is special case, valid only up to LT/QN order.
       ! This is eqn. 22-24, [Webb03]
       DO edgenum = 1,6
          beg_1 = COMPLEX_GAUSS_QUAD_LINE_VBF(elem,edgenum,(2),time_derivative_order,(2),edge_lengths(edgenum))
          Meg_11 = 1.0_SP/3.0_SP ! Computed analytically beforehand.
          coefficients(6+edgenum) = beg_1 / Meg_11
       END DO
    END IF

    ! LT/QN terms.
    IF (order_value.EQ.2) THEN
       ! Calculate GRADIENT_LAMBDA and [V] for repeated later use:
       grad_lam = GRADIENT_LAMBDA(elem,.true.) 
       DO vcount1 = 1,4
          DO vcount2 = 1,4
             V_matrix(vcount1,vcount2,1:3) = &
                  CROSS_PRODUCT(grad_lam(vcount1,1:3),grad_lam(vcount2,1:3))
          END DO
       END DO
       ! Match face rot functions. Eqns. 29-31 [Webb03] 
       DO facenum = 1,4
	  face_normal = ELEMENT_NORMAL_VECTOR(elem,facenum)
          tempfacenodes = LOCAL_FACENODES(facenum)
	  num_qpoints = 3 ! degree of precision 2, sufficient for two first order functions. 
          Mfr(1,1) =  MJI(3,3,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(1,2) =  MJI(3,4,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(2,1) =  MJI(4,3,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          Mfr(2,2) =  MJI(4,4,V_matrix,face_normal,tempfacenodes,num_qpoints) 
          ! Find bfr:
          bfr(1) = BJ(3,time_derivative_order,elem,V_matrix,face_normal,tempfacenodes,num_qpoints)
          bfr(2) = BJ(4,time_derivative_order,elem,V_matrix,face_normal,tempfacenodes,num_qpoints)

	  ! Solve using Cramer's rule for speed:
	  A1 = Mfr
	  A2 = Mfr
	  A1(1,1) = bfr(1)
	  A1(2,1) = bfr(2)
	  A2(1,2) = bfr(1)
	  A2(2,2) = bfr(2)
          coefficients(12+facenum)  = DET_DIM2(A1)/DET_DIM2(Mfr)
          coefficients(16+facenum)  = DET_DIM2(A2)/DET_DIM2(Mfr)
       END DO
    END IF

    ! QT/QN and higher order terms.
    IF (order_value.GE.3.OR.&
         .NOT.mixed_value.AND.order_value.GE.2) THEN
       STOP 'Only up to LT/QN elements implemented in INTERPOLATE_FUNCTION'
    END IF

  CONTAINS

    FUNCTION BJ(VBF_j,time_derivative_order,elem,V_matrix,normal,local_face_nodes,num_qpoints)
      USE nrtype
      USE material_properties
      IMPLICIT NONE
      !*******************************************************************************
      ! This function computes the equivalent of eqn(31) using quadrature. 
      ! The result is returned WITHOUT the area, since this will also be 
      ! omitted when evaluating the equivalent of (30). 
      !*******************************************************************************
      COMPLEX(SPC) :: BJ
      INTEGER(I4B), INTENT(IN) ::  VBF_j,time_derivative_order,elem,num_qpoints
      REAL(SP), DIMENSION(4,4,3),INTENT(IN) :: V_matrix        ! See [eqn(7),S&P 96]
      INTEGER(I4B), DIMENSION(3), INTENT(IN)  :: local_face_nodes
      REAL(SP), DIMENSION(3), INTENT(IN) :: normal
      REAL(SP), DIMENSION(4) :: lambda             ! Simplex coordinates of point
      REAL(SP), DIMENSION(3) :: curl_Njfr
      COMPLEX(SPC), DIMENSION(3) :: curl_E
      REAL(SP) :: N_dot_curl_Njfr
      COMPLEX(SPC) :: N_dot_curl_E
      REAL(SP), DIMENSION(3) :: tempcoord
      INTEGER(I4B) :: quad_point,inode,i1,i2,i3
      REAL(SP) :: Y_char, x, y, z
      BJ = 0.0_SPC
      DO quad_point = 1,num_qpoints
         ! Associate the correct 3D simplex coordinates with the 2D integration:
         lambda = 0.0
         DO inode = 1,3
            lambda(local_face_nodes(inode)) = quad_tri_rules(num_qpoints)%rule(quad_point,inode)
         END DO

         i1 = local_face_nodes(1)
         i2 = local_face_nodes(2)
         i3 = local_face_nodes(3)
         curl_Njfr = CURL_VBF(VBF_j,lambda,V_matrix,i1,i2,i3)
         N_dot_curl_Njfr = DOT_PRODUCT(normal,curl_Njfr )

         ! Find \curl E; this is computed from the time derivative of H. (\curl e = - dB/dt). Hence the TD source must 
         ! be called with _derivative_order+1 
         Y_char = REAL(eps_r(HOMOG_MEDIUM)/mu_r(HOMOG_MEDIUM))/Z_zero ! Real-valued characteristic admittance of 
         ! homogeneous background medium. 

         ! Find values for x,y and z corresponding to sample points (all coordinates used). 
         tempcoord(1:3) = XYZ_COORDINATES(elem,lambda)
         x = tempcoord(1)
         y = tempcoord(2)
         z = tempcoord(3)

         SELECT CASE (time_derivative_order)
         CASE(0) 
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* j*2.0_SP*pi*frequency*FD_INC_FIELD(x,y,z,'H',Y_char)
         CASE(1)
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* (-2.0_SP*pi*frequency)**2*FD_INC_FIELD(x,y,z,'H',Y_char)
         CASE(2)
            curl_E = -mu_0*mu_r(HOMOG_MEDIUM)* (j*2.0_SP*pi*frequency)**3*FD_INC_FIELD(x,y,z,'H',Y_char)
         CASE DEFAULT
            STOP 'IE: BJ called with incorrect <derivative_order>'
         END SELECT
         N_dot_curl_E = DOT_PRODUCT(normal,curl_E )
         ! Add the contribution of this quadrature point:
         BJ = BJ + quad_tri_rules(num_qpoints)%rule(quad_point,4) * N_dot_curl_Njfr*N_dot_curl_E
      END DO

      !Scaling  by area not needed. 

    END FUNCTION BJ

  END SUBROUTINE COMPLEX_INTERPOLATE_FUNCTION



  FUNCTION COMPLEX_GAUSS_QUAD_LINE_VBF(element_num,local_edge_num,&
       VBF_type,time_derivative_order,num_qpoints,ell)

    USE boundary_conditions
    USE frequency_data
    USE geometry
    USE problem_info
    USE quad_tables
    USE FD_scat_source, ONLY: FD_INC_FIELD
    IMPLICIT NONE
    !*******************************************************************************
    ! Real-valued Gaussian quadrature along an edge for 
    ! the product of a vector basis function of type VBF_type and the incident field 
    ! of time derivative order as specified. 
    ! This is the incident field in the 4th term in (12.80) (Jin 2nd), for scattered field analysis.
    ! 
    ! Note that this function supports ONLY the Webb99 type elements, since some analytical preprocessing
    ! used here relies on this. The routine supports coefficient matching as described in JP Webb,
    ! "Matching a Given Field using hierarchal basis functions", Electromagnetics, 2003 (to appear). 
    !
    ! 1, 2, 3 and 4 point rules are available. 
    !
    ! 04 June 2003: First version. DBD.
    ! 
    ! Note that this function is grouped under basis_function, otherwise circular
    ! referencing results. 
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: element_num,local_edge_num,VBF_type,time_derivative_order,num_qpoints
    REAL(SP), INTENT(IN) :: ell        ! Edge length

    COMPLEX(SPC) :: COMPLEX_GAUSS_QUAD_LINE_VBF
    REAL(SP), DIMENSION(3) :: tempcoord
    REAL(SP) :: x,y,z                             ! 3D coordinates. 
    REAL(SP), DIMENSION(3) :: temp_vec1           ! Temporary vector. 
    COMPLEX(SPC), DIMENSION(3) :: temp_vec2_c     ! Temporary vector. 
    REAL(SP) :: N2eg,N1eg                         ! See eqn. (24) Webb03.  
    INTEGER(I4B) :: ii,inode                      ! Misc. counter.
    REAL(SP), DIMENSION(4) :: lambda              ! Simplex coordinates of point
    REAL(SP), DIMENSION(4,3) :: nabla_lambda      ! Gradients of simplex
    ! coordinates. Constant
    ! within element.
    INTEGER(I4B), DIMENSION(2) :: tempedgenodes   ! local nodes of the edge

    ! Test for correct data:
    IF (SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: Edge length scaling not supported in GAUSS_QUAD_LINE_VBF.'
    END IF
    IF (ELEMENT_TYPE.NE.3) THEN 
       STOP 'IE: Only Webb99 type elements supported in GAUSS_QUAD_LINE_VBF.'
    END IF
    IF (VBF_type.NE.2) THEN 
       STOP 'IE: GAUSS_QUAD_LINE_VBF only supports edge functions of grad type.'
    END IF

    ! Get the (normalized) gradients of the simplex coordinates.
    nabla_lambda =  GRADIENT_LAMBDA(element_num,.true.)

    COMPLEX_GAUSS_QUAD_LINE_VBF = 0.0 ! Initialize
    tempedgenodes = LOCAL_EDGENODES(local_edge_num)

    ! for debugging

    !if (element_num.EQ.27.AND.local_edge_num.EQ.6) then
    !write(fileout,*) 'Working in element 27, edge 6'
    !end if


    ! Perform quadrature:
    DO ii = 1,num_qpoints

       ! Associate the correct 3D simplex coordinates with the 1D integration:
       lambda = 0.0
       DO inode = 1,2
          lambda(tempedgenodes(inode)) = quad_line_rules(num_qpoints)%rule(ii,inode)
       END DO

       ! Find values for x,y and z corresponding to sample points (all coordinates used). 
       tempcoord(1:3) = XYZ_COORDINATES(element_num,lambda)
       x = tempcoord(1)
       y = tempcoord(2)
       z = tempcoord(3)

       temp_vec1 = EDGE_UNIT_VECTOR(tempedgenodes(1),tempedgenodes(2),element_num)

       SELECT CASE (time_derivative_order)
       CASE(0) 
          temp_vec2_c = FD_INC_FIELD(x,y,z,'E')
       CASE(1) 
          temp_vec2_c = j*2.0_SP*pi*frequency*FD_INC_FIELD(x,y,z,'E')
       CASE(2)
          temp_vec2_c = -(2.0_SP*pi*frequency)**2*FD_INC_FIELD(x,y,z,'E')
       CASE DEFAULT
          STOP 'IE: GAUSS_QUAD_LINE_VBF called with incorrect <derivative_order>'
       END SELECT

       N2eg = lambda(tempedgenodes(1)) 
       N1eg = lambda(tempedgenodes(2))

       ! Add weighted function value.  
       COMPLEX_GAUSS_QUAD_LINE_VBF = COMPLEX_GAUSS_QUAD_LINE_VBF +                           &
            quad_line_rules(num_qpoints)%rule(ii,3) * & 
            (N2eg-N1eg)* DOT_PRODUCT(CMPLX(temp_vec1),temp_vec2_c)

    END DO ! Quadrature loop

    ! Above computes essentially b^{eg}_j, j=1, in Webb03, (24). Scale by edge length when finished:
    COMPLEX_GAUSS_QUAD_LINE_VBF = COMPLEX_GAUSS_QUAD_LINE_VBF * ell

  END FUNCTION COMPLEX_GAUSS_QUAD_LINE_VBF
  !*******************************************************************************


  SUBROUTINE POINT_EVALUATE_FACE_FUNCTIONS(elem,face_num,xyz,num_vbfs,term1,term2)
    USE geometry
    USE nrtype
    IMPLICIT NONE   
    !*******************************************************************************
    ! Calculates the values of the edge/face functions at the point <xyz>/<lam>, in the 
    ! local face <face_num> of element <elem>. It is assumed that the face and evaluation
    ! point lies in the z=0 plane. t1 and t2 refer to (del dot ^z cross basis functions)
    ! and (^z cross basis functions). Values for the lowest order set of <num_vbfs> 
    ! vector basis functions are returned.
    !
    ! 2000-11-22: Created. MMB.
    ! 2001-09-27: Added save capability and xyz/lam co-ords option. MMB.
    ! 2002-05-11: Changed to call by a number of vector basis functions, rather
    !             than specific arguments for CTLN/LTQN etc. . MMB.
    ! 2002-05-17: Added Webb99 elements to 2nd order complete. MMB.
    ! 2002-05-20: MAJOR BUG found!! <z_cross_grad_lam(count1,1:3)> was incorrectly
    !             <z_cross_grad_lam(1,1:3)> in the data setup loop. MMB.
    !*******************************************************************************
    SAVE
    INTEGER(I4B), INTENT(IN) :: elem
    INTEGER(I4B), INTENT(IN) :: face_num
    REAL(SP), DIMENSION(3), INTENT(IN) :: xyz
    INTEGER(I4B), INTENT(IN) :: num_vbfs
    REAL(SP), DIMENSION(num_vbfs), INTENT(OUT) :: term1
    REAL(SP), DIMENSION(num_vbfs,3), INTENT(OUT) :: term2

    INTEGER(I4B) :: old_elem = 0
    INTEGER(I4B) :: count1
    REAL(SP), DIMENSION(4) :: simp,cart
    REAL(SP), DIMENSION(4,4) :: simat
    INTEGER(I4B), DIMENSION(2) :: tempnodes
    INTEGER(I4B), DIMENSION(3) :: tempedges
    INTEGER(I4B), DIMENSION(3) :: facenodes
    REAL(SP), DIMENSION(3) :: edge_lengths,vec1,vec2,tempvec
    REAL(SP), DIMENSION(4,3) :: grad_lam,z_cross_grad_lam
    REAL(SP), DIMENSION(3) :: t1_f1_coeff,t1_f2_coeff

    IF_OLD_NEW: IF (old_elem.NE.elem) THEN ! new element-specific data must be calculated

       ! Save for future use:
       old_elem = elem

       ! Calculate the simplex and vertices matrices:
       CALL SIMPLEX_COEFFICIENTS(elem,simat(1:4,4),simat(1:4,1), &
            simat(1:4,2),simat(1:4,3))

       ! Calculate the three edge lengths for repeated later use:
       tempedges = LOCAL_FACEEDGES(face_num)
       DO count1 = 1,3 ! 3 edges
          tempnodes = LOCAL_EDGENODES(tempedges(count1))
          tempnodes = elements(elem)%nodes(tempnodes) ! global numbers for length fuction
          edge_lengths(count1) = T_LENGTH(tempnodes(1),tempnodes(2))
       END DO

       ! Calculate for simplex assignment and F1,F2 related calculations:
       facenodes = LOCAL_FACENODES(face_num)

       ! Calculate GRADIENT_LAMBDA of the face simplex co-ordinates for repeated later use:
       grad_lam = GRADIENT_LAMBDA(elem,.true.) 

       ! Calculate (^z times <grad_lam>):
       DO count1 = 1,4
	  z_cross_grad_lam(count1,1:3) = (/ -grad_lam(count1,2), grad_lam(count1,1), 0.0 /)
       END DO

       ! Calculate t1 F1,F2 coefficients for repeated, later use: (see MMB's own notes (21 Nov 2000))
       t1_f1_coeff(1) =                                        &
            3.0 * simat(facenodes(1),1) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(2),2) -    &
            simat(facenodes(2),1)*simat(facenodes(3),2) )    
       t1_f1_coeff(2) =                                        &
            3.0 * simat(facenodes(1),2) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(2),2) -    &
            simat(facenodes(2),1)*simat(facenodes(3),2) )    
       t1_f1_coeff(3) =                                        &
            2.0 * simat(facenodes(1),4) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(2),2) -    &
            simat(facenodes(2),1)*simat(facenodes(3),2) )    &
            + simat(facenodes(1),2) *                            &
            ( simat(facenodes(2),4)*simat(facenodes(3),1) -    &
            simat(facenodes(3),4)*simat(facenodes(2),1) )    &
            + simat(facenodes(1),1) *                            &
            ( simat(facenodes(3),4)*simat(facenodes(2),2) -    &
            simat(facenodes(2),4)*simat(facenodes(3),2) )
       t1_f2_coeff(1) =                                        &
            3.0 * simat(facenodes(2),1) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(1),2) -    &
            simat(facenodes(1),1)*simat(facenodes(3),2) )    
       t1_f2_coeff(2) =                                        &
            3.0 * simat(facenodes(2),2) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(1),2) -    &
            simat(facenodes(1),1)*simat(facenodes(3),2) )    
       t1_f2_coeff(3) =                                        &
            2.0 * simat(facenodes(2),4) *                      &
            ( simat(facenodes(3),1)*simat(facenodes(1),2) -    &
            simat(facenodes(1),1)*simat(facenodes(3),2) )    &
            + simat(facenodes(2),2) *                            &
            ( simat(facenodes(1),4)*simat(facenodes(3),1) -    &
            simat(facenodes(3),4)*simat(facenodes(1),1) )    &
            + simat(facenodes(2),1) *                            &
            ( simat(facenodes(3),4)*simat(facenodes(1),2) -    &
            simat(facenodes(1),4)*simat(facenodes(3),2) )

    END IF IF_OLD_NEW

    ! Calculate the simplex co-ords of the evaluation point:
    cart(1:3) = xyz
    cart(4) = 1.0
    simp = MATMUL(simat,cart) ! simplex co-ordinates of the point xyz

    ! <num_vbfs> must at least be =1:
    IF (num_vbfs.LT.1) STOP 'IE: <num_vbfs> less than 1 in POINT_EVALUATE_FACE_FUNCTIONS.'

    ! Evaluate E1:
    DO count1 = 1,3 ! 3 edges         
       IF (num_vbfs.LT.count1) RETURN
       tempnodes = LOCAL_EDGENODES(tempedges(count1))
       vec1 = simp(tempnodes(1)) * grad_lam(tempnodes(2),1:3)
       vec2 = simp(tempnodes(2)) * grad_lam(tempnodes(1),1:3)
       tempvec =  vec1 - vec2 
       term2(count1,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(count1) = 2.0*( simat(tempnodes(2),1)*simat(tempnodes(1),2) -  &
            simat(tempnodes(2),2)*simat(tempnodes(1),1) )
       IF (SCALE_BY_EDGE_LENGTH) THEN
          term2(count1,1:3) = edge_lengths(count1) * term2(count1,1:3)
          term1(count1)     = edge_lengths(count1) * term1(count1)
       END IF
    END DO

    ! Evaluate E2:
    DO count1 = 4,6 ! 3 edges
       IF (num_vbfs.LT.count1) RETURN
       tempnodes = LOCAL_EDGENODES(tempedges(count1-3))
       vec1 = simp(tempnodes(1)) * grad_lam(tempnodes(2),1:3)
       vec2 = simp(tempnodes(2)) * grad_lam(tempnodes(1),1:3)
       tempvec = vec1 + vec2
       term2(count1,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(count1) = 0.0
       IF (SCALE_BY_EDGE_LENGTH) THEN
          term2(count1,1:3) = edge_lengths(count1-3) * term2(count1,1:3)
          term1(count1)     = edge_lengths(count1-3) * term1(count1)
       END IF
    END DO

    ! Evaluate F1: (see SUBROUTINE VBF for definitions of F1 and F2)
    IF (num_vbfs.LT.7) RETURN
    SELECT CASE(ELEMENT_TYPE)
    CASE(1) ! Savage
       tempvec = simp(facenodes(1)) *                             &
            ( simp(facenodes(2))*grad_lam(facenodes(3),1:3) -  &
            simp(facenodes(3))*grad_lam(facenodes(2),1:3) )
       term2(7,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(7)     = & ! see MMB's own notes (21 Nov 2000)
            t1_f1_coeff(1)*xyz(1) + t1_f1_coeff(2)*xyz(2) + t1_f1_coeff(3)
    CASE(3) ! Webb99
       tempvec =       simp(facenodes(2))*simp(facenodes(3))*grad_lam(facenodes(1),1:3)  &
            + simp(facenodes(1))*simp(facenodes(3))*grad_lam(facenodes(2),1:3)  &
            - 2.0*simp(facenodes(1))*simp(facenodes(2))*grad_lam(facenodes(3),1:3)
       term2(7,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(7)     = &
            DOT_PRODUCT(z_cross_grad_lam(facenodes(1),1:3),                                                &
            simp(facenodes(2))*grad_lam(facenodes(3),1:3) + simp(facenodes(3))*grad_lam(facenodes(2),1:3)) &
            +     DOT_PRODUCT(z_cross_grad_lam(facenodes(2),1:3),                                                &
            simp(facenodes(1))*grad_lam(facenodes(3),1:3) + simp(facenodes(3))*grad_lam(facenodes(1),1:3)) &
            - 2.0*DOT_PRODUCT(z_cross_grad_lam(facenodes(3),1:3),                                                &
            simp(facenodes(1))*grad_lam(facenodes(2),1:3) + simp(facenodes(2))*grad_lam(facenodes(1),1:3))
    CASE DEFAULT
       STOP 'IE: element type not supported in POINT_EVALUATE_FACE_FUNCTIONS.'
    END SELECT

    ! Evaluate F2:
    IF (num_vbfs.LT.8) RETURN
    SELECT CASE(ELEMENT_TYPE)
    CASE(1) ! Savage
       tempvec = simp(facenodes(2)) *                                 &
            ( simp(facenodes(1))*grad_lam(facenodes(3),1:3) -  &
            simp(facenodes(3))*grad_lam(facenodes(1),1:3) )
       term2(8,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(8)     = & ! see MMB's own notes (21 Nov 2000)
            t1_f2_coeff(1)*xyz(1) + t1_f2_coeff(2)*xyz(2) + t1_f2_coeff(3)
    CASE(3) ! Webb99
       tempvec =       simp(facenodes(3))*simp(facenodes(1))*grad_lam(facenodes(2),1:3)  &
            + simp(facenodes(2))*simp(facenodes(1))*grad_lam(facenodes(3),1:3)  &
            - 2.0*simp(facenodes(2))*simp(facenodes(3))*grad_lam(facenodes(1),1:3)
       term2(8,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(8)     = &
            DOT_PRODUCT(z_cross_grad_lam(facenodes(2),1:3),                                                &
            simp(facenodes(3))*grad_lam(facenodes(1),1:3) + simp(facenodes(1))*grad_lam(facenodes(3),1:3)) &
            +     DOT_PRODUCT(z_cross_grad_lam(facenodes(3),1:3),                                                &
            simp(facenodes(2))*grad_lam(facenodes(1),1:3) + simp(facenodes(1))*grad_lam(facenodes(2),1:3)) &
            - 2.0*DOT_PRODUCT(z_cross_grad_lam(facenodes(1),1:3),                                                &
            simp(facenodes(2))*grad_lam(facenodes(3),1:3) + simp(facenodes(3))*grad_lam(facenodes(2),1:3))
    CASE DEFAULT
       STOP 'IE: element type not supported in POINT_EVALUATE_FACE_FUNCTIONS.'
    END SELECT

    ! Evaluate E3:
    DO count1 = 9,11 ! 3 edges
       IF (num_vbfs.LT.count1) RETURN
       tempnodes = LOCAL_EDGENODES(tempedges(count1-8))
       SELECT CASE(ELEMENT_TYPE)
       CASE(3) ! Webb99
          tempvec =   simp(tempnodes(2))*(2.0*simp(tempnodes(1))-simp(tempnodes(2)))*grad_lam(tempnodes(1),1:3) &
               + simp(tempnodes(1))*(simp(tempnodes(1))-2.0*simp(tempnodes(2)))*grad_lam(tempnodes(2),1:3)
          term2(count1,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
          term1(count1) = &
               2.0*DOT_PRODUCT(z_cross_grad_lam(tempnodes(1),1:3),                                                &
               simp(tempnodes(1))*grad_lam(tempnodes(2),1:3) + simp(tempnodes(2))*grad_lam(tempnodes(1),1:3)  &
               - simp(tempnodes(2))*grad_lam(tempnodes(2),1:3))                                               &
               + 2.0*DOT_PRODUCT(z_cross_grad_lam(tempnodes(2),1:3),                                                &
               simp(tempnodes(1))*grad_lam(tempnodes(1),1:3) - simp(tempnodes(2))*grad_lam(tempnodes(1),1:3)  &
               - simp(tempnodes(1))*grad_lam(tempnodes(2),1:3))
       CASE DEFAULT
          STOP 'IE: element type not supported in POINT_EVALUATE_FACE_FUNCTIONS.'
       END SELECT
    END DO

    ! Evaluate F3:
    IF (num_vbfs.LT.12) RETURN
    SELECT CASE(ELEMENT_TYPE)
    CASE(3) ! Webb99
       tempvec =   simp(facenodes(1))*simp(facenodes(2))*grad_lam(facenodes(3),1:3)  &
            + simp(facenodes(2))*simp(facenodes(3))*grad_lam(facenodes(1),1:3)  &
            + simp(facenodes(3))*simp(facenodes(1))*grad_lam(facenodes(2),1:3)
       term2(12,1:3) = (/ -tempvec(2), tempvec(1), 0.0 /) ! Carry out the ^z cross-operation
       term1(12)     = &
            DOT_PRODUCT(z_cross_grad_lam(facenodes(3),1:3),                                                &
            simp(facenodes(1))*grad_lam(facenodes(2),1:3) + simp(facenodes(2))*grad_lam(facenodes(1),1:3)) &
            + DOT_PRODUCT(z_cross_grad_lam(facenodes(2),1:3),                                                &
            simp(facenodes(1))*grad_lam(facenodes(3),1:3) + simp(facenodes(3))*grad_lam(facenodes(1),1:3)) &
            + DOT_PRODUCT(z_cross_grad_lam(facenodes(1),1:3),                                                &
            simp(facenodes(2))*grad_lam(facenodes(3),1:3) + simp(facenodes(3))*grad_lam(facenodes(2),1:3))
    CASE DEFAULT
       STOP 'IE: element type not supported in POINT_EVALUATE_FACE_FUNCTIONS.'
    END SELECT

    IF (num_vbfs.GT.12) STOP 'IE: Higher order VBFs not implemented in POINT_EVALUATE_FACE_FUNCTIONS.'

  END SUBROUTINE POINT_EVALUATE_FACE_FUNCTIONS
  !*******************************************************************************







  FUNCTION VBF(vbf_type,lambda,grad_lambda,node1,node2,node3)
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! Vector basis function for node1 and node2 evaluated at simplex coordinates
    ! (lambda(1),lambda(2),lambda(3),lambda(4)).
    ! Webb99 elements added 24 Feb 02 DBD. 
    ! Note that scaling by the edge length (if appropriate) is done elsewhere. 
    ! The numbering conventions for the type of basis function is:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2. (Up to and including LT/QN)
    ! Future extensions will include:
    ! 5 -> e3, 6 -> f3 (Up to and including QT/QN)
    ! "Dummy" code at present to return zeros.

    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::   vbf_type 
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
    ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda 
    ! Gradient of simplex coordinates
    INTEGER(I4B), INTENT(IN) ::   node1,node2
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 ! Only required for face-based
    ! functions.
    REAL(SP), DIMENSION (3) :: VBF                   

    INTEGER(i4b) :: tmpi


    IF (vbf_type.NE.1.AND.vbf_type.NE.2.AND.vbf_type.NE.5.AND..NOT.PRESENT(node3)) THEN
       STOP 'IE: Internal error in FUNCTION VBF, third node required'
    END IF


!!$    PRINT*, 'lambda: ', lambda
!!$
!!$    PRINT*, 'grad_lambda: '
!!$    DO tmpi = 1, SIZE(grad_lambda, 1)
!!$       PRINT*, grad_lambda(tmpi, :)
!!$    END DO
!!$
!!$    PRINT*, 'node1: ', node1, ' node2: ', node2


    SELECT CASE(vbf_type) 
    CASE (1) ! e1 
       VBF(1:3) = lambda(node1)*grad_lambda(node2,1:3) - & 
            lambda(node2)*grad_lambda(node1,1:3)
    CASE (2) ! e2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,3)  ! Savage98 and Webb99 elements. 
          VBF(1:3) = lambda(node1)*grad_lambda(node2,1:3) + & 
               lambda(node2)*grad_lambda(node1,1:3)
       CASE (2)    ! Andersen and Volakis elements. 
          VBF(1:3) = ( lambda(node1) - lambda(node2) ) * & 
               ( lambda(node1)*grad_lambda(node2,1:3) - & 
               lambda(node2)*grad_lambda(node1,1:3) )
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE (3)! f1
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2)  ! Savage98 and A&V elements. 
          VBF(1:3) =   lambda(node1)* & 
               ( lambda(node2)*grad_lambda(node3,1:3) - & 
               lambda(node3)*grad_lambda(node2,1:3) )
       CASE (3)  ! Webb99 elements.
          VBF(1:3) =   lambda(node2)*lambda(node3)*grad_lambda(node1,1:3) & 
               +lambda(node1)*lambda(node3)*grad_lambda(node2,1:3) & 
	       -2.0_SP*lambda(node1)*lambda(node2)*grad_lambda(node3,1:3) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE (4) !f2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2)  ! Savage98 and A&V elements. 
          VBF(1:3) =   lambda(node2)* & 
               ( lambda(node1)*grad_lambda(node3,1:3) - & 
               lambda(node3)*grad_lambda(node1,1:3) )
       CASE (3)  ! Webb99 elements.
          VBF(1:3) =   lambda(node3)*lambda(node1)*grad_lambda(node2,1:3) & 
               +lambda(node2)*lambda(node1)*grad_lambda(node3,1:3) & 
	       -2.0_SP*lambda(node2)*lambda(node3)*grad_lambda(node1,1:3)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE (5) ! e3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only. 
          VBF(1:3) = lambda(node2)*(2.0_SP*lambda(node1)-lambda(node2))* grad_lambda(node1,1:3) &
               - lambda(node1)*(2.0_SP*lambda(node2)-lambda(node1))* grad_lambda(node2,1:3)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE (6) ! f3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only.
          VBF(1:3) =   lambda(node2)*lambda(node3)*grad_lambda(node1,1:3) & 
               +lambda(node1)*lambda(node3)*grad_lambda(node2,1:3) & 
               +lambda(node1)*lambda(node2)*grad_lambda(node3,1:3) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE DEFAULT
       STOP 'IE: Invalid vector basis function type'
    END SELECT
  END FUNCTION VBF
  !*******************************************************************************


  FUNCTION VBF_S(vbf_type,lambda,grad_lambda,normal,node1,node2,node3)
    USE math_tools, ONLY: CROSS_PRODUCT
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! Surface vector basis function for node1 and node2 evaluated at simplex coordinates
    ! (lambda(1),lambda(2),lambda(3),lambda(4)). 
    ! This function computes the surface vbf's by taking the surface normal 
    ! cross producted with the gradient of the (3D) simplex coordinates; 
    ! it has been shown that 
    ! this projection into  the plane of a triangular element is identical
    ! to the gradient of the equivalent 2D simplex coordinates. 
    ! [Lab Notebook, DB Davidson, p.40-41]
    !
    ! Webb99 elements added 24 Feb 02 DBD. 
    ! Note that scaling by the edge length (if appropriate) is done elsewhere. 
    !
    ! Future extensions will include:
    ! 5 -> e3, 6 -> f3 (Up to and including QT/QN)
    ! "Dummy" code at present to return zeros.

    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) ::   vbf_type 
    ! Type of vector basis function:
    ! 1 -> e1, 2 -> e2, 3 -> f1, 4->f2.
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
    ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda 
    ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(3), INTENT(IN) :: normal ! Outward unit normal on surface
    INTEGER(I4B), INTENT(IN) ::   node1,node2
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  node3 ! Only required for face-based
    ! functions.
    REAL(SP), DIMENSION (3) :: VBF_S
    REAL(SP), DIMENSION(4,3) :: nXgrad_lambda 
    ! normal cross gradient of simplex coordinates
    INTEGER(I4B) :: count ! counter

    IF (vbf_type.NE.1.AND.vbf_type.NE.2.AND.vbf_type.NE.5.AND..NOT.PRESENT(node3)) THEN
       STOP 'IE: Internal error in FUNCTION VBF_S, third node required'
    END IF

    DO count = 1,4 ! Four simplex coordinates
       nXgrad_lambda(count,1:3) = CROSS_PRODUCT(normal,grad_lambda(count,1:3))
    END DO

    SELECT CASE(vbf_type) 
    CASE (1) ! e1
       VBF_S(1:3) = lambda(node1)*nXgrad_lambda(node2,1:3) - & 
            lambda(node2)*nXgrad_lambda(node1,1:3)
    CASE (2) ! e2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,3) ! Savage98 and Webb99 elements.   
          VBF_S(1:3) = lambda(node1)*nXgrad_lambda(node2,1:3) + & 
               lambda(node2)*nXgrad_lambda(node1,1:3)
       CASE (2)   ! Andersen & Volakis elements.
          VBF_S(1:3) = ( lambda(node1) - lambda(node2) ) * & 
               ( lambda(node1)*nXgrad_lambda(node2,1:3) - & 
               lambda(node2)*nXgrad_lambda(node1,1:3) )

       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF_S'
       END SELECT
    CASE (3)! f1
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2)  ! Savage98 and A&V elements. 
          VBF_S(1:3) = ( lambda(node1)) * & 
               ( lambda(node2)*nXgrad_lambda(node3,1:3) - & 
               lambda(node3)*nXgrad_lambda(node2,1:3) )
       CASE (3)  ! Webb99 elements.
          VBF_S(1:3) = ( lambda(node2)*lambda(node3)*nXgrad_lambda(node1,1:3) & 
               +lambda(node1)*lambda(node3)*nXgrad_lambda(node2,1:3) & 
	       -2.0_SP*lambda(node1)*lambda(node2)*nXgrad_lambda(node3,1:3) )
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF_S'
       END SELECT
    CASE (4) !f2
       SELECT CASE(ELEMENT_TYPE) 
       CASE (1,2)  ! Savage98 and A&V elements. 
          VBF_S(1:3) = ( lambda(node2)) * & 
               ( lambda(node1)*nXgrad_lambda(node3,1:3) - & 
               lambda(node3)*nXgrad_lambda(node1,1:3) )
       CASE (3)  ! Webb99 elements.
          VBF_S(1:3) = ( lambda(node3)*lambda(node1)*nXgrad_lambda(node2,1:3) & 
               +lambda(node2)*lambda(node1)*nXgrad_lambda(node3,1:3) & 
	       -2.0_SP*lambda(node2)*lambda(node3)*nXgrad_lambda(node1,1:3) )
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF_S'
       END SELECT
    CASE (5) ! e3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only.
          VBF_S(1:3) = lambda(node2)*(2.0_SP*lambda(node1)-lambda(node2))* & 
               nXgrad_lambda(node1,1:3) &
               - lambda(node1)*(2.0_SP*lambda(node2)-lambda(node1))* &
               nXgrad_lambda(node2,1:3)
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE (6) ! f3
       SELECT CASE(ELEMENT_TYPE) 
       CASE (3)  ! Webb99 elements only.
          VBF_S(1:3) = lambda(node2)*lambda(node3)*nXgrad_lambda(node1,1:3) & 
               +lambda(node1)*lambda(node3)*nXgrad_lambda(node2,1:3) & 
               +lambda(node1)*lambda(node2)*nXgrad_lambda(node3,1:3) 
       CASE DEFAULT
          STOP 'IE: Invalid element type in FUNCTION VBF'
       END SELECT
    CASE DEFAULT
       STOP 'IE: Invalid vector basis function type'
    END SELECT
  END FUNCTION VBF_S
  !*******************************************************************************

  FUNCTION PARENT_VBF(vbf_type,ii,u,v,w)
    USE nrtype
    USE problem_info
    IMPLICIT NONE
    !*******************************************************************************
    ! The vector basis functions on a parent "unitary" tetrahedron. 
    ! Edge length scaling is not supported.
    ! DBD, 04 August 05. Implemented only for Webb99 basis functions. 
    !*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: vbf_type ! 
    INTEGER(I4B), INTENT(IN) :: ii       ! 
    REAL(SP),INTENT(IN) :: u,v,w         ! Local coordinates on parent triangle (associated with u1,u2,u3 or x,y,z or xi,eta,zeta etc)  
    REAL(SP), DIMENSION (3) :: PARENT_VBF

    IF (SCALE_BY_EDGE_LENGTH) THEN
       STOP 'IE: FUNCTION PARENT_VBF does not support this element type with edge length scaling'
    END IF


    SELECT CASE(vbf_type) 
    CASE (1) ! e1
       SELECT CASE(ii) ! edge number
       CASE (1) 
          PARENT_VBF = (/  1.0_SP-v-w, u, u  /)
       CASE (2) 
          PARENT_VBF = (/  v, 1.0_SP-u-w , v /)
       CASE (3) 
          PARENT_VBF = (/  w , w , 1.0_SP-u-v/)
       CASE (4) 
          PARENT_VBF = (/  -v , u , 0.0_SP   /)
       CASE (5) 
          PARENT_VBF = (/  -w , 0.0_SP, u   /)
       CASE (6) 
          PARENT_VBF = (/  0.0_SP, -w, v   /)
       END SELECT
    CASE (2) ! e2
       SELECT CASE(ii) 
       CASE (1) 
          PARENT_VBF = (/  1.0_SP-2.0_SP*u-v-w, -u, -u  /)
       CASE (2) 
          PARENT_VBF = (/  -v, 1.0_SP-u-2.0_SP*v-w , -v /)
       CASE (3) 
          PARENT_VBF = (/  -w , -w , 1.0_SP-u-v-2.0_SP*w/)
       CASE (4) 
          PARENT_VBF = (/   v , u , 0.0_SP   /)
       CASE (5) 
          PARENT_VBF = (/   w , 0.0_SP, u   /)
       CASE (6) 
          PARENT_VBF = (/  0.0_SP,  w, v   /)
       END SELECT
    CASE (3)! f1
       IF(ELEMENT_TYPE.NE.3) THEN
          STOP 'IE: Only Webb99-type basis functions supported in FUNCTION PARENT_VBF'
       END IF
       SELECT CASE(ii) ! face number
       CASE (1) 
          PARENT_VBF = (/ v-2*u*v-v**2-v*w , 2*u**2-2*u+u*v+2*u*w, -u*v /)
       CASE (2) 
          PARENT_VBF = (/ -2*u*w+w-v*w-w**2, -u*w, -2*u+2*u**2+2*u*v+u*w/)
       CASE (3) 
          PARENT_VBF = (/ -v*w, w-u*w-2*v*w-w**2, -2*v+2*u*v+2*v**2+v*w/)
       CASE (4) 
          PARENT_VBF = (/ v*w, u*w, -2*u*v/)
       END SELECT
    CASE (4)! f2
       IF(ELEMENT_TYPE.NE.3) THEN
          STOP 'IE: Only Webb99-type basis functions supported in FUNCTION PARENT_VBF'
       END IF
       SELECT CASE(ii) ! face number
       CASE (1) 
          PARENT_VBF = (/ v+u*v-v**2-v*w , u - u**2 + u*v-u*w, +2*u*v /)
       CASE (2) 
          PARENT_VBF = (/ w+u*w-v*w-w**2, 2*u*w, u-u**2-u*v+u*w/)
       CASE (3) 
          PARENT_VBF = (/ 2*v*w, w-u*w+v*w-w**2, v-u*v-v**2+v*w/)
       CASE (4) 
          PARENT_VBF = (/ -2*v*w, u*w, +u*v/)
       CASE DEFAULT
          STOP 'IE: Invalid face number in FUNCTION PARENT_VBF'
       END SELECT
    CASE (5) ! e3
       STOP 'UNIMPLEMTED IN PARENT_VBF'
    CASE (6) ! f3
       STOP 'UNIMPLEMTED IN PARENT_VBF'
    CASE DEFAULT
       STOP 'IE: Invalid or unimplemented vector basis function type in FUNCTION PARENT_VBF'
    END SELECT
  END FUNCTION PARENT_VBF
  !*******************************************************************************



END MODULE basis_function

