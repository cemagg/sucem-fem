!**********************************************************************
!* 
!* MODULE twoform_basisfun
!* 
!* Description:
!* 
!*   Routines that deal with two-form (div-conforming) basis functions.
!*   At present, only the lowest order face functions are 
!*   implemented. 
!*
!*   3-D two-form basis functions have DOFs associated with faces and
!*   volumes (unlike one-forms, that have edge-related DOFs too).
!*
!* Interface Design:
!*
!*   The user calls functions, that return the value of the given basis
!*   function at the given point. The functions are named as:
!*
!*   (operator)_two_vbf_(DOF geom assoc)
!*   
!*  eg. curl_two_vbf_face() would evaluate the curl of the two-form 
!*  face-associated basis function.
!* 
!**********************************************************************


MODULE twoform_vbf
  USE nrtype
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: twoform_vbf_face, convert_gradlambda

CONTAINS

  !**********************************************************************
  !* 
  !* FUNCTION convert_gradlambda(grad_lambda, facenodes)
  !* 
  !* Input parameters:
  !*   
  !*   grad_lambda, real, dimension(4,3)
  !*                     -> The gradients of the four simplex co-ordinates
  !*   facenodes, integer, dimension(3)
  !*                     -> Element local numbers that make up the face
  !*   
  !* Output:
  !*   
  !*   Returns the cross-products of pairs of simplex co-ordinates for
  !*   use as a covariant base for a single face, within a given element.
  !*   
  !* Description:
  !*   
  !*   The elemental properties are contained in grad_lambda, which is constant
  !*   within a given element. Supposing a face is defined by nodes i, j, k,
  !*   convert_gradlambda returns:
  !*                              Denotes Cross Product
  !*                                     |
  !*   convert_gradlambda(1,:) = grad(lambda_j) X grad(lambda_k)
  !*   convert_gradlambda(2,:) = grad(lambda_k) X grad(lambda_i)
  !*   convert_gradlambda(3,:) = grad(lambda_i) X grad(lambda_j)
  !*   
  !**********************************************************************
  
  FUNCTION convert_gradlambda(grad_lambda, facenodes)
    USE math_tools, ONLY: CROSS_PRODUCT
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda
    INTEGER(i4b), DIMENSION(3), INTENT(in) :: facenodes
    REAL(SP), DIMENSION(3,3) :: convert_gradlambda
!!!
!!! Local Variables
!!!
    INTEGER(i4b) :: i,a,b       ! Loop var
    ! Temporar gradients variable
    REAL(SP), DIMENSION(3) :: grad_a, grad_b
!!!
    DO i=1,3                    ! Loop over facenodes
       a = MODULO(i+1,3) 
       b = MODULO(a+1,3)
       IF (a == 0)  a = 3
       IF (b == 0)  b = 3
       grad_a = grad_lambda(a,:) 
       grad_b = grad_lambda(b,:)
       convert_gradlambda(i, :) = CROSS_PRODUCT(grad_a, grad_b)
    END DO

  END FUNCTION convert_gradlambda


  !**********************************************************************
  !* 
  !* FUNCTION twoform_vbf_face(vbf_type, lambda, grad_lambda, facenodes)
  !* 
  !* Input parameters:
  !*   
  !*   vbf_type, integer -> partial "order" of face func, for now only
  !*                        1 -> LT/CN face func
  !*   lambda, real(SP), dimension (4)
  !*                     -> The values of the 4 simplex co-ords
  !*   grad_lambda, real(SP), dimension(4,3)
  !*                     -> the covariant basis to use for basis function 
  !*                        expansion
  !*   facenodes, integer, dimension(3)
  !*                     -> The three (element local) node numbers of
  !*                        the face under consideration                  
  !* 
  !* Output:
  !* 
  !*   The value of the basis function at the requested point
  !* 
  !* Description:
  !* 
  !*   This function is just a higher-level case-statement that calls
  !*   the individual basis function for the appropriate order of face
  !*   function. It also does not check that the values of lambda are
  !*   valid (eg. outside element, not adding up to 1, etc.)
  !* 
  !*   Note that for rectilinear elements, unscaled_cov_basis is constant
  !*   throughout the whole element, and is calculated by taking the 
  !*   cross product of the gradients of the simplex co-ordinates. See,
  !*   eg. the face basis functions defined in [Lee97] sec IV. For a more
  !*   general overview, see [botha05].
  !*   
  !*   The covariant basis function is considered to be unscaled, since 
  !*   the defintion of of covariant vectors require the cross-product
  !*   of the contravariant vectors to be normalised by the Jacobian.
  !*   In [botha05]'s notation, we are dealing with the f_i vectors.
  !*   
  !**********************************************************************
  
  FUNCTION twoform_vbf_face(vbf_type, lambda, grad_lambda, facenodes)
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    ! Partial "order" of the face function
    INTEGER(I4B), INTENT(IN) ::   vbf_type 
    ! Simplex coordinates of point.
    REAL(SP), DIMENSION(4), INTENT(IN) :: lambda   
    ! Gradient of simplex coordinates
    REAL(SP), DIMENSION(4,3),INTENT(IN) :: grad_lambda
    INTEGER(i4b), INTENT(in), DIMENSION(3) :: facenodes
    REAL(SP), DIMENSION(3) :: twoform_vbf_face
!!!
!!! Local Variables
!!!
    ! Face Local simplex co-ords
    REAL(SP), DIMENSION(3) :: lambda_face
    ! Face Local basis
    REAL(SP), DIMENSION(3,3) :: basis_face
    INTEGER(i4b) :: i,j         ! Loop variables
!!!
    lambda_face = lambda(facenodes) ! Assign face-local co-ordinates
    basis_face = convert_gradlambda(grad_lambda, facenodes) ! Get face local basis vectors.
    SELECT CASE(vbf_type)
    CASE(1)
       twoform_vbf_face = two_vbf_face_1(lambda_face, basis_face)
    CASE default
       PRINT*, "Unknown face order vbf_type = ", vbf_type, " in MODULE twoform_basisfun, FUNCTION two_vbf_face"
       STOP 
    END SELECT
  END FUNCTION twoform_vbf_face

  !**********************************************************************
  !* 
  !* FUNCTION two_vbf_face_1(lambda, basis)
  !* 
  !* Input parameters:
  !*   
  !*   lambda, real(SP), dimension (3)                                    
  !*                     -> The values of the 3 face-local simplex co-ords           
  !*   basis, real(SP), dimension(3,3)                       
  !*                     -> the covariant basis to use for basis function 
  !*                        expansion                                     
  !* Output:
  !*   
  !*   The 2-form face function evaluated at lambda, as described in, eg.
  !*   [Lee97], sec IV.
  !*   
  !* Description:
  !*   
  !*   
  !*   
  !**********************************************************************
  
  FUNCTION two_vbf_face_1(lambda, basis)
    IMPLICIT NONE
!!!
!!! Interface Variables
!!!
    ! Simplex coordinates of point, using face-local coordinate numbers.
    REAL(SP), DIMENSION(3), INTENT(IN) :: lambda   
    REAL(SP), DIMENSION(3,3),INTENT(IN) :: basis
    REAL(SP), DIMENSION(3) :: two_vbf_face_1
!!!
!!! Local Variables
!!!
    INTEGER(i4b) :: i           ! Loop variable
!!!
    two_vbf_face_1 = 0          ! Initialise to 0
    
    DO i=1,3
       two_vbf_face_1 = two_vbf_face_1 + lambda(i)*basis(i,:)
    END DO
    two_vbf_face_1 = two_vbf_face_1*2
  END FUNCTION two_vbf_face_1
  
END MODULE twoform_vbf
