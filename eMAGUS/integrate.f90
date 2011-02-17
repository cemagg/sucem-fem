MODULE integrate
  USE nrtype
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: integrate_tet

CONTAINS

  FUNCTION integrate_tet(int_func, order)
    USE quad_tables
    IMPLICIT NONE
!!!
!!! Interface Variables/Functions
!!!
    REAL(SP) :: integrate_tet
    INTERFACE
       FUNCTION int_func(lambda)
         USE nrtype
         REAL(SP), DIMENSION(4), INTENT(in) :: lambda
         REAL(SP) :: int_func
       END FUNCTION int_func
    END INTERFACE
    INTEGER(i4b), INTENT(in) :: order ! Integration order
!!!
!!! Local Variables
!!! 
    INTEGER(i4b) :: num_cpts    ! Number of cubature points to use
    INTEGER(i4b) :: i_point     ! Integration point loop var
    REAL(SP) :: weight          ! Weight of an integration  point
    REAL(SP), DIMENSION(4) :: lambda ! Integration point simplex co-ordinates
    REAL(DP) :: accumulator
!!!
    SELECT CASE(order)
    CASE(1:2)
       num_cpts = 4
    CASE(3:4)
       num_cpts = 11
    CASE(5:6)
       num_cpts = 24
    CASE default
       PRINT*, "Unimplemented cubature order ", order, &
            " in module integrate, function integrate_tet"
       STOP 
    END SELECT

    accumulator = 0             ! Initialise the integral to 0

    DO i_point = 1, num_cpts
       lambda = cube_tet_rules(num_cpts)%rule(i_point,1:4)
       weight = cube_tet_rules(num_cpts)%rule(i_point,5)
       accumulator = accumulator + weight * int_func(lambda)
    END DO

    integrate_tet = accumulator
    
  END FUNCTION integrate_tet


END MODULE integrate
