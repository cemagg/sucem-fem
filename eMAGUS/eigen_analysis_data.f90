MODULE eigen_analysis_data
  USE nrtype
  USE problem_info
  SAVE 

  INTEGER(I4B) first_eigenmode, num_eigenmodes, num_spurious_eigenmodes
                                             ! Self-explanatory
  REAL(SP) first_eigen_approx                ! Approximation of eigenvalue
                                             ! (alternate to specifying number
                                             ! of eigenvalue or estimate.
  LOGICAL(LGT) PREDICT_SPURIOUS_EIGENMODES   ! Flag to trigger theoretical
                                             ! prediction.
CONTAINS
  
  SUBROUTINE eigen_data_defaults
    IMPLICIT NONE
    PREDICT_SPURIOUS_EIGENMODES = .FALSE.
    first_eigenmode = 1
    num_eigenmodes = 1
    first_eigen_approx = 0.0_SP
  END SUBROUTINE eigen_data_defaults

  SUBROUTINE eigen_data_check
    IMPLICIT NONE
    IF (first_eigenmode < 0) THEN
       PRINT*, &
            "Invalid negative value of FIRST_EIGENMODE in module eigen_analysis_data"
       STOP
    END IF

    IF (COMPUTE_EIGENMODES .AND. .NOT.(PREDICT_SPURIOUS_EIGENMODES)) THEN
       IF (FIRST_EIGEN_APPROX <= 0) THEN
          PRINT*, &
               "Invalid eigenvalue approximation specified in module eigen_analysis_data"
          STOP
       END IF
    END IF
    
  END SUBROUTINE eigen_data_check
    
END MODULE eigen_analysis_data
!***********************************************************************
