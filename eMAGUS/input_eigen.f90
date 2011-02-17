!*******************************************************************************
! MODULE input_eigen : Subroutines and related data structures for 
!                      reading the eigensolution control elements from
!                      the input file namelist
!*******************************************************************************

MODULE input_eigen
  USE nrtype
  USE problem_info, ONLY: OPTIONLENGTH
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_eigen_namelist, test_input_eigen_print
  
  CHARACTER(OPTIONLENGTH) :: first_mode_select
  INTEGER(I4B) :: mode_index, num_modes
  REAL(SP) :: eigen_estimate

  NAMELIST/EIGEN/ first_mode_select, mode_index, num_modes, eigen_estimate

CONTAINS

!*******************************************************************************
!   SUBROUTINE read_eigen_namelist(infileunit)
!
! Reads the eigen namelist from unit infileunit, and then
! puts the data in the eigen_analysis_data module structures
!*******************************************************************************

  SUBROUTINE read_eigen_namelist(infileunit)
    USE namelist_tools
    USE eigen_analysis_data
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios
    
    CALL eigen_input_defaults
    CALL eigen_data_defaults
    
    READ(unit=infileunit, nml=eigen,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)
    
    first_eigenmode = mode_index 
    first_eigen_approx = eigen_estimate
    SELECT CASE(first_mode_select)
    CASE('estimate')
       COMPUTE_EIGENMODES = .TRUE.
    CASE('predict')
       COMPUTE_EIGENMODES = .TRUE.
       PREDICT_SPURIOUS_EIGENMODES = .TRUE.
    CASE('none')
       COMPUTE_EIGENMODES = .FALSE.
       PREDICT_SPURIOUS_EIGENMODES = .FALSE.
    CASE default
       PRINT*, "Invalid value &EIGEN first_mode_select=", TRIM(first_mode_select)
       STOP
    END SELECT

    CALL eigen_data_check
    
  END SUBROUTINE read_eigen_namelist

  SUBROUTINE eigen_input_defaults
    first_mode_select = 'none'
    num_modes = 0
    mode_index = 1
    eigen_estimate = 0.0_SP
  END SUBROUTINE eigen_input_defaults

  SUBROUTINE test_input_eigen_print
    IMPLICIT NONE
    
    PRINT*,'Namelist: eigen'
    PRINT*,'first_mode_select     : ', TRIM(first_mode_select)
    PRINT*,'mode_index            : ', mode_index
    PRINT*,'num_modes             : ', num_modes
    PRINT*,'eigen_estimate        : ', eigen_estimate
END SUBROUTINE test_input_eigen_print

END MODULE input_eigen
