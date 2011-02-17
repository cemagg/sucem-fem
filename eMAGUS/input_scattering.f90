!*******************************************************************************
! MODULE input_scattering : Subroutines and related data structures for 
!                      reading the scattering control elements from
!                      the input file namelist
!*******************************************************************************

MODULE input_scattering
  USE nrtype
  USE problem_info, ONLY: OPTIONLENGTH
  USE scattering_analysis_data
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_scattering_namelist, test_input_scattering_print
  
  CHARACTER(OPTIONLENGTH) :: scat_type

  NAMELIST/SCATTER/ SCAT_FIELD, TEST_INTERPOLATE_FIELD, scat_type, &
       sph_radius

CONTAINS

!*******************************************************************************
!   SUBROUTINE read_scattering_namelist(infileunit)
!
! Reads the scattering namelist from unit infileunit, and then
! puts the data in the scattering_analysis_data module structures
!*******************************************************************************

  SUBROUTINE read_scattering_namelist(infileunit)
    USE namelist_tools
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios
    
    CALL scattering_input_defaults
    CALL scattering_data_defaults
    
    READ(unit=infileunit, nml=scatter,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)

    SELECT CASE(scat_type)
    CASE('pec_sphere')  
       PEC_SPHERE_SCAT = .TRUE.
    CASE('pec_general') 
       PEC_GENERAL_SCAT = .TRUE.
    CASE('penetrable')  
       PENETRABLE_SCAT = .TRUE.
    CASE default
       PRINT*, "Invalid &scatter scat_type=", TRIM(scat_type)
    END SELECT

    CALL test_input_scattering_print
    CALL scattering_data_check
  END SUBROUTINE read_scattering_namelist

  SUBROUTINE scattering_input_defaults
    scat_type = 'none'
  END SUBROUTINE scattering_input_defaults

  SUBROUTINE test_input_scattering_print
    IMPLICIT NONE
    
    PRINT*,'Namelist: scattering'
    PRINT*,'SCAT_FIELD                 : ', SCAT_FIELD
    PRINT*,'TEST_INTERPOLATE_FIELD     : ', TEST_INTERPOLATE_FIELD
    PRINT*,'PEC_SPHERE_SCAT            : ', PEC_SPHERE_SCAT
    PRINT*,'PEC_GENERAL_SCAT           : ', PEC_GENERAL_SCAT
    PRINT*,'PENETRABLE_SCAT            : ', PENETRABLE_SCAT
    PRINT*,'sph_radius                 : ', sph_radius
END SUBROUTINE test_input_scattering_print

END MODULE input_scattering
