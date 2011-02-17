!*******************************************************************************
! MODULE input_problem_info : Subroutines and related data structures for 
!                             reading the general problem  control elements from
!                             the input file namelist
!*******************************************************************************
MODULE input_problem_info
  USE nrtype
  USE problem_info
  USE quad_tables, ONLY: gauss_points
  USE matrix, ONLY: residual_norm, max_iter_factor, restart_GMRES,       &
       matrix_defaults
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_problem_namelists

  CHARACTER (OPTIONLENGTH) :: solution_type, element_order, solver
  INTEGER(i4b) :: tri_intg_order

  NAMELIST/PROBLEM/ solution_type, element_order, EXECUTE, CUBATURE,     &
       OUTPUT_ELEMENT_DATA, ON_SCREEN_REPORTING,             &
       SPARSE, USE_PRE_CONDITIONER, solver, OUTPUT_ELEMENT_DATA,         &
       OUTPUT_ELEMENT_SHAPE, tri_intg_order

  CHARACTER(OPTIONLENGTH) :: elemental_modeling

  NAMELIST/GEOMETRICAL/ elemental_modeling

  NAMELIST/MATRIX_SOLUTION/ residual_norm, max_iter_factor, restart_GMRES

CONTAINS
  
  SUBROUTINE read_problem_namelists(infileunit)
    USE namelist_tools
    USE input_eigen
    USE input_scattering
    USE input_frequency

    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios
    LOGICAL(LGT) :: frequency_domain ! Flag that is set if frequency info
                                     ! is to be read
    frequency_domain = .FALSE.
    
    CALL problem_input_defaults ! Set up defaults for this module
    CALL problem_info_defaults  ! Set up default values for module problem_info
    CALL matrix_defaults        ! Set up default values for module matrix
    
    READ(unit=infileunit, nml=problem,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)

    IF (SPARSE) THEN 
       SOLVER_TYPE = 1          ! Default to BICG for sparse matrices
    ELSE
       SOLVER_TYPE = 0          ! LU Decomp for full matrices
    END IF

    READ(unit=infileunit, nml=matrix_solution,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)
    
    READ(unit=infileunit, nml=geometrical,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)

    SELECT CASE(solution_type)
    CASE('real eigen')
       REAL_EIGEN_ANALYSIS = .TRUE.
       analysis_counter = 1
       CALL read_eigen_namelist(infileunit)
    CASE('complex eigen')
       CMPLX_EIGEN_ANALYSIS = .TRUE. 
       analysis_counter = 1
       CALL read_eigen_namelist(infileunit)
    CASE('cbaa')
       CBAA_ANALYSIS = .TRUE.
       frequency_domain = .TRUE.
    CASE('guided wave')
       GW_ANALYSIS = .TRUE.
       frequency_domain = .TRUE.
    CASE('fd scattering')
       FD_SCAT_ANALYSIS = .TRUE.
       CALL read_scattering_namelist(infileunit)
       frequency_domain = .TRUE.
    CASE('time domain')
       TD_ANALYSIS = .TRUE.
    CASE('coupled time domain')
       COUPLED_TD_ANALYSIS = .TRUE.
    CASE('write mesh')
       WRITE_MESH = .TRUE.
    CASE default
       PRINT*, 'Unknown solution type ', solution_type
       STOP
    END SELECT
    ! Read the frequency namelist for FD analyses
    IF (frequency_domain) CALL read_frequency_namelist(infileunit)

    SELECT CASE(solver)
    CASE('LU')
       SOLVER_TYPE = 0
    CASE('BICG')
       SOLVER_TYPE = 1
    CASE('CG')
       SOLVER_TYPE = 2
    CASE('QMR')
       SOLVER_TYPE = 3
    CASE('GMRES')
       SOLVER_TYPE = 4
    CASE('Experimental')
       SOLVER_TYPE = 5
    CASE default
       PRINT*, 'Unknown Solver Type, ', solver
       STOP
    END SELECT
    

    SELECT CASE(element_order)
    CASE('CTLN') 
       HIERARCHAL_ORDER = 1
       MIXED_ORDER_FLAG = .TRUE.
    CASE('LTLN') 
       HIERARCHAL_ORDER = 1
       MIXED_ORDER_FLAG = .FALSE.
    CASE('LTQN') 
       HIERARCHAL_ORDER = 2 
       MIXED_ORDER_FLAG = .TRUE.
    CASE('QTQN') 
       HIERARCHAL_ORDER = 2 
       MIXED_ORDER_FLAG = .FALSE.
    CASE default
       PRINT*, 'Unknown element order ', element_order
       STOP
    END SELECT

    SELECT CASE(elemental_modeling)
    CASE('simplex')
       ! Do nothing, since CURVILINEAR and ALL_CURVILINEAR are initialised to false
    CASE('curvilinear')
       CURVILINEAR = .TRUE.
    CASE('all curvilinear')
       CURVILINEAR = .TRUE.
       ALL_CURVILINEAR = .TRUE.
    CASE default
       PRINT*, 'Unknown elemental_modeling ', elemental_modeling
    END SELECT

    SELECT CASE(tri_intg_order)
    CASE(1)
       gauss_points = 1
    CASE(2)
       gauss_points = 3
    CASE(3)
       gauss_points = 4
    CASE(4)
       gauss_points = 6
    CASE(5)
       gauss_points = 7
    CASE(7)
       gauss_points = 13
    CASE(10)
       gauss_points = 25
    CASE default
       PRINT*, 'Unimplemented triangle integration order ', tri_intg_order
       STOP
    END SELECT
       
  END SUBROUTINE read_problem_namelists

  SUBROUTINE problem_input_defaults

    solution_type = 'none'
    element_order = 'none'
    elemental_modeling = 'simplex'
    tri_intg_order = 4          ! Default to 4th order triangle integration
    solver = 'CG'               ! Default to using the Conjugate Gradient solver
    element_order = 'CTLN'      ! Default to using CT/LN elements
  END SUBROUTINE problem_input_defaults

END MODULE input_problem_info
