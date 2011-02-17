MODULE problem_info
  USE nrtype
  USE debugvar
  IMPLICIT NONE
!***********************************************************************
! This module contains values that are constant for a specific problem.
! Created: 2001-09-21, MMB.
! Last changed: 2002-02-20, DBD.
! 2002-03-17: Added ADAPTIVE_ANALYSIS flag. MMB.
! 2002-05-09: Added MAX_PRESENT_... . MMB.
! 2002-03-24: Added UNBOUNDED_TD_ANALYSIS. DBD.
!***********************************************************************
  SAVE

  ! Flags for analysis type
  LOGICAL(LGT) REAL_EIGEN_ANALYSIS
  LOGICAL(LGT) CMPLX_EIGEN_ANALYSIS 
  LOGICAL(LGT) CBAA_ANALYSIS
  LOGICAL(LGT) GW_ANALYSIS
! Added DBD 24 Mar 03
  LOGICAL(LGT) TD_ANALYSIS
! End DBD 24 Mar 03
! Added DBD 16 Dec 03
  LOGICAL(LGT) FD_SCAT_ANALYSIS
! End added DBD 16 Dec 03
  LOGICAL(LGT) COUPLED_TD_ANALYSIS          ! Coupled time domain formulation
  LOGICAL(lgt) WRITE_MESH                   ! Write out mesh, no calculations
  INTEGER(I4B) HIERARCHAL_ORDER             ! Order of hierarchal element
! Added DBD 20 Feb 02
  LOGICAL(LGT) MIXED_ORDER_FLAG             ! Flag to use mixed order elements.
! End additions
  INTEGER(I4B) ELEMENT_TYPE                 ! Type of (higher-order) element.
  LOGICAL(LGT) CUBATURE                     ! Flag to compute elemental matrix
                                            ! contributions numercially.
  LOGICAL(LGT) CURVILINEAR                  ! Flag to use curvilinear elements.
  LOGICAL(LGT) ALL_CURVILINEAR              ! Flag to make ALL elements curvilinear.
  LOGICAL(LGT) BANDRENUM                    ! Flag to renumber mesh.
  LOGICAL(LGT) BANDRENUM_STORE              ! Flag to store the S and T matrices
                                            ! in banded format  
  LOGICAL(LGT) SCALE_BY_EDGE_LENGTH         ! Flag to scale elements by edge length
                                            ! also primarily for testing.
  LOGICAL(LGT) GW_COMPUTE_S_PARAMS          ! Flag to trigger computation of 
                                            ! either S parameters or 
                                            ! field (within guided wave module). 
  LOGICAL(LGT) PEC_ground                   ! Flag for PEC gnd @ z=0.
  LOGICAL(LGT) USE_PRE_CONDITIONER          ! Flag for pre-conditioner.



  ! Program version information:
  CHARACTER(24),  PARAMETER :: VERSION_NUMBER = '1.00--neilen-0.0' 
  CHARACTER(8),  PARAMETER :: REVISION_NUMBER = 'patch-14' ! 1-99 supported. 
  CHARACTER(8), PARAMETER :: LAST_CHANGED = '20050422'

  ! Following controls output data:
  LOGICAL(LGT) OUTPUT_ELEMENT_DATA
  LOGICAL(LGT)  OUTPUT_ELEMENT_SHAPE

  ! Following controls execution:
  LOGICAL(LGT) EXECUTE
   
  ! Following controls sparse solution:
  LOGICAL(LGT) SPARSE
   
  ! solver type:
  INTEGER(I4B) SOLVER_TYPE 
   
  ! and on-screen reporting:
  LOGICAL(LGT) ON_SCREEN_REPORTING
   
  LOGICAL(LGT) COMPUTE_EIGENMODES
  ! Flag indicating whether the Whitney elements approximation (Gong95)
  ! or the general, dominant mode coax approach must be used (this flag
  ! is set internally):
  LOGICAL(LGT) COAX_WHITNEY

  ! Flag indicating whether adaptive FEM related actions must be
  ! performed, as specified in the AD card (thus it indicates the
  ! presence of an AD card):
  LOGICAL(LGT) ADAPTIVE_ANALYSIS

  ! Different analysis can occur in one *.fek file, these valiables 
  ! specify which one is currently processed and their total:
  INTEGER(I4B) :: analysis_counter,NUM_analysis

  ! Follwing are assigned in subroutine MESHIN:
  INTEGER(I4B) MAX_NODES      ! Maximum number of nodes
  INTEGER(I4B) MAX_ELEMENTS   ! Maximum number of elements
  INTEGER(I4B) MAX_EDGES      ! Maximum number of edges.
  INTEGER(I4B) MAX_FACES      ! Maximum number of faces.   
  INTEGER(I4B) MAX_BOUNDARIES ! Maximum number of different BC's.
  INTEGER(I4B) MAX_PORTS      ! Maximum number of ports.


  ! Some hard dimensions:

  INTEGER(I4B), PARAMETER :: ELEM_TET_MATRIX_SIZE = 30 
                                                ! Maximum dimension of elemental
						! tetrahedral matrix. 
						! Presently dimensioned for up to 
						! complete order 2nd order elements
						! ie QT/QN.

  INTEGER(I4B), PARAMETER :: ELEM_TRI_MATRIX_SIZE = 12
                                                ! Maximum dimension of elemental
						! triangular matrix. 
						! Presently dimensioned for up to 
						! complete order 2nd order elements
						! ie QT/QN.

  INTEGER(I4B), PARAMETER :: MAX_PRESENT_ORDER = 2
  LOGICAL(LGT), PARAMETER :: MAX_PRESENT_MIXED = .FALSE.
                                                ! Maximum polynomial order presently
						! implemented. (Currently = QT/QN)


  INTEGER(I4B), PARAMETER :: MAX_QUAD_TRI_RULES = 25

                                                ! Maximum size triangular quadrature
												! rule currently implemented.


  INTEGER(I4B), PARAMETER :: FILENAMELENGTH=128 ! Maximum length of filename.
                                                ! Note: will fail on some OS's eg
                                                ! DOS!


  INTEGER(I4B), PARAMETER :: RUN_LABEL_LINE_LENGTH=78 ! Max length of a run label
                                                      ! comment line.
  INTEGER(I4B), PARAMETER :: RUN_LABEL_NUM_LINES   =5   ! Max number of run label 
                                                      ! comment lines. 
  INTEGER(I4B), PARAMETER :: OPTIONLENGTH = 60


CONTAINS

  SUBROUTINE problem_info_defaults
    REAL_EIGEN_ANALYSIS = .FALSE.
    CMPLX_EIGEN_ANALYSIS = .FALSE.
    CBAA_ANALYSIS = .FALSE.
    GW_ANALYSIS = .FALSE.
    TD_ANALYSIS = .FALSE.
    FD_SCAT_ANALYSIS = .FALSE.
    COUPLED_TD_ANALYSIS = .FALSE.
    WRITE_MESH = .FALSE.
    HIERARCHAL_ORDER = 1
    MIXED_ORDER_FLAG = .TRUE.
    ELEMENT_TYPE = 3   ! Webb99 elements
    CUBATURE = .TRUE.
    CURVILINEAR = .FALSE.
    ALL_CURVILINEAR = .FALSE.
    BANDRENUM = .FALSE.   ! Controls renumbering to reduce matrix bandwidth.
    BANDRENUM_STORE = .FALSE.
    SCALE_BY_EDGE_LENGTH = .FALSE.
    GW_COMPUTE_S_PARAMS = .FALSE.
    PEC_ground = .FALSE.  
    USE_PRE_CONDITIONER = .FALSE.
    OUTPUT_ELEMENT_DATA = .FALSE.
    EXECUTE = .TRUE.
    SPARSE = .TRUE.
    SOLVER_TYPE = 0 ! LU
    ON_SCREEN_REPORTING = .TRUE. ! On-screen monitoring of iterative solver.
    COMPUTE_EIGENMODES = .FALSE.
    COAX_WHITNEY = .FALSE.      ! Coax formulation
    ADAPTIVE_ANALYSIS = .FALSE.
    analysis_counter = 0            ! Used for allocating analysis numbers
    NUM_analysis = 0
    MAX_NODES      = 0
    MAX_ELEMENTS   = 0
    MAX_EDGES      = 0
    MAX_FACES      = 0
    MAX_BOUNDARIES = 0
    MAX_PORTS      = 0
  END SUBROUTINE problem_info_defaults

END MODULE problem_info


