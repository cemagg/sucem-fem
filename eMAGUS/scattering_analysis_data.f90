MODULE scattering_analysis_data
  USE nrtype
  USE problem_info
  SAVE 

  LOGICAL(LGT) SCAT_FIELD       ! Flag indicating scattered field computations.
  LOGICAL(LGT) TEST_INTERPOLATE_FIELD ! Flag for special test mode for 
                                ! interpolated field. 
  LOGICAL(LGT) PEC_SPHERE_SCAT  ! Flag indicating PEC sphere scatterer,found by 
                                ! radius search.
  LOGICAL(LGT) PEC_GENERAL_SCAT ! Flag indicating general PEC scatterer,found 
                                ! by label search.
  LOGICAL(LGT) PENETRABLE_SCAT  ! Flag indicating general penetrable scatterer,
                                ! found by label search.
  REAL(SP) sph_radius           ! Radius of target sphere. Relevant only in 
                                ! scattering analysis.
CONTAINS
  
  SUBROUTINE scattering_data_defaults
    IMPLICIT NONE
    SCAT_FIELD = .FALSE.
    TEST_INTERPOLATE_FIELD = .FALSE.
    PEC_SPHERE_SCAT  = .FALSE.
    PEC_GENERAL_SCAT = .FALSE.
    PENETRABLE_SCAT  = .FALSE.
    sph_radius = 0.0
    
  END SUBROUTINE scattering_data_defaults

!***********************************************************************
! SUBROUTINE scattering_data_check: Checks for inconsistant scattering
!                                   info.
!***********************************************************************
  SUBROUTINE scattering_data_check
    IMPLICIT NONE

    IF (PEC_SPHERE_SCAT.AND.(sph_radius <= 0)) THEN
       STOP 'Spherical scatterer defined with zero radius'
    END IF
    
  END SUBROUTINE scattering_data_check
    
END MODULE scattering_analysis_data
