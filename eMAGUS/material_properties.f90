MODULE material_properties
  USE nrtype

  ! Material related data structures:
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: eps_r    ! Relative permittivity
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: eps_r_xx ! Diagonally 
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: eps_r_yy ! anisotropic 
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: eps_r_zz ! relative permittivity
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: mu_r     ! Relative permeability
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: mu_r_xx  ! Diagonally 
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: mu_r_yy  ! anisotropic 
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: mu_r_zz  ! relative permeability
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: material_defined 
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: material_type 
  ! Maximum number of different materials including internal PML materials.
  INTEGER(i4b) :: max_materials 

  INTEGER(I4B) :: DI_counter,NUM_DI_cards 
  TYPE DI_card_record
    INTEGER(I4B) :: medium                          ! Medium number 
    INTEGER(I4B) :: material_type                   ! Type of material
    REAL(SP) :: eps_r_prime                         ! Real relative permittivity
    REAL(SP) :: mu_r_prime                          ! Real relative permeability
    REAL(SP) :: sigma                               ! Electric loss 
    REAL(SP) :: tandmue                             ! Magnetic loss 
    REAL(SP) :: tandeps                             ! Electric loss (alternate specification)
    INTEGER(I4B) :: analysis_no                     ! analysis number that this card belongs to
  END TYPE DI_card_record

!! DIdata never deallocated? NM

  TYPE(DI_card_record), DIMENSION(:), ALLOCATABLE :: DIdata
  
!! Obsolete, but retained for possible later use.
!! Well, actually it's still used by meshin_old NM
  INTEGER(I4B) :: MA_counter,NUM_MA_cards 
  TYPE MA_card_record
    INTEGER(I4B) :: label                           ! Material label
    INTEGER(I4B) :: material_type                   ! Type of material
    COMPLEX(SPC) :: eps_r                           ! isotropic eps_r
    COMPLEX(SPC) :: mu_r                            ! isotropic mu_r
    COMPLEX(SPC), DIMENSION(3,3) :: eps_r_tensor    ! anisotropic eps_r
    COMPLEX(SPC), DIMENSION(3,3) :: mu_r_tensor     ! anisotropic mu_r
    INTEGER(I4B) :: analysis_no                     ! analysis number that this 
                                                    ! card belongs to
  END TYPE MA_card_record
  TYPE(MA_card_record), DIMENSION(:), ALLOCATABLE :: MAdata
! End obsolete code

CONTAINS
  
  SUBROUTINE material_properties_allocate(no_materials)
!********************************************************************
!* Allocate memory for no_materials number of material types, and 
!* sets the module variable max_materials
!********************************************************************
    INTEGER(I4B), INTENT(in) :: no_materials
    max_materials = no_materials
    ALLOCATE (eps_r(0:max_materials))
    ALLOCATE (eps_r_xx(0:max_materials))
    ALLOCATE (eps_r_yy(0:max_materials))
    ALLOCATE (eps_r_zz(0:max_materials))
    ALLOCATE (mu_r(0:max_materials))
    ALLOCATE (mu_r_xx(0:max_materials))
    ALLOCATE (mu_r_yy(0:max_materials))
    ALLOCATE (mu_r_zz(0:max_materials))
    ALLOCATE (material_type(0:max_materials+7)) ! To allow space for PML's 
                                                ! properties. 
    ALLOCATE (material_defined(0:max_materials))
  END SUBROUTINE material_properties_allocate

!!!
!!! SUBROUTINE material_properties_initialise:
!!!
!!! Initialises the materials once they have been allocated
!!!
  SUBROUTINE material_properties_init
    IMPLICIT NONE
    ! First initialize ALL materials as free space, isotropic and undefined.
    ! Note that PML materials are handled differently.
    eps_r(0:MAX_MATERIALS)    = (1.0_SP,0.0_SP)
    mu_r(0:MAX_MATERIALS)     = (1.0_SP,0.0_SP)
    eps_r_xx(0:MAX_MATERIALS) = (1.0_SP,0.0_SP)
    eps_r_yy(0:MAX_MATERIALS) = (1.0_SP,0.0_SP)
    eps_r_zz(0:MAX_MATERIALS) = (1.0_SP,0.0_SP)
    mu_r_xx(0:MAX_MATERIALS)  = (1.0_SP,0.0_SP)
    mu_r_yy(0:MAX_MATERIALS)  = (1.0_SP,0.0_SP)
    mu_r_zz(0:MAX_MATERIALS)  = (1.0_SP,0.0_SP)
    material_type(0:MAX_MATERIALS) = 1 ! isotropic
    material_defined(0:MAX_MATERIALS) = .FALSE.

  END SUBROUTINE material_properties_init
  
  
  SUBROUTINE material_properties_deallocate
    DEALLOCATE (eps_r)
    DEALLOCATE (eps_r_xx)
    DEALLOCATE (eps_r_yy)
    DEALLOCATE (eps_r_zz)
    DEALLOCATE (mu_r)
    DEALLOCATE (mu_r_xx)
    DEALLOCATE (mu_r_yy)
    DEALLOCATE (mu_r_zz)
    DEALLOCATE (material_type)
    IF(ALLOCATED(MAdata))  DEALLOCATE (MAdata) ! MA card parameters
    IF(ALLOCATED(material_defined))  DEALLOCATE (material_defined) ! 

  END SUBROUTINE material_properties_deallocate

  SUBROUTINE material_properties_errcheck
    USE problem_info
    USE output_error

    INTEGER(i4b) :: i_mat

    ! error checking:
    IF (REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS) THEN 
      DO i_mat = 0, MAX_MATERIALS
        IF ( (AIMAG(eps_r(i_mat)).GT.EPS).OR.(AIMAG(mu_r(i_mat)).GT.EPS) ) THEN
          IF (REAL_EIGEN_ANALYSIS) THEN
            CALL ERROR_FEMFEKO(1,4042,int1=i_mat)
          END IF
        END IF
        IF (material_type(i_mat).NE.1) THEN 
          CALL ERROR_FEMFEKO(1,4043,int1=i_mat,int2=material_type(i_mat))
        END IF
      END DO
    END IF
  END SUBROUTINE material_properties_errcheck


END MODULE material_properties
