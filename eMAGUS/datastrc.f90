! Last changed 21 Feb 2002 DBD - quad_tables_tet removed (now in module quad_tab).
! 2002-03-21: Added frac_upgrade variable in MODULE adaptive_fem. MMB.


!***********************************************************************
! Define data structures and interfaces
!***********************************************************************


MODULE adaptive_fem
  USE nrtype
  SAVE

  ! Storage for element volume and face residuals:
  REAL(SP), DIMENSION(:), ALLOCATABLE :: residuals_elements
  REAL(SP), DIMENSION(:), ALLOCATABLE :: residuals_faces

  ! Storage for elemental total error endicators/estimators:
  REAL(SP), DIMENSION(:), ALLOCATABLE :: error_indicators

  ! Storage for indices to increasing order of elemental error 
  ! indicators/estimators:
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: sorted_order

  ! Storage for element and face maximum diameters:
  REAL(SP), DIMENSION(:), ALLOCATABLE :: diameters_elements
  REAL(SP), DIMENSION(:), ALLOCATABLE :: diameters_faces

  ! Storage for error level flags:
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: error_level_flags

  ! A data structure to store the requested adaptive actions:
  TYPE AD_card_record
    INTEGER(I4B) :: file_orders   ! 0/1/2 = none/write/read the element orders to/from a file
    INTEGER(I4B) :: analysis_type ! 0/1 = none/explicit residual error analysis
    INTEGER(I4B) :: upgrade_type  ! 0/1/2/3/4 = none(CT/LN) / LT/LN / LT/QN / QT/QN / evenly graded
    REAL(SP)     :: relative_res  ! relative weight of element volume to face contributions to the residual estimator
    REAL(SP)     :: frac_upgrade  ! Fraction of elements to upgrade to higher order
  END TYPE AD_card_record
  TYPE(AD_card_record) :: ADdata

END MODULE adaptive_fem
!***********************************************************************

MODULE bandwidth
  USE nrtype
  SAVE
  INTEGER(I4B) kl, ku         ! number of super and subdiagonals 
END MODULE bandwidth  
!***********************************************************************

MODULE boundary_conditions
  USE nrtype 
  SAVE

  INTEGER(I4B) num_ports                              ! Number of ports         

  ! Port specification:
  TYPE PORT
    INTEGER(I4B) BC_number                            ! Index to relevant BC.
    INTEGER(I4B) type                                 !  = 1 TE mode
                                                      !  <>1 Unimplemented.

    COMPLEX(SPC) excitation                           ! Relative excitation 
                                                      ! (Re,Im)
    REAL(SP) a                                        ! a and 
    REAL(SP) b                                        ! b dimensions of wg.
    REAL(SP) :: x0,y0                                 ! origin of the port geometry
    REAL(SP) z_coord                                  ! z coordinate of port
    REAL(SP), DIMENSION(3) :: normal                  ! unit normal pointing away from the mesh
    REAL(SP), DIMENSION(3) :: tangent                 ! unit tangential vector
                                                      ! defining positive modal direction .
    REAL(SP), DIMENSION(4,4) :: T                     ! Geometric transformation
                                                      ! matrix for the port. 
    REAL(SP), DIMENSION(4,4) :: T_inv                 ! Geometric inverse 
                                                      ! transformation
                                                      ! matrix for the port. 
  END TYPE PORT
  TYPE(PORT), DIMENSION(:), ALLOCATABLE :: ports


! Added DBD 26 March 2003
  INTEGER(I4B) num_ABCs                               ! Number of ABCs
  INTEGER(I4B) AB_counter

  TYPE ABC
    INTEGER(I4B) label                                ! Unique ABC identifier.
    INTEGER(I4B) BC_num                               ! The associated BC card. 
    INTEGER(I4B) type                                 !  = 1 1st order 
                                                      !  <>1 Unimplemented.
    REAL(SP), DIMENSION(3) :: normal                  ! unit normal pointing away from the mesh (outwards)
	REAL(SP) Yc                                       ! Characteristic admittance of unbounded medium
	REAL(SP) r                                        ! Radius of spherical ABC (used by 2nd order ABC).
!    REAL(SP), DIMENSION(3) :: tangent                 ! unit tangential vector
!                                                      ! defining positive modal direction .
!    REAL(SP), DIMENSION(4,4) :: T                    ! Geometric transformation
!                                                     ! matrix for the ABC. 
!    REAL(SP), DIMENSION(4,4) :: T_inv                ! Geometric inverse 
!                                                     ! transformation
!                                                     ! matrix for the ABC. 
  END TYPE ABC
  TYPE(ABC), DIMENSION(:), ALLOCATABLE :: ABCs

  TYPE CubicBox                                       ! This is a cubix box with sides parallel to the 
    CHARACTER(5) S1,S2                                ! coordinate axes - 1st type of box on QU card. 
  END TYPE CubicBox                                   ! S1 and S2 are the names of points making opposite corners. 
  TYPE(CubicBox) :: ABC_Box

  INTEGER(I4B) num_DPoints
  INTEGER(I4B) DP_counter
  TYPE DPoint
    CHARACTER(5) name                                 ! Name of point
    REAL(SP), DIMENSION(3)::  coords                  ! Coordinates of point x,y,z
  END TYPE DPoint
  TYPE(DPoint), DIMENSION(:), ALLOCATABLE :: DPoints

! End added DBD 26 March 2003

! Added DBD 20 May 2003
!  TYPE Sphere                                         ! This is a sphere as defined in the KU card.
!    CHARACTER(5) S1,S2,S3
!	REAL (SP), DIMENSION(3) :: centre                 ! x,y,z coordinates of centre of sphere
!	REAL (SP) radius                                  ! radius of sphere
!  END TYPE Sphere                                    
!  TYPE(Sphere) :: Spherical_Boundary
! End DBD addition 20 May 2003


  INTEGER(I4B) :: A9_counter,NUM_A9_cards 

  TYPE A9_card_record
    INTEGER(I4B) :: anfl                    ! New/additional excitation
    INTEGER(I4B) :: port_num                ! Associated port
    INTEGER(I4B) :: port_type               ! Type of port excitation
    REAL(SP) :: re_excitation               ! RE(relative excitation)
    REAL(SP) :: im_excitation               ! IM(relative excitation)
! Start DBD additions 20 Mar 2001
    REAL(SP), DIMENSION(3) :: normal        ! Outward directed normal vector
    REAL(SP), DIMENSION(3) :: tangent       ! Tangential vector (see port)
! End DBD additions 20 Mar 2001
    INTEGER(I4B) :: analysis_no             ! analysis number that this card belongs to
  END TYPE A9_card_record
  TYPE(A9_card_record), DIMENSION(:), ALLOCATABLE :: A9data


END MODULE boundary_conditions

!***********************************************************************

MODULE CBAA_data
  USE nrtype 
  SAVE                                 
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ap_elnumbers        ! aperture element numbers
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: which_local_edges ! aperture element 
                                                                 ! local egdes
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: which_local_face    ! aperture element 
                                                                 ! local face num.
  INTEGER(I4B) :: num_apelements                                 ! number of aperture elements
  INTEGER(I4B) :: ap_dof                                         ! number of aperture dof's
  LOGICAL(LGT) :: CBAA_FMM_storage
  LOGICAL(LGT) :: CBAA_FMM_debug

  ! FMM parameters:
  REAL(SP), DIMENSION(4) :: ap_dims     ! extreme co-ordinates of aperture [min(x,y) max(x,y)]
  REAL(SP) :: gr_dim,D_max              ! FMM group dimensions
  INTEGER(I4B) :: gr_numx               ! Number of global FMM groups in x-direction
  INTEGER(I4B) :: gr_num                ! Number of occupied FMM groups
  INTEGER(I4B) :: max_groups            ! Max. number of groups possible
  REAL(SP) :: k0_fmm                    ! used for FMM allocation and interaction determination
  INTEGER(I4B) :: L_max                 ! max. order added to FMM translation term
  REAL(SP) :: FarInteract_min           ! min. FMM far interaction distance


  ! FMM and sparse related storage:
  INTEGER(I4B) :: num_k_dirs                                  ! Number of dirictions: spherical integration
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: k_dirs             ! unit, direction vectors
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: CBAA_BE_mat    ! Full BE matrix
  COMPLEX(SPC), ALLOCATABLE, DIMENSION(:,:) :: temp_BE_mat    ! FMM debug
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: CBAA_BE_val      ! sparse Z'
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: CBAA_BE_rowind   ! sparse Z'
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: CBAA_BE_colind   ! sparse Z'
  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: CBAA_groupmat       ! T in array(k) of VTVt
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: CBAA_elemvec_term1    ! V in array(k) of VTVt
  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: CBAA_elemvec_term2  ! V in array(k) of VTVt
                                                                     ! (terms in Jin, eq9.136)
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: edge_grnum       ! Global group numbers of 
                                                              ! aperture dof's
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: gr_matindex      ! Relates global group number
                                                              ! to an index in [T]
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: gredge_rowind    ! Sparse admin for [Vt]
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: gredge_colind    ! Sparse admin for [Vt]

END MODULE CBAA_data

!***********************************************************************

MODULE coax_feed
  ! Variables describing coax feed to the CBAA, FEM volume:
  ! Last update: 8 June 2000 - MMB : Prettied up
  ! 19 Feb 2001: MMB : created AX card for multiple coax excitations specified within the *.fek file.
  USE nrtype
  SAVE

  INTEGER(I4B) :: AX_counter,NUM_AX_cards 

  TYPE AX_card_record
    INTEGER(I4B) :: anfl                     ! new/additional source = 0/1
    INTEGER(I4B) :: coaxdir                  ! coax direction (x/y/z:1/2/3)
    REAL(SP), DIMENSION(3) :: coaxcentre     ! xyz of aperture centre 
    REAL(SP), DIMENSION(3) :: coaxend        ! xyz of end of extended inner conductor
    REAL(SP) :: coax_a                       ! inner conductor radius
    REAL(SP) :: coax_b                       ! outer conductor radius
    REAL(SP) :: coax_Iabs                    ! magnitude at aperture of ingoing coax current wave
    REAL(SP) :: coax_Iphase                  ! phase at aperture of ingoing coax current wave
    REAL(SP) :: coaxlen                      ! length of inner conductor into the FEM volume
    REAL(SP) :: coax_eps                     ! eps_r of coax filling
    REAL(SP) :: coax_mu                      ! mu_r of coax filling
    INTEGER(I4B), DIMENSION(12) :: coaxedges ! free edges in the coax aperture (6 edges)
    INTEGER(I4B) :: coaxedges_num            ! number of free edges in aperture (6)
    COMPLEX(SPC) :: coax_zin1,coax_zin2      ! input impedance
    REAL(SP), DIMENSION(3) :: normal         ! normal, unit vector away from FEM volume
    REAL(SP) :: Z_c                          ! Characteristic impedance
    COMPLEX(SPC) :: V_plus                   ! Amplitude of volage wave into FEM volume (incident)
    COMPLEX(SPC) :: V_minus                  ! Amplitude of volage wave out of FEM volume (reflected)
    COMPLEX(SPC) :: V_tot                    ! Amplitude of total volage
  END TYPE AX_card_record

  TYPE(AX_card_record), DIMENSION(:), ALLOCATABLE :: AXdata

END MODULE coax_feed

!***********************************************************************

MODULE far_field_data
  USE nrtype
  SAVE
  INTEGER(I4B) :: FF_counter,NUM_FF_cards 

  TYPE FF_card_record
    REAL(SP) :: Theta0
    REAL(SP) :: Phi0
    REAL(SP) :: DTheta
    REAL(SP) :: DPhi
    INTEGER(I4B) :: FFreq
    INTEGER(I4B) :: NTheta
    INTEGER(I4B) :: NPhi
    INTEGER(I4B) :: rige
    INTEGER(I4B) :: analysis_no     ! analysis number that this card belongs to
  END TYPE FF_card_record

  TYPE(FF_card_record), DIMENSION(:), ALLOCATABLE :: FFdata

END MODULE far_field_data

!***********************************************************************

! Start DBD additions 3 Apr 2001.
MODULE gw_data
  USE nrtype
  SAVE
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: S_params
  INTEGER(I4B) num_wg_modes_m, num_wg_modes_n 
                              ! Maximum waveguide indices.
  INTEGER(I4B) mode_m,mode_n  ! Parameters of current analysis 
END MODULE gw_data
!***********************************************************************
! End DBD additions 3 Apr  2001.

MODULE inc_field
  ! See FEKO A0 card for full definitions and explanatory figures.
  USE nrtype
  SAVE

  INTEGER(I4B) :: A0_counter,NUM_A0_cards 

  TYPE A0_card_record
    INTEGER(I4B) anfl               ! New/additional excitation.
    INTEGER(I4B) NThei              ! Number of theta angles.
    INTEGER(I4B) NPhii              ! Number of phi angles.
    REAL(SP) EiR1                   ! Magnitude of Einc (V/m).
    REAL(SP) EiR2                   ! Phase of Einc  (degrees)
    REAL(SP) EiR3                   ! Theta angle of incidence (degrees).
    REAL(SP) EiR4                   ! Phi angle of incidence (degrees).
    REAL(SP) EiR5                   ! Eta angle of polarization (degrees).
    REAL(SP) DThei                  ! Delta theta (degrees).
    REAL(SP) DPhii                  ! Delta phi (degrees).
    INTEGER(I4B) :: analysis_no     ! analysis number that this card belongs to
  END TYPE A0_card_record

  TYPE(A0_card_record), DIMENSION(:), ALLOCATABLE :: A0data

END MODULE inc_field

!***********************************************************************

MODULE matrix
  USE nrtype
  SAVE
  INTEGER(I4B) dof ! degrees of freedom
! DBD addition 6 August 2004
  INTEGER(I4B) pre_dof ! prescribed degrees of freedom
! DBD addition 6 August 2004
  ! The routines that assemble the system matrix, allocate 
  ! real/complex storage respectively.
  REAL(SP),     DIMENSION(:,:), ALLOCATABLE :: s_mat, t_mat
! DBD addition 25 March 2003
  REAL(SP),     DIMENSION(:,:), ALLOCATABLE :: A_mat,B_mat,C_mat
! End DBD addition 25 March 2003
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: s_mat_c, t_mat_c
  REAL(SP),     DIMENSION(:,:), ALLOCATABLE :: sb_mat,tb_mat
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: x_vec_c
! DBD addition 07 August 2004
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: x_pre_vec_c ! Prescribed dof's.
! End DBD addition 07 August 2004
! DBD addition 26 March 2003 and 25-29 April 2003 and 16 July 2003
!  REAL(SP),     DIMENSION(:),   ALLOCATABLE :: x_vec
  REAL(SP),     DIMENSION(:),   ALLOCATABLE :: coupled_td_x_vec
  REAL(SP),     DIMENSION(:),   ALLOCATABLE :: u_nplus1,u_n,u_nmin1,&
                                               f_nplus1,f_n,f_nmin1,b_vec,&
                                               temp_vec_real,&
                                               psi_x_n,psi_y_n,psi_z_n,&
                                               psi_x_nplus1,psi_y_nplus1,psi_z_nplus1

  REAL(DP),     DIMENSION(:),   ALLOCATABLE :: b_vec_DP,x_vec_DP

! End DBD addition 26 March 2003 and 25-29 April 2003 and 16 July 2003
  ! The following are allocated in EIGEN_SYSMAT, either real 
  ! or complex, depending on analysis requested.
  REAL(SP),     DIMENSION(:),   ALLOCATABLE :: eigenvalues
  REAL(SP),     DIMENSION(:,:), ALLOCATABLE :: eigenvectors
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: iwrk
  REAL(SP),     DIMENSION(:),   ALLOCATABLE :: work
  REAL(DP),     DIMENSION(:),   ALLOCATABLE :: rwrk
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: eigenvalues_c
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: eigenvectors_c
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: work_c
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: ipiv
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: A_mat_c
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: full_x_vec_c
  REAL(SP) :: residual_norm      ! Provided for possible interative solvers.
  REAL(SP) :: max_iter_factor    ! Ditto.
  INTEGER(I4B) :: restart_GMRES  ! Ditto.

  ! The full RHS vector including all prescribed edges.
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: b_vec_c
  ! Following are exclusively for sparse matrix solvers:
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: diagonal, lower_tri
  COMPLEX(SPC), DIMENSION(:),   ALLOCATABLE :: workvec1,workvec2

  ! degree of freedom indices (moved from geometry MMB 2001 Sept 25)
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: dof_type      ! Type of hierarchal dof
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_e1 ! d.o.f. "pointers"
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_e2 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_e3 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_f1 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_f2 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_f3 !      ditto.

  ! Following added by DBD 6 Aug 04; 
  ! prescribed "degree of freedom" indices 
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: pre_dof_type      ! Type of hierarchal dof
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_e1 ! d.o.f. "pointers"
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_e2 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_e3 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_f1 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_f2 !      ditto.
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: renumbered_pre_f3 !      ditto.


 ! This is exclusively for ARPACK solver
  INTEGER(I4B) :: mode               ! ARPACK mode 2,3,4 or 5
  INTEGER(I4B) :: maxitr             ! max. iterations, suggested value is 300
  LOGICAL(LGT) :: rvec               ! flag to compute Ritz vectors  ! Corrected DBD 13 Aug 04
  REAL(SP)     :: sigma, sigma2      ! shift for modes 3-5
  CHARACTER    which*2               ! set to 'LM', 'BE' etc. to indicate 
                                     ! which part of eigenvalue spectrum is desired
  CHARACTER       bmat
  INTEGER(I4B)  :: nev               ! number of eigenvalues requested
  INTEGER(I4B)  :: ncv               ! Represents the dimension of the Lanczos basis 
                                     ! constructed by ssaupd for OP.
  INTEGER(I4B) :: nconv_eigenvalues  ! the number of converged eigenvalues (only
                                     ! relevant when ARPACK was used)
  REAL(SP),     DIMENSION(:), ALLOCATABLE :: S_min_sigT ! added by JPS 22 Jul 05. Sparse matrix S-sigma*T
                                                        ! Required for ARPACK 
                                                        ! shift invert mode linear system solve.
 

 !This is exclusively for sparse solver
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: Asparse_c, Ssparse_c, Tsparse_c
! Changed DBD 24 April 2003
  REAL(SP),     DIMENSION(:), ALLOCATABLE :: Asparse, Bsparse, Csparse, Ssparse, Tsparse, LUval
  
  REAL(DP),     DIMENSION(:), ALLOCATABLE :: Asparse_skyline 
  REAL(DP),     DIMENSION(:), ALLOCATABLE :: Asparse_DP 
											 
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: IAUdiag
  INTEGER(I4B) :: num_nonzeros_sky

! End changed DBD 24 April 2003
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: col_ind, row_ind
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ucolind, urowind
  REAL(SP),     DIMENSION(:), ALLOCATABLE :: Upperval, Lowerval, CRSval
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: upper_rowind, upper_colind
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: CRSupper_row, CRSupper_col

CONTAINS
!*******************************************************************************
!* SUBROUTINE matrix_defaults
!* 
!* Initialises _some_ of the matrix module data members
!*
!*******************************************************************************
  SUBROUTINE matrix_defaults
    IMPLICIT NONE
    
    ! Defaults for matrix storage and linear algebra
    residual_norm = 1.0E-4_SP ! Default residual norm for linear algebra
    max_iter_factor = 2.0_SP  ! Default factor (times dof) specifying maximum iterations.
    restart_GMRES = 0
  END SUBROUTINE matrix_defaults


END MODULE matrix
!***********************************************************************

MODULE near_field_data
  USE nrtype
  SAVE
  INTEGER(I4B) :: FE_counter,NUM_FE_cards 
  LOGICAL(LGT) FIELD_AVERAGING  
  REAL(SP) :: eps_field_averaging

  TYPE FE_card_record
    REAL(SP) :: delta_x
    REAL(SP) :: delta_y
    REAL(SP) :: delta_z
    REAL(SP) :: x0
    REAL(SP) :: y0
    REAL(SP) :: z0
    INTEGER(I4B) :: felkor
    INTEGER(I4B) :: feltyp
    INTEGER(I4B) :: n_x
    INTEGER(I4B) :: n_y
    INTEGER(I4B) :: n_z
    INTEGER(I4B) :: analysis_no     ! analysis number that this card belongs to
  END TYPE FE_card_record
  TYPE(FE_card_record), DIMENSION(:), ALLOCATABLE :: FEdata

END MODULE near_field_data
!***********************************************************************

MODULE probe_feed
! Variables for CBAA current probe exitation card (A8 card):
! Last update: 8 June 2000 - MMB : Created

  USE nrtype
  SAVE
  INTEGER(I4B) :: A8_counter,NUM_A8_cards 
  LOGICAL(LGT) :: probe_present

  TYPE A8_card_record
    INTEGER(I4B) :: anfl          ! new/additional source = 0/1
    REAL(SP) :: prabs
    REAL(SP) :: prphase           ! current magnitude and phase
    REAL(SP) :: prlen             ! length of probe
    REAL(SP) :: prx0
    REAL(SP) :: pry0
    REAL(SP) :: prz0              ! xyz of starting point (current defined as
                                  ! flowing away from this point)
    REAL(SP) :: prrad             ! radius of probe
    INTEGER(I4B) :: prdir         ! x/y/z : 1/2/3 direction of probe
    COMPLEX(SPC) :: prvolt        ! resulting voltage as calculated
    INTEGER(I4B) :: analysis_no   ! analysis number that this card belongs to
  END TYPE A8_card_record

  TYPE(A8_card_record), DIMENSION(:), ALLOCATABLE :: A8data

END MODULE probe_feed

!***********************************************************************

MODULE renumbering
  USE nrtype
  SAVE
  REAL(SP), DIMENSION(3) :: length_array, lengths, min_array, max_array 
  INTEGER(I4B) :: maxpos, minpos, midpos, debug_var_int    
  LOGICAL(LGT) :: ascending
END MODULE renumbering

!***********************************************************************

MODULE unit_numbers
  USE nrtype
  INTEGER(I4B),PARAMETER ::INFILE=9,MESHIN=10,FILEOUT=11,FNFOUT=12, & 
                           TESTFILE = 13,FEKO_OUT=14, AUXINFILE=15, &
                           MESHOUT = 16
  ! Input, mesh and output file numbers

END MODULE unit_numbers
!***********************************************************************


