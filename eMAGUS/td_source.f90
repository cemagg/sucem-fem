! Last changed 21 Feb 2002, DBD. Some typos corrected. 
! Tetrahedral rules added. 

MODULE TD_source
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This module contains datastractures and routines
! related to the time domain source and time domain 
! code operation. 
!
! Because some of the routines calling this are re-cycled versions of 
! guided wave analysis routines, and it is not desirable to change 
! their arguments, a number of parameters for the "high-level" function 
! TD_E_INC are passed via the module and not
! via the argument list (eg angles etc). 
! Created: 24 Mar 2003, DBD.
!*******************************************************************************
   SAVE

   INTEGER (I4B) :: timestep                    
   REAL(SP) :: delta_t,Newmark_beta
   REAL(SP) :: sigma_pulse,offset_pulse
   INTEGER(I4B) :: num_timesteps
   INTEGER(I4B) :: TD_source_type
   REAL(SP) :: t_pw,p_pw,e_pw                             ! Polarization properties of the incident plane wave     
                                                          ! theta,phi,eta (all radians)  
   REAL(SP) :: mag_pw                                     ! Magnitude of incident plane wave     
  
   REAL(SP), DIMENSION(3) :: O1,O2,O3,O4,O1p,O2p,O3p,O4p  ! See Taflove and Hagness, 
                                                          ! "Computational Electrodynamics: the FDTD", 2nd edition, 
                                                          ! Artech House 2000,  Fig. 5.14. "p" is "prime".
                                                          

CONTAINS 

FUNCTION TD_INC_FIELD(x,y,z,field_type,Y_char)
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the field a point (x,y,z) at timestep t for the 
! the plane wave characterized by propagation angles theta (t_pw) and phi (phi_pw),
! polarization angle eta (t_pw) and magnitude mag_pw.
!
! The routine returns either the electric field, or the accompanying magnetic field WITHOUT the correct
! Yc scaling. (This is done to simplify subseqeunt operations). 
! 
! An option has now been added to include scaling by optional parameter Y_char for the H field case. 
!*******************************************************************************
  REAL(SP), INTENT(IN) :: x,y,z                        ! Coordinates of field point. 
  REAL(SP), DIMENSION(3) :: TD_INC_FIELD               ! x, y and z components of field at point and time. 
  CHARACTER, INTENT(IN) :: field_type                  ! Electric or magnetic (without Yc term) field

  REAL(SP) delay, time, timefunction 
  REAL(SP), INTENT(IN), OPTIONAL :: Y_char             ! Characteristic admittance.

  delay = TD_INC_FIELD_DELAY(x,y,z) 
  time = timestep*delta_t - delay
  SELECT CASE (TD_source_type) 
  CASE (1)
   timefunction  = GAUSSDER(time)
  CASE DEFAULT
	STOP 'Unimplemented time domain source'
  END SELECT 

  SELECT CASE (field_type) 
  CASE ('E','e')
    TD_INC_FIELD = timefunction * mag_pw * TD_E_INC_POLARIZATION() 
  CASE ('H','h')
    TD_INC_FIELD = timefunction * mag_pw * TD_H_INC_POLARIZATION() 
	IF (PRESENT(Y_char)) THEN
      TD_INC_FIELD  = TD_INC_FIELD * Y_char
	END IF
  CASE DEFAULT
	STOP 'Unimplemented field type'
  END SELECT 
END FUNCTION TD_INC_FIELD


FUNCTION TD_D_BY_DT_INC_FIELD(x,y,z,field_type,Y_char)
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the time differential of the field at point (x,y,z) at timestep t for the 
! the plane wave characterized by propagation angles theta (t_pw) and phi (phi_pw),
! polarization angle eta (t_pw) and magnitude mag_pw.
!
! The routine returns either the electric field, or the accompanying magnetic field WITHOUT the correct
! Yc scaling. (This is done to simplify subseqeunt operations) but with the correct (minus) sign.
! 
! An option has now been added to include scaling by optional parameter Y_char for the H field case. 
!*******************************************************************************
!  REAL(SP), INTENT(IN) :: timestep                    
  REAL(SP), INTENT(IN) :: x,y,z                        ! Coordinates of field point. 
!  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pw              ! Polarization properties of the incident plane wave     
!                                                      ! theta,phi,eta (all radians)  
!  REAL(SP), INTENT(IN) :: mag_pw                      ! Magnitude of incident plane wave     
  REAL(SP), DIMENSION(3) :: TD_D_BY_DT_INC_FIELD       ! x, y and z components of field at point and time. 
  CHARACTER, INTENT(IN) :: field_type                  ! Electric or magnetic (without Yc term) field
  REAL(SP), INTENT(IN), OPTIONAL :: Y_char             ! Characteristic admittance.

  REAL(SP) delay, time, d_by_dt_timefunction 

  delay = TD_INC_FIELD_DELAY(x,y,z) 
  time = timestep*delta_t - delay
  SELECT CASE (TD_source_type) 
  CASE (1)
   d_by_dt_timefunction  =   D_BY_DT_GAUSSDER(time)
  CASE DEFAULT
	STOP 'Unimplemented time domain source'
  END SELECT 

  SELECT CASE (field_type) 
  CASE ('E','e')
    TD_D_BY_DT_INC_FIELD = d_by_dt_timefunction * mag_pw * TD_E_INC_POLARIZATION() 
  CASE ('H','h')
    TD_D_BY_DT_INC_FIELD = - d_by_dt_timefunction * mag_pw * TD_H_INC_POLARIZATION() 
	IF (PRESENT(Y_char)) THEN
      TD_D_BY_DT_INC_FIELD  = TD_D_BY_DT_INC_FIELD * Y_char
	END IF
  CASE DEFAULT
	STOP 'Unimplemented field type'
  END SELECT 
END FUNCTION TD_D_BY_DT_INC_FIELD


!******************************************************************************


FUNCTION TD_D2_BY_DT2_INC_FIELD(x,y,z,field_type,Y_char)
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the second time differential of the field at point (x,y,z) at timestep t for the 
! the plane wave characterized by propagation angles theta (t_pw) and phi (phi_pw),
! polarization angle eta (t_pw) and magnitude mag_pw.
!
! The routine returns either the electric field, or the accompanying magnetic field WITHOUT the correct
! Yc scaling. (This is done to simplify subseqeunt operations).
! 
! An option has now been added to include scaling by optional parameter Y_char for the H field case.
!*******************************************************************************
!  REAL(SP), INTENT(IN) :: timestep                    
  REAL(SP), INTENT(IN) :: x,y,z                        ! Coordinates of field point. 
!  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pw              ! Polarization properties of the incident plane wave     
!                                                      ! theta,phi,eta (all radians)  
!  REAL(SP), INTENT(IN) :: mag_pw                      ! Magnitude of incident plane wave     
  REAL(SP), DIMENSION(3) :: TD_D2_BY_DT2_INC_FIELD       ! x, y and z components of field at point and time. 
  CHARACTER, INTENT(IN) :: field_type                  ! Electric or magnetic (without Yc term) field
  REAL(SP), INTENT(IN), OPTIONAL :: Y_char             ! Characteristic admittance.

  REAL(SP) delay, time, d2_by_dt2_timefunction

  delay = TD_INC_FIELD_DELAY(x,y,z) 
  time = timestep*delta_t - delay
  SELECT CASE (TD_source_type) 
  CASE (1)
   d2_by_dt2_timefunction  =   D2_BY_DT2_GAUSSDER(time)
  CASE DEFAULT
	STOP 'Unimplemented time domain source'
  END SELECT 

  SELECT CASE (field_type) 
  CASE ('E','e')
    TD_D2_BY_DT2_INC_FIELD = d2_by_dt2_timefunction * mag_pw * TD_E_INC_POLARIZATION() 
  CASE ('H','h')
    TD_D2_BY_DT2_INC_FIELD = d2_by_dt2_timefunction * mag_pw * TD_H_INC_POLARIZATION() 
	IF (PRESENT(Y_char)) THEN
      TD_D2_BY_DT2_INC_FIELD  = TD_D2_BY_DT2_INC_FIELD * Y_char
	END IF
  CASE DEFAULT
	STOP 'Unimplemented field type'
  END SELECT 
END FUNCTION TD_D2_BY_DT2_INC_FIELD

!******************************************************************************


FUNCTION TD_D3_BY_DT3_INC_FIELD(x,y,z,field_type,Y_char)
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the third time differential of the field at point (x,y,z) at timestep t for the 
! the plane wave characterized by propagation angles theta (t_pw) and phi (phi_pw),
! polarization angle eta (t_pw) and magnitude mag_pw.
!
! The routine returns either the electric field, or the accompanying magnetic field WITHOUT the correct
! Yc scaling. (This is done to simplify subseqeunt operations).
! 
! An option has now been added to include scaling by optional parameter Y_char for the H field case.
!*******************************************************************************
!  REAL(SP), INTENT(IN) :: timestep                    
  REAL(SP), INTENT(IN) :: x,y,z                        ! Coordinates of field point. 
!  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pw              ! Polarization properties of the incident plane wave     
!                                                      ! theta,phi,eta (all radians)  
!  REAL(SP), INTENT(IN) :: mag_pw                      ! Magnitude of incident plane wave     
  REAL(SP), DIMENSION(3) :: TD_D3_BY_DT3_INC_FIELD       ! x, y and z components of field at point and time. 
  CHARACTER, INTENT(IN) :: field_type                  ! Electric or magnetic (without Yc term) field
  REAL(SP), INTENT(IN), OPTIONAL :: Y_char             ! Characteristic admittance.

  REAL(SP) delay, time, d3_by_dt3_timefunction

  delay = TD_INC_FIELD_DELAY(x,y,z) 
  time = timestep*delta_t - delay
  SELECT CASE (TD_source_type) 
  CASE (1)
   d3_by_dt3_timefunction  =   D3_BY_DT3_GAUSSDER(time)
  CASE DEFAULT
	STOP 'Unimplemented time domain source'
  END SELECT 

  SELECT CASE (field_type) 
  CASE ('E','e')
    TD_D3_BY_DT3_INC_FIELD = d3_by_dt3_timefunction * mag_pw * TD_E_INC_POLARIZATION() 
  CASE ('H','h')
    TD_D3_BY_DT3_INC_FIELD = d3_by_dt3_timefunction * mag_pw * TD_H_INC_POLARIZATION() 
	IF (PRESENT(Y_char)) THEN
      TD_D3_BY_DT3_INC_FIELD  = TD_D3_BY_DT3_INC_FIELD * Y_char
	END IF
  CASE DEFAULT
	STOP 'Unimplemented field type'
  END SELECT 
END FUNCTION TD_D3_BY_DT3_INC_FIELD

!******************************************************************************

SUBROUTINE TD_INC_FIELD_SETUP 
  USE boundary_conditions
  USE geometry
  USE inc_field
  USE nrtype
  USE output_error
  USE problem_info
  use scattering_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! This routine sets up the points required by the algorithm in 
! Taflove and Hagness, "Computational Electrodynamics: the FDTD", 2nd edition, 
! Artech House 2000 in Section 5.8 is used.

! It also returns the properties of the plane wave from the relevant A9 card.
!  
! It is presently assumed that accompanying this 
! card are six ABC's numbered 1-6 and forming 
! a rectangular box aligned with the x-y-z axes. 
! Some error checking is also done. 

! It also computes the peak value of the time shape specified 
! for normalization purposes. 
!*******************************************************************************
! Original routine: D B Davidson, April 2003.
! Extended:         DBD, 20 May 2003. 
!*******************************************************************************
 
  ! Assign eight points O1 - O4 and O1' - O4'.
  INTEGER(I4B) ii
  LOGICAL (LGT) S1_found,S2_found,S3_found
  REAL(SP), DIMENSION(3) :: peak_field 
  REAL(SP), DIMENSION(3) :: sphereS2,sphereS3
  REAL(SP) temp1,temp2

  S1_found = .false.
  S2_found = .false.


  ! First, assign opposite points corresponding to S1 and S2.
  DO ii = 1,num_DPoints
    IF(ABC_Box%S1.EQ.DPoints(ii)%name) THEN
      O1(1:3) = DPoints(ii)%coords
      S1_found = .true.
	ELSE IF (ABC_Box%S2.EQ.DPoints(ii)%name) THEN
      O3p(1:3) = DPoints(ii)%coords
	  S2_found = .true.
    END IF
  END DO

  ! Check that S1 and S2 were found.
  IF(.NOT.S1_found.OR..NOT.S2_found) THEN
    CALL ERROR_FEMFEKO(1,4091)
  END IF


  ! Check that S1 is closer to origin the S2 by checking distance from origin. 
  IF(DOT_PRODUCT(O1,O1).GE.DOT_PRODUCT(O3p,O3p)) THEN
    CALL ERROR_FEMFEKO(1,4092)
  END IF
  
  ! Now build up other points - first on "bottom" (see Fig. 5.14b, T&H):
  O2 = O1
  O2(1) = O3p(1)
  O3 = O2
  O3(2) = O3p(2)
  O4 = O1
  O4(2) = O3p(2)
  ! and now on "top".
  O1p = O1
  O1p(3) = O3p(3)
  O2p = O2
  O2p(3) = O3p(3)
  O4p = O4
  O4p(3) = O3p(3)
  IF(DEBUG_TD) THEN
    print *,'O1,O2,O3,O4,O1p,O2p,O3p,O4p'
    print *,O1,O2,O3,O4,O1p,O2p,O3p,O4p
  END IF
  
  IF(NUM_A0_cards.NE.1) THEN
    CALL ERROR_FEMFEKO(1,4093)
  END IF
  mag_pw = A0data(1)%EiR1
  IF (ABS(A0data(1)%EiR2).GT.EPS) THEN
    CALL ERROR_FEMFEKO(1,4094)
  END IF
  t_pw   = A0data(1)%EiR3 * D2R
  p_pw   = A0data(1)%EiR4 * D2R
  e_pw   = A0data(1)%EiR5 * D2R


  WRITE(FILEOUT,'(//,16X,(A))') 'EXCITATION BY INCIDENT PLANE ELECTROMAGNETIC WAVE'

  WRITE(FILEOUT,'(/,1X,(A))') 'Number of excitation:   N =      1 (only one permitted)'
  IF(TD_SOURCE_TYPE.EQ.1) THEN 
    WRITE(FILEOUT,'(1X,(A))') 'Differentiated Gaussian'
    WRITE(FILEOUT,'(1X,(A))') 'pulse used'
    WRITE(FILEOUT,'(1X,(A),T26,(A),G12.4)') 'Pulse parameters:','SIGMA =',sigma_pulse
  ELSE 
    STOP 'Unknown pulse type'
  ENDIF
  WRITE(FILEOUT,'(1X,T26,(A),G12.4)')     'PULSE DELAY = ',offset_pulse
  WRITE(FILEOUT,'(1X,(A),T26,(A),F7.2,2X,(A),F7.2)') 'Direction of incidence:','THETA =',t_pw*r2d,'PHI =',p_pw*r2d 
  WRITE(FILEOUT,'(1X,(A),T26,(A))')       'Polarisation:','LINEAR (only)'
  WRITE(FILEOUT,'(1X,(A),T26,(A))')       'Axial ratio:','V   =    0.0000 (only)'
  WRITE(FILEOUT,'(1X,(A),T26,(A),F7.2)')  'Dir. of polarisation:','ETA =',e_pw*r2d
  WRITE(FILEOUT,'(1X,(A),T26,(A),G12.5)') 'Direction of propag.','k^0X =',-SIN(t_pw)*COS(p_pw)
  WRITE(FILEOUT,'(1X,T26,(A),G12.5)')     'k^0Y =',-SIN(t_pw)*SIN(p_pw)
  WRITE(FILEOUT,'(1X,T26,(A),G12.5)')     'k^0Z =',-COS(t_pw)
  peak_field = TD_E_INC_POLARIZATION()
  WRITE(FILEOUT,'(1X,(A),T26,(A),F8.2)')  'Field strength in V/m:','E0X = ',peak_field(1)*mag_pw
  WRITE(FILEOUT,'(1X,(A),T26,(A),F8.2)')  '(Phase meangingless','E0Y = ',peak_field(2)*mag_pw
  WRITE(FILEOUT,'(1X,(A),T26,(A),F8.2)')  'in time domain):','E0Z = ',peak_field(3)*mag_pw


  
  SCAT_FIELD_TEST: IF(SCAT_FIELD) THEN
    WRITE(FILEOUT,'(//,16X,(A),I3)') 'SCATTERED FIELD FORMULATION. HOMOGENEOUS EMBEDDING MATERIAL NUMBER: ',HOMOG_MEDIUM
  ELSE
    WRITE(FILEOUT,'(//,16X,(A))') 'TOTAL FIELD FORMULATION.' 
  END IF SCAT_FIELD_TEST

  !SCAT_FIELD_TEST: IF(SCAT_FIELD) THEN
  !  S1_found = .false.
  !  S2_found = .false.
  !  S3_found = .false.
  !  ! Find parameters of spherical boundary.
  !  DO ii = 1,num_DPoints
  !    IF(Spherical_Boundary%S1.EQ.DPoints(ii)%name) THEN
  !      Spherical_Boundary%centre(1:3) = DPoints(ii)%coords
  !      S1_found = .true.
  ! 	  ELSE IF (Spherical_Boundary%S2.EQ.DPoints(ii)%name) THEN
  !       sphereS2(1:3) = DPoints(ii)%coords
  !       S2_found = .true.
  !	  ELSE IF (Spherical_Boundary%S3.EQ.DPoints(ii)%name) THEN
  !      sphereS3(1:3) = DPoints(ii)%coords
  !      S3_found = .true.
  !  END IF
  ! END DO
  !  ! Check that S1, S2 and S3 were found.
  !  IF(.NOT.S1_found.OR..NOT.S2_found.OR..NOT.S3_found) THEN
  !    CALL ERROR_FEMFEKO(1,4098)
  !  END IF
  !  temp1 = SQRT(DOT_PRODUCT(sphereS2(1:3)-Spherical_Boundary%centre(1:3),sphereS2(1:3)-Spherical_Boundary%centre(1:3)))
  !  temp2 = SQRT(DOT_PRODUCT(sphereS3(1:3)-Spherical_Boundary%centre(1:3),sphereS3(1:3)-Spherical_Boundary%centre(1:3)))
  !  IF (ABS(temp1-temp2).GE.EPS) THEN
  !    CALL ERROR_FEMFEKO(1,4098)
  !  ELSE 
  !	  Spherical_Boundary%radius = temp1
  !  END IF
  
  !  WRITE(FILEOUT,'(//,16X,(A))') 'SCATTERED/TOTAL FIELD INTERFACE: SPHERICAL BOUNDARY ENTIRELY ENCLOSING SCATTERER:'
  !  WRITE(FILEOUT,'(/,1X,(A,3(G12.5,1X)))') 'Centre of spherical boundary:',Spherical_Boundary%centre
  !  WRITE(FILEOUT,'(/,1X,(A,G12.5))')       'Radius of spherical boundary:',Spherical_Boundary%radius
  !END IF SCAT_FIELD_TEST


END SUBROUTINE TD_INC_FIELD_SETUP


!******************************************************************************


FUNCTION TD_INC_FIELD_DELAY(x,y,z) 
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the delay time  [in seconds] at point x,y,z on the FE
! boundary assuming a rectangular box aligned with the x-y-z axes as boundary.
! 
! To compute the relative delays, the algorithm in 
! Taflove and Hagness, "Computational Electrodynamics: the FDTD", 2nd edition, 
! Artech House 2000 in Section 5.8 is used with some very minor modifications.
! In the FDTD case, the phase velocity at an angle is known; this is not so in 
! an irregular FE mesh, so the phase velocity is assumed to be that of free space
!*******************************************************************************
  
  REAL(SP), INTENT(IN) :: x,y,z                          ! Coordinates of field point. 
  REAL(SP) TD_INC_FIELD_DELAY                            ! Delay in seconds relative to first incidence on FE boundary.
  REAL(SP), DIMENSION(3) :: r_comp                       ! Position vector from appropriate origin to field point. 
                                                         ! [T&H,eqn.(5.60),p.221], but in m, not "grid cells"
  REAL(SP), DIMENSION(3) :: k_inc                        ! Unit incident wavevector 
  REAL(SP), DIMENSION(3) :: local_origin                 ! Local origin - see [T&H,p.222]
  REAL(SP), DIMENSION(3) :: field_point                  ! Point (x,y,z)
  REAL(SP) d                                             ! Distance d [T&H,eqn.(5.41),p.221], but in m, not "grid cells"
  INTEGER (I4B) ABC_num,ii

  field_point(1) = x
  field_point(2) = y
  field_point(3) = z


! Locate local origins: (see p.222 Taflove and Hagness)
! Note that because FEKO defines the angles differently, the origins differ.

  IF (0.LE.t_pw.AND.t_pw.LE.PI/2) THEN
    IF (0.LE.p_pw.AND.p_pw.LE.PI/2) THEN
      local_origin  = O3p
	ELSE IF (PI/2.LT.p_pw.AND.p_pw.LE.PI) THEN
      local_origin  = O4p
	ELSE IF (PI.LT.p_pw.AND.p_pw.LE.3.0_SP*PI/2.0_SP) THEN
      local_origin  = O1p
	ELSE IF (3.0_SP*PI/2.0_SP.LT.p_pw.AND.p_pw.LT.2.0_SP*PI) THEN
      local_origin  = O2p
	ELSE
	  STOP 'Invalid incident angle phi: range [0,360)'
	END IF
  ELSE IF (PI/2.LT.t_pw.AND.t_pw.LT.PI) THEN
    IF (0.LE.p_pw.AND.p_pw.LE.PI/2) THEN
      local_origin  = O3
	ELSE IF (PI/2.LT.p_pw.AND.p_pw.LE.PI) THEN
      local_origin  = O4
	ELSE IF (PI.LT.p_pw.AND.p_pw.LE.3.0_SP*PI/2.0_SP) THEN
      local_origin  = O1
	ELSE IF (3.0_SP*PI/2.0_SP.LT.p_pw.AND.p_pw.LT.2.0_SP*PI) THEN
      local_origin  = O2
	END IF
	  STOP 'Invalid incident angle phi: range [0,360)'
  ELSE
	  STOP 'Invalid incident angle theta: range [0,180)'
  END IF 

  ! Following is [T&H,eqn.(5.60),p.221], FEKO convention has unit vector k in OPPOSITE direction to T&H, hence minus signs. 
  k_inc(1) = - sin(t_pw)*cos(p_pw) 
  k_inc(2) = - sin(t_pw)*sin(p_pw) 
  k_inc(3) = - cos(t_pw)

  r_comp(1:3) =  field_point(1:3) - local_origin(1:3) ! Note that this is the distance in m, not in "grid cells"
  d = DOT_PRODUCT(k_inc,r_comp)
  TD_INC_FIELD_DELAY = d/c_0 ! Delay in seconds.

END FUNCTION TD_INC_FIELD_DELAY


!******************************************************************************


FUNCTION TD_E_INC_POLARIZATION()  
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the Cartesian components for the incident 
! time-domain plane wave polarization. 
!*******************************************************************************
!  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pw                 ! Polarization properties of the incident plane wave     
!                                                         ! theta,phi,eta (all radians)  
  
  REAL(SP), DIMENSION(3) :: TD_E_INC_POLARIZATION        ! x, y and z components planewave polarized as spedified. 

  TD_E_INC_POLARIZATION(1) = -COS(e_pw)*COS(t_pw)*COS(p_pw)-SIN(e_pw)*SIN(p_pw)
  TD_E_INC_POLARIZATION(2) = -COS(e_pw)*COS(t_pw)*SIN(p_pw)+SIN(e_pw)*COS(p_pw)
  TD_E_INC_POLARIZATION(3) = COS(e_pw)*SIN(t_pw)
END FUNCTION TD_E_INC_POLARIZATION


FUNCTION TD_H_INC_POLARIZATION()  ! (t_pw,p_pw,e_pw)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the Cartesian components of the associated
! magnetic field for the incident time-domain plane wave polarization. 
!*******************************************************************************
!  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pw                 ! Polarization properties of the incident plane wave     
!                                                         ! theta,phi,eta (all radians)  
  
  REAL(SP) e_pwH
  REAL(SP), DIMENSION(3) :: TD_H_INC_POLARIZATION        ! x, y and z components planewave polarized as spedified. 

  e_pwH = e_pw + PI/2.0_SP
  TD_H_INC_POLARIZATION(1) = -COS(e_pwH)*COS(t_pw)*COS(p_pw)-SIN(e_pwH)*SIN(p_pw)
  TD_H_INC_POLARIZATION(2) = -COS(e_pwH)*COS(t_pw)*SIN(p_pw)+SIN(e_pwH)*COS(p_pw)
  TD_H_INC_POLARIZATION(3) = COS(e_pwH)*SIN(t_pw)
END FUNCTION TD_H_INC_POLARIZATION

!******************************************************************************


FUNCTION GAUSSDER(t)
  USE nrtype
  USE output_error
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the Gaussian derivative pulse at time t, scaled so that 
! its peak value is unity for convenience. 
!*******************************************************************************
  REAL(SP) GAUSSDER
  REAL(SP) , INTENT(IN) :: t
  GAUSSDER =  - EXP(0.5)*(t-offset_pulse)/sigma_pulse * EXP(-(t-offset_pulse)**2/(2.0_SP*sigma_pulse**2))  
!  write(fileout,*) 'gaussder(t)',gaussder, 't=',t
END FUNCTION GAUSSDER


FUNCTION D_BY_DT_GAUSSDER(t)
  USE nrtype
   USE output_error
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the time derivative of the Gaussian derivative pulse 
! above at time t.
!*******************************************************************************
  REAL(SP) D_BY_DT_GAUSSDER
  REAL(SP) , INTENT(IN) :: t
  D_BY_DT_GAUSSDER =  - EXP(0.5)/sigma_pulse*( 1 - ((t-offset_pulse)/sigma_pulse)**2) * & 
                      EXP(-(t-offset_pulse)**2/(2.0_SP*sigma_pulse**2))
!  write(fileout,*) 'D_t_gaussder(t)',d_by_dt_gaussder, 't=',t
END FUNCTION D_BY_DT_GAUSSDER


FUNCTION D2_BY_DT2_GAUSSDER(t)
  USE nrtype
   USE output_error
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the second time derivative of the Gaussian derivative pulse 
! above at time t.
!*******************************************************************************
  REAL(SP) D2_BY_DT2_GAUSSDER
  REAL(SP) , INTENT(IN) :: t
  D2_BY_DT2_GAUSSDER =  + EXP(0.5)/(sigma_pulse**3)*(t-offset_pulse)*( 3.0_SP - ((t-offset_pulse)/sigma_pulse)**2) * & 
                      EXP(-(t-offset_pulse)**2/(2.0_SP*sigma_pulse**2))
!  write(fileout,*) 'D2_t2_gaussder(t)',d2_by_dt2_gaussder, 't=',t
END FUNCTION D2_BY_DT2_GAUSSDER


FUNCTION D3_BY_DT3_GAUSSDER(t)
  USE nrtype
   USE output_error
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the third time derivative of the Gaussian derivative pulse 
! above at time t.
!*******************************************************************************
  REAL(SP) D3_BY_DT3_GAUSSDER
  REAL(SP) , INTENT(IN) :: t
  D3_BY_DT3_GAUSSDER =  + EXP(0.5)/(sigma_pulse**3) * & 
                       ( 3.0_SP - 6.0_SP*((t-offset_pulse)/sigma_pulse)**2 + ((t-offset_pulse)/sigma_pulse)**4 ) * & 
                       EXP(-(t-offset_pulse)**2/(2.0_SP*sigma_pulse**2))
END FUNCTION D3_BY_DT3_GAUSSDER

! Add additional time sources and their analytical derivatives here....


END MODULE TD_source




