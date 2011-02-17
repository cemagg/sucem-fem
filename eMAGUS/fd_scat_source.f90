MODULE FD_scat_source
  USE nrtype
  USE scattering_analysis_data
  IMPLICIT NONE
!*******************************************************************************
! This module contains datastractures and routines
! related to the frequency domain source for scattering computations. 
!
! Parts of the routines are based on the those in TD_source, developed
! by DBD during research stay at TU Delft, 2003. 
! Created: 16 Dec 2003, DBD.
!*******************************************************************************
   SAVE

   REAL(SP) :: t_pw,p_pw,e_pw                             ! Polarization properties of the incident plane wave     
                                                          ! theta,phi,eta (all radians)  
   REAL(SP) :: mag_pw,phase_pw                            ! Magnitude [V/m] and phase (radians) of incident plane wave     
 

INTERFACE FD_INC_FIELD
  MODULE PROCEDURE FD_INC_FIELD
END INTERFACE


INTERFACE FD_INC_FIELD_SETUP 
  MODULE PROCEDURE FD_INC_FIELD_SETUP 
END INTERFACE



CONTAINS 


FUNCTION FD_INC_FIELD(x,y,z,field_type,Y_char)
  USE boundary_conditions
  USE frequency_data
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the field a point (x,y,z) for a monochromatic uniform 
! plane wave, wavenumber k. The wave is characterized by propagation angles theta (t_pw) and phi (phi_pw),
! polarization angle eta (t_pw) and magnitude mag_pw.
!
! The routine returns either the electric field, or the accompanying magnetic field WITHOUT the correct
! Yc scaling - this cancels in some cases - or includes scaling by optional parameter Y_char to return the 
! correct magnitude of the H field if required. 
!*******************************************************************************
  REAL(SP), INTENT(IN) :: x,y,z                        ! Coordinates of field point. 
  COMPLEX(SPC), DIMENSION(3) :: FD_INC_FIELD           ! x, y and z components of field at point . 
  CHARACTER, INTENT(IN) :: field_type                  ! Electric or magnetic (without Yc term) field

  REAL(SP) kdotr 
  REAL(SP), INTENT(IN), OPTIONAL :: Y_char             ! Characteristic admittance.

  kdotr = FD_INC_FIELD_PHASEFACTOR(k0,x,y,z)                          ! 

  SELECT CASE (field_type) 
  CASE ('E','e')
    FD_INC_FIELD = mag_pw * EXP(j*phase_pw) * FD_E_INC_POLARIZATION() * EXP(-j*kdotr)
  CASE ('H','h')
    FD_INC_FIELD = mag_pw * EXP(j*phase_pw) * FD_H_INC_POLARIZATION() * EXP(-j*kdotr)
	IF (PRESENT(Y_char)) THEN
      FD_INC_FIELD  = FD_INC_FIELD * Y_char
	END IF
  CASE DEFAULT
	STOP 'Unimplemented field type'
  END SELECT 
END FUNCTION FD_INC_FIELD
!******************************************************************************

SUBROUTINE FD_INC_FIELD_SETUP 
  USE boundary_conditions
  USE frequency_data
  USE geometry
  USE inc_field
  USE nrtype
  USE output_error
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine sets up the incident field. 
!*******************************************************************************
! Original routine: D B Davidson, April 2003.
! Extended:         DBD, 20 May 2003. 
!*******************************************************************************
   
  REAL(SP), DIMENSION(3) :: polarize
 
  IF(NUM_A0_cards.NE.1) THEN
    CALL ERROR_FEMFEKO(1,4093)
  END IF
  mag_pw   =  A0data(1)%EiR1
  phase_pw =  A0data(1)%EiR2 * D2R
  t_pw     =  A0data(1)%EiR3 * D2R
  p_pw     =  A0data(1)%EiR4 * D2R
  e_pw     =  A0data(1)%EiR5 * D2R


  WRITE(FILEOUT,'(//,16X,(A))') 'EXCITATION BY INCIDENT PLANE ELECTROMAGNETIC WAVE'

  WRITE(FILEOUT,'(/,1X,(A))') 'Number of excitation:   N =      1 (only one permitted)'
  WRITE(FILEOUT,'(1X,(A),T26,(A),G11.5)') 'Frequency in Hz:','FREQ = ',frequency
  WRITE(FILEOUT,'(1X,(A),T26,(A),G11.5)') 'Wavelength in m:','LAMBDA = ',c_0/frequency
!  IF(TD_SOURCE_TYPE.EQ.1) THEN 
!    WRITE(FILEOUT,'(1X,(A))') 'Differentiated Gaussian'
!    WRITE(FILEOUT,'(1X,(A))') 'pulse used'
!    WRITE(FILEOUT,'(1X,(A),T26,(A),G12.4)') 'Pulse parameters:','SIGMA =',sigma_pulse
!  ELSE 
!    STOP 'Unknown pulse type'
!  ENDIF
!  WRITE(FILEOUT,'(1X,T26,(A),G12.4)')     'PULSE DELAY = ',offset_pulse
  WRITE(FILEOUT,'(1X,(A),T26,(A),F7.2,2X,(A),F7.2)') 'Direction of incidence:','THETA =',t_pw*r2d,'PHI =',p_pw*r2d 
  WRITE(FILEOUT,'(1X,(A),T26,(A))')       'Polarisation:','LINEAR (only)'
  WRITE(FILEOUT,'(1X,(A),T26,(A))')       'Axial ratio:','V   =    0.0000 (only)'
  WRITE(FILEOUT,'(1X,(A),T26,(A),F7.2)')  'Dir. of polarisation:','ETA =',e_pw*r2d
  WRITE(FILEOUT,'(1X,(A),T26,(A),G12.5)') 'Direction of propag.','BETA0X =',-SIN(t_pw)*COS(p_pw)
  WRITE(FILEOUT,'(1X,T26,(A),G12.5)')     'BETA0Y =',-SIN(t_pw)*SIN(p_pw)
  WRITE(FILEOUT,'(1X,T26,(A),G12.5)')     'BETA0Z =',-COS(t_pw)
  polarize = FD_E_INC_POLARIZATION()
  WRITE(FILEOUT,'(1X,(A),T26,(A),F12.6,T49,(A),F8.2)')  'Field strength in V/m:','|E0X| = ',mag_pw*polarize(1), & 
                                                        'ARG(E0X) = ',phase_pw*R2D
  WRITE(FILEOUT,'(1X,(A),T26,(A),F12.6,T49,(A),F8.2)')  '(Phase in deg.)','|E0Y| = ',mag_pw*polarize(2),'ARG(E0Y) = ',phase_pw*R2D
  WRITE(FILEOUT,'(1X,(A),T26,(A),F12.6,T49,(A),F8.2)')  ' ','|E0Z| = ',mag_pw*polarize(3),'ARG(E0Z) = ',phase_pw*R2D


  IF (ABCs(1)%r.GT.0.0_SP) THEN
    IF(ABCs(1)%type.EQ.1) THEN
      WRITE(FILEOUT,'(//,16X,(A),F7.2)') 'FIRST ORDER ABC APPLIED AT EXTERNAL BOUNDARY, AT RADIUS',ABCs(1)%r
    ELSE IF(ABCs(1)%type.EQ.-2) THEN
      WRITE(FILEOUT,'(//,16X,(A),F7.2)') 'SECOND ORDER ABC (CURL TERM ONLY) APPLIED AT EXTERNAL BOUNDARY, AT RADIUS',ABCs(1)%r
    ELSE IF(ABCs(1)%type.EQ.2) THEN
      WRITE(FILEOUT,'(//,16X,(A),F7.2)') 'SECOND ORDER ABC (BOTH CURL AND DIV TERM) APPLIED AT EXTERNAL BOUNDARY, AT RADIUS',&
	        ABCs(1)%r
    ELSE IF(ABCs(1)%type.EQ.3) THEN
      WRITE(FILEOUT,'(//,16X,(A),F7.2)') & 
	       'SECOND ORDER ABC (DIV TERM EXPERIMENTAL FORMULATION) APPLIED AT EXTERNAL BOUNDARY, AT RADIUS',ABCs(1)%r
    ELSE
      STOP 'IE: FD_INC_FIELD_SETUP'
    END IF
  ELSE 
    IF(ABCs(1)%type.EQ.1) THEN
      WRITE(FILEOUT,'(//,16X,(A))') 'FIRST ORDER ABC APPLIED AT EXTERNAL BOUNDARY'
    ELSE IF(ABCs(1)%type.EQ.-2) THEN
      WRITE(FILEOUT,'(//,16X,(A))') 'SECOND ORDER ABC (CURL TERM ONLY) APPLIED AT EXTERNAL BOUNDARY'
    ELSE IF(ABCs(1)%type.EQ.2) THEN
      WRITE(FILEOUT,'(//,16X,(A))') 'SECOND ORDER ABC (BOTH CURL AND DIV TERM) APPLIED AT EXTERNAL BOUNDARY'
    ELSE IF(ABCs(1)%type.EQ.3) THEN
      WRITE(FILEOUT,'(//,16X,(A))') 'SECOND ORDER ABC (DIV TERM EXPERIMENTAL FORMULATION) APPLIED AT EXTERNAL BOUNDARY'
    ELSE
      STOP 'IE: FD_INC_FIELD_SETUP'
    END IF
  END IF


  SCAT_FIELD_TEST: IF(SCAT_FIELD) THEN
     WRITE(FILEOUT,'(//,16X,(A),I3)') & 
          'SCATTERED FIELD FORMULATION, WITH HARD-CODED PEC SCATTERER. HOMOGENEOUS EMBEDDING MATERIAL NUMBER: ',HOMOG_MEDIUM
     IF(PEC_SPHERE_SCAT) THEN
        WRITE(FILEOUT,'(//,16X,(A),G12.5)') 'PEC SCATTERER FOUND BASED ON RADIUS SEARCH, RADIUS = ', sph_radius
        IF(CURVILINEAR) THEN
           WRITE(FILEOUT,'(//,16X,(A),G12.5)') 'SPHERICAL PEC SCATTERER TREATED WITH CURVILINEAR ELEMENTS, RADIUS = ', sph_radius
           WRITE(FILEOUT,'(//,16X,(A),G12.5)') 'CODE NOT IMPLEMENTED YET FOR SCATTERED FIELD FORMULATION: EXECUTION STOPPED'
           STOP 
        END IF
     ELSE IF(.NOT.PEC_SPHERE_SCAT.AND.CURVILINEAR) THEN
        WRITE(FILEOUT,'(//,16X,(A),G12.5)') &
             'ONLY SPHERICAL PEC SCATTERERS CAN BE TREATED WITH CURVILINEAR ELEMENTS: ERROR IN INPUT FILE'
        STOP
     ELSE IF(PENETRABLE_SCAT) THEN
        WRITE(FILEOUT,'(//,16X,(A))') 'PENETRABLE SCATTERER NOT AVAILABLE IN SCATTERED FIELD FORMULATION: ERROR IN INPUT FILE'
        STOP
     ELSE                       
        WRITE(FILEOUT,'(//,16X,(A))') 'PEC SCATTERER FOUND BASED ON LABEL SEARCH'
     END IF
     IF (TEST_INTERPOLATE_FIELD) THEN
        WRITE(FILEOUT,'(//,16X,(A))') 'WARNING: SPECIAL TEST MODE ACTIVATED. RESULTS ARE NEGATIVE OF INTERPOLATED INCIDENT FIELD.'
     END IF
  ELSE
     WRITE(FILEOUT,'(//,16X,(A))') 'TOTAL FIELD FORMULATION.' 
     IF(PEC_SPHERE_SCAT.AND.CURVILINEAR) THEN
        WRITE(FILEOUT,'(//,16X,(A))') 'PEC SCATTERER FOUND BASED ON LABEL SEARCH'
        WRITE(FILEOUT,'(//,16X,(A),G12.5)') 'SPHERICAL PEC SCATTERER TREATED WITH CURVILINEAR ELEMENTS, RADIUS = ', sph_radius
     ELSE IF(.NOT.PEC_SPHERE_SCAT.AND.CURVILINEAR) THEN
        WRITE(FILEOUT,'(//,16X,(A),G12.5)') &
             'ONLY SPHERICAL PEC SCATTERERS CAN BE TREATED WITH CURVILINEAR ELEMENTS: ERROR IN INPUT FILE'
        STOP
     ELSE IF(PEC_GENERAL_SCAT) THEN
        WRITE(FILEOUT,'(//,16X,(A))') 'PEC SCATTERER FOUND BASED ON LABEL SEARCH'
     ELSE IF(PENETRABLE_SCAT) THEN
        WRITE(FILEOUT,'(//,16X,(A))') 'PENETRABLE SCATTERER FOUND BASED ON LABEL SEARCH'
     END IF

  END IF SCAT_FIELD_TEST


END SUBROUTINE FD_INC_FIELD_SETUP


!******************************************************************************


FUNCTION FD_INC_FIELD_PHASEFACTOR(k,x,y,z) 
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the phase factor, \vec{k} \cdot \vec{r}, appearing in 
! the plane wave expression e^{-j\vec{k} \cdot \vec{r}}. 
!
!*******************************************************************************
  
  REAL(SP), INTENT(IN) :: k,x,y,z                        ! Wavenumber; coordinates of field point. 
  REAL(SP) FD_INC_FIELD_PHASEFACTOR                      ! See above. 
  REAL(SP), DIMENSION(3) :: k_inc                        ! Unit incident wavevector 
  REAL(SP), DIMENSION(3) :: field_point                  ! Point (x,y,z)
  INTEGER (I4B) ABC_num,ii

  field_point(1) = x
  field_point(2) = y
  field_point(3) = z


  ! compute k dot r - k as in above; r assuming some origin. Make origin 0 for present.....

  k_inc(1) = - k*sin(t_pw)*cos(p_pw) 
  k_inc(2) = - k*sin(t_pw)*sin(p_pw) 
  k_inc(3) = - k*cos(t_pw)

  FD_INC_FIELD_PHASEFACTOR = DOT_PRODUCT(k_inc,field_point)
END FUNCTION FD_INC_FIELD_PHASEFACTOR


!******************************************************************************


FUNCTION FD_E_INC_POLARIZATION()  
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the Cartesian components for the incident 
! frequency-domain plane wave polarization. (Identical to TD code). 
!*******************************************************************************
  
  REAL(SP), DIMENSION(3) :: FD_E_INC_POLARIZATION        ! x, y and z components planewave polarized as spedified. 

  FD_E_INC_POLARIZATION(1) = -COS(e_pw)*COS(t_pw)*COS(p_pw)-SIN(e_pw)*SIN(p_pw)
  FD_E_INC_POLARIZATION(2) = -COS(e_pw)*COS(t_pw)*SIN(p_pw)+SIN(e_pw)*COS(p_pw)
  FD_E_INC_POLARIZATION(3) = COS(e_pw)*SIN(t_pw)
END FUNCTION FD_E_INC_POLARIZATION


FUNCTION FD_H_INC_POLARIZATION()  ! (t_pw,p_pw,e_pw)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the Cartesian components of the associated
! magnetic field for the incident frequency-domain plane wave polarization. 
! (Identical to TD code).
!*******************************************************************************
  
  REAL(SP) e_pwH
  REAL(SP), DIMENSION(3) :: FD_H_INC_POLARIZATION        ! x, y and z components planewave polarized as spedified. 

  e_pwH = e_pw + PI/2.0_SP
  FD_H_INC_POLARIZATION(1) = -COS(e_pwH)*COS(t_pw)*COS(p_pw)-SIN(e_pwH)*SIN(p_pw)
  FD_H_INC_POLARIZATION(2) = -COS(e_pwH)*COS(t_pw)*SIN(p_pw)+SIN(e_pwH)*COS(p_pw)
  FD_H_INC_POLARIZATION(3) = COS(e_pwH)*SIN(t_pw)
END FUNCTION FD_H_INC_POLARIZATION

!******************************************************************************




END MODULE FD_scat_source




