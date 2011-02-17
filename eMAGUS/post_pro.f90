! During SGI port, approx line 1712, temporarily edited call to 
! SPSORT. DBD 8 Feb 05

! 2002-05-07: Added ERM_LOCAL_SOLVE. MMB.
! Last changed 1 May 2002
! New routine GW_OUTPUT_FIELDS added. 

SUBROUTINE GW_OUTPUT_S_PARAMS(freqcount)
  USE boundary_conditions
  USE math_tools, ONLY: PHASE
  USE frequency_data
  USE geometry
  USE gw_data
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Output the S parameters. 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: freqcount
  REAL(SP) lambda,lambda_g                  ! Free space and guide wavelengths
  REAL(SP) sigma                            ! Running sum.
  INTEGER(I4B) :: iport,jport               ! counters

  lambda   = c_0/frequency
  lambda_g = lambda/SQRT(1-(lambda/(2.0_SP*ports(1)%a))**2) 
  WRITE(FILEOUT,'(//,20X,A,/)') 'GUIDED WAVE RESULTS'
  WRITE(FILEOUT,'(1X,A,I8/)') & 
    'S parameters and conservation of energy check'

  WRITE(FILEOUT,'(1X,A,2X,A))') & 
       'Frequency','Guide wavelength'
  WRITE(FILEOUT,'(3X,A,10X,A)') '[GHz]','[mm]'
  WRITE(FILEOUT,'(2(G11.5,2X))') frequency/1.0E+9_SP, lambda_g*1.0E3_SP 
! Following contains MATLAB work-arounds, should be re-formatted for 
! full WinFEKO integration.
  WRITE(FILEOUT,'(/20X,A)') 'Complex valued format: Re, Im. MATLAB format'
! Note that freqcount runs from 0, an offset of 1 is needed. 
  WRITE(FILEOUT,'(A,I4,A)') 'S_FEM(:,:,',freqcount+1,') = ['
  DO iport = 1,num_ports
    DO jport = 1,num_ports
      WRITE(FILEOUT,'(10(E11.5,A,E11.5,4X))',ADVANCE='NO') REAL(S_params(iport,jport)),&
               '+j*',AIMAG(S_params(iport,jport))
    END DO
    WRITE(FILEOUT,'(A)') ' ' ! To advance line.  
  END DO
  WRITE(FILEOUT,'(A)') '];'

  WRITE(FILEOUT,'(/20X,A)') 'Magnitude/phase (degrees) format' 
  DO iport = 1,num_ports
    DO jport = 1,num_ports
      WRITE(FILEOUT,'(10(G11.5,2X))',ADVANCE='NO') ABS(S_params(iport,jport)),&
                                                   PHASE(S_params(iport,jport))
    END DO
    WRITE(FILEOUT,'(A)') ' ' ! To advance line.  
  END DO

  WRITE(FILEOUT,'(/20X,A)') 'Sum of |S_ij|^2 - by column'
  DO jport = 1,num_ports
    sigma = 0.0_SP
    DO iport = 1,num_ports
      sigma = sigma + ABS(S_params(iport,jport))**2
    END DO
      WRITE(FILEOUT,'(G11.5,2X))',ADVANCE='NO') sigma
  END DO
  WRITE(FILEOUT,'(A)') ' ' ! To advance line.  

END SUBROUTINE GW_OUTPUT_S_PARAMS
!*******************************************************************************

SUBROUTINE GW_OUTPUT_FIELDS(freqcount)
  USE boundary_conditions
  USE frequency_data
  USE gw_data
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Output some headers for near fields when using GW analysis.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: freqcount
  REAL(SP) lambda,lambda_g                  ! Free space and guide wavelengths

  lambda   = c_0/frequency
  lambda_g = lambda/SQRT(1-(lambda/(2.0_SP*ports(1)%a))**2) 
  WRITE(FILEOUT,'(//,20X,A,/)') 'GUIDED WAVE RESULTS'
  WRITE(FILEOUT,'(1X,A,I8/)') & 
    'Near field results'

  WRITE(FILEOUT,'(1X,A,2X,A))') & 
       'Frequency','Guide wavelength'
  WRITE(FILEOUT,'(3X,A,10X,A)') '[GHz]','[mm]'
  WRITE(FILEOUT,'(2(G11.5,2X))') frequency/1.0E+9_SP, lambda_g*1.0E3_SP 
  WRITE(FILEOUT,'(A)') ' ' ! To advance line.  
END SUBROUTINE GW_OUTPUT_FIELDS
!*******************************************************************************


SUBROUTINE OUTPUT_FE_RESULTS(feltyp,n_x,n_y,n_z,felkor,x0,y0,z0,delta_x,delta_y,delta_z,i_mode)
  USE fd_scat_source
  USE feminterface, ONLY: FEM_FIELDCALC
  USE cbaa_sys, ONLY: CBAA_FIELDCALC
  USE math_tools, ONLY: PHASE
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE td_source
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Write the data requested by the FEdata(cardnumber) record to the output file, 
! according to the type of analysis. (CBAA and GW are treated the same, but EIG 
! requires extra attention.)
! Extended to only write block headers once for TD analysis.
! Extended to compute total field in TD analysis case when needed.
! Extended DBD 12 Aug 04 to compute either scattered or total fields for FD formulation.
! 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: feltyp,felkor                 ! Field type and co-ordinate system to be used
  INTEGER(I4B), INTENT(IN) :: n_x,n_y,n_z                   ! numbers of points in the various directions
  REAL(SP), INTENT(IN) :: x0,y0,z0,delta_x,delta_y,delta_z  ! positional starting point and increments
  INTEGER(I4B), INTENT(IN), OPTIONAL :: i_mode              ! Eigenmode to evaluate in the EIG case
  LOGICAL(LGT), SAVE :: first_call = .true.

  INTEGER(I4B) :: xcount,ycount,zcount           ! counters
  REAL(SP) :: x_ob,y_ob,z_ob                     ! FE evaluation co-ordinates
  COMPLEX(SPC), DIMENSION(3) :: E_value,H_value  ! Field values (xyz)
  REAL(SP), DIMENSION(3) :: Re_E_inc             ! Incident field values (xyz)
  COMPLEX(SPC), DIMENSION(3) :: E_inc            ! Incident field values (xyz)

  ! Write the data block header:

  SELECT CASE (feltyp)
  CASE (-8,8)
    IF(.NOT.TD_ANALYSIS) THEN
      IF(FD_SCAT_ANALYSIS) THEN
   	    IF (feltyp.EQ.-8) THEN
		  IF (.NOT.TEST_INTERPOLATE_FIELD) THEN
            WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE SCATTERED ELECTRIC FIELD STRENGTH in V/m'
          ELSE 
            WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE INTERPOLATED INCIDENT ELECTRIC FIELD STRENGTH in V/m'
          END IF
	    ELSE 
          WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE TOTAL ELECTRIC FIELD STRENGTH in V/m'
        END IF 
	  ELSE 
        WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE ELECTRIC FIELD STRENGTH in V/m'
        WRITE (FILEOUT,'(/,20X,A)')  '      in the FEM/free space region(s)'
        WRITE (FILEOUT,'(/,A20,19X,3(A13,13X))') 'LOCATION','EX','EY','EZ'
        WRITE (FILEOUT,'(9(A13))') 'X/m','Y/m','Z/m','magn.','phase','magn.','phase','magn.','phase'
      END IF
    ELSE IF (TD_ANALYSIS.AND.first_call) THEN
	  IF (feltyp.EQ.-8) THEN
        WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE SCATTERED ELECTRIC FIELD STRENGTH in V/m'
	  ELSE 
        WRITE (FILEOUT,'(//,20X,A)') 'VALUES OF THE TOTAL ELECTRIC FIELD STRENGTH in V/m'
      END IF 
      WRITE (FILEOUT,'(/,20X,A)')  '      in the FEM/free space region(s)'
      WRITE (FILEOUT,'(/,A10,3X,A20,19X,3(A13))') 'TIME STEP', 'LOCATION','EX','EY','EZ'
      WRITE (FILEOUT,'(7(A13))') '* delta_t', 'X/m','Y/m','Z/m','value', 'value',' value'
	  first_call = .false.
	END IF
  CASE DEFAULT
    STOP 'IE: Unimplemented/Illegal field type request'
  END SELECT
  

  ! Write out the requested field values:
  DO zcount = 0,n_z-1
    z_ob = z0 + REAL(zcount)*delta_z ! observation point
    DO ycount = 0,n_y-1
      y_ob = y0 + REAL(ycount)*delta_y ! observation point
      DO xcount = 0,n_x-1
        x_ob = x0 + REAL(xcount)*delta_x ! observation point
                      
        IF (REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS) &
                           CALL FEM_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value,i_mode)
        IF (CBAA_ANALYSIS) CALL CBAA_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value)
        IF (GW_ANALYSIS)   CALL FEM_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value)
! Added DBD 17 Dec 03
        IF (FD_SCAT_ANALYSIS)   CALL FEM_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value)
! End added DBD 17 Dec 03        
		IF (TD_ANALYSIS)   CALL FEM_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value)
                    
        SELECT CASE(feltyp)
        CASE(-8) ! Electric field
		  IF(TD_ANALYSIS.AND.SCAT_FIELD) THEN 
            WRITE (FILEOUT,'(1X,I8,4X,6(1X,E12.5))') timestep, x_ob, y_ob, z_ob,     &
              REAL(E_value(1)), &
              REAL(E_value(2)), &
              REAL(E_value(3))
		  ELSE IF(FD_SCAT_ANALYSIS.AND.SCAT_FIELD) THEN ! Scattered field computed and required.
            WRITE (FILEOUT,'(9(1X,E12.5))') x_ob, y_ob, z_ob,     &
              ABS(E_value(1)), PHASE(E_value(1)), &
              ABS(E_value(2)), PHASE(E_value(2)), &
              ABS(E_value(3)), PHASE(E_value(3))
		  ELSE IF(FD_SCAT_ANALYSIS.AND..NOT.SCAT_FIELD) THEN ! Total field computed, scattered field required.
            E_inc=FD_INC_FIELD(x_ob,y_ob,z_ob,'E')
			E_value = E_value-E_inc ! Find scattered field.			  			  
            WRITE (FILEOUT,'(9(1X,E12.5))') x_ob, y_ob, z_ob,     &
              ABS(E_value(1)), PHASE(E_value(1)), &
              ABS(E_value(2)), PHASE(E_value(2)), &
              ABS(E_value(3)), PHASE(E_value(3))
		  ELSE 
            CALL ERROR_FEMFEKO(0,4903)
		  END IF 
        CASE(8) ! Electric field
		  IF(TD_ANALYSIS) THEN 
            IF (.NOT.SCAT_FIELD) THEN ! FEM fields are already total fields, can be written directly. 
              WRITE (FILEOUT,'(1X,I8,4X,6(1X,E12.5))') timestep,x_ob, y_ob, z_ob,     &
                REAL(E_value(1)), &
                REAL(E_value(2)), &
                REAL(E_value(3))
            ELSE ! FEM fields are scattered fields, add incident fields to each component. 
              Re_E_inc=TD_INC_FIELD(x_ob,y_ob,z_ob,'E')			  			  
              WRITE (FILEOUT,'(1X,I8,4X,6(1X,E12.5))') timestep, x_ob, y_ob, z_ob,     &
                REAL(E_value(1))+Re_E_inc(1), &
                REAL(E_value(2))+Re_E_inc(2), &
                REAL(E_value(3))+Re_E_inc(3)
            END IF
		  ELSE IF(FD_SCAT_ANALYSIS.AND.SCAT_FIELD) THEN ! Scattered  field computed, total field required.
            E_inc=FD_INC_FIELD(x_ob,y_ob,z_ob,'E')
			E_value = E_value+E_inc ! Find total field.			  			  
            WRITE (FILEOUT,'(9(1X,E12.5))') x_ob, y_ob, z_ob,     &
              ABS(E_value(1)), PHASE(E_value(1)), &
              ABS(E_value(2)), PHASE(E_value(2)), &
              ABS(E_value(3)), PHASE(E_value(3))
		  ELSE ! Total field computed and required.
            WRITE (FILEOUT,'(9(1X,E12.5))') x_ob, y_ob, z_ob,     &
              ABS(E_value(1)), PHASE(E_value(1)), &
              ABS(E_value(2)), PHASE(E_value(2)), &
              ABS(E_value(3)), PHASE(E_value(3))
		  END IF
        END SELECT
 
      END DO
    END DO
  END DO

END SUBROUTINE OUTPUT_FE_RESULTS
!*******************************************************************************

SUBROUTINE OUTPUT_FF_RESULTS(ntheta,nphi,theta0,phi0,dtheta,dphi,data_type,Einc_abs)
  USE cbaa_sys, ONLY: CBAA_FIELDCALC
  USE frequency_data
  USE math_tools, ONLY: PHASE
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Write requested FF data to the output file. The number of points, and their positions,
! is supplied as parameters, as well as the type of data to calculate.
! MMB 19 Jan 2001
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: ntheta,nphi        ! numbers of points in the theta an phi directions
  REAL(SP), INTENT(IN) :: theta0,phi0,           &
                          dtheta,dphi            ! starting values and increments in radians
  INTEGER(I4B), INTENT(IN) :: data_type          ! 1=RCS/2=directivity/3=gain/4=none
  REAL(SP), INTENT(IN), OPTIONAL :: Einc_abs     ! ABS of the incident plane wave

  INTEGER(I4B) :: cardk,cardl,i_mode,count1       ! Counters
  INTEGER(I4B) :: theta_total,phi_total           ! Counters
  REAL(SP) :: theta_ob,phi_ob                     ! angles describing the FF evaluation point
  REAL(SP) :: x_ob,y_ob,z_ob                      ! FE evaluation co-ordinates
  COMPLEX(SPC), DIMENSION(3) :: E_value,H_value   ! Field values (xyz)
  COMPLEX(SPC), DIMENSION(3) :: E_value_ptr       ! Field values (phi theta radius components)
  REAL(SP), DIMENSION(3) :: intensity_p_t_tot,    &
                            log10_intensity       ! Gain or directivity, phi, theta and total components
  REAL(SP) :: RCS                                 ! calculated Radar Cross Section
  REAL(SP) :: ff_radius                           ! radius where field will be evaluated
  REAL(SP), DIMENSION(3) :: unit_vec              ! unit vector
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)       ! The complex constant j

  ff_radius = 1000.0  ! assign a suitable value, this may change in order to accomodate
                      ! far-field criterion (> 2*D^2, see for example Stutzman & Thiele)

  ! Write the information header:
  WRITE (FILEOUT,'(///,20X,A)') 'VALUES OF THE SCATTERED ELECTRIC FIELD STRENGTH IN THE FAR FIELD in V'
  WRITE (FILEOUT,'(20X,A)')    '              Factor e^(-j*BETA*R)/R not considered'
  SELECT CASE(data_type)
  CASE (1) ! RCS
    WRITE (FILEOUT,'(/,3(A20,6X),2(A26,13X))') 'LOCATION','ETHETA','EPHI', &
                                               'scattering cross sect.','POLARISATION'
    WRITE (FILEOUT,'(6(A13),A26,13X,3(A13))') 'THETA','PHI','magn.','phase','magn.','phase', &
                                              'in m*m','axial r.','angle','direction'
  CASE (2) ! directivity
    WRITE (FILEOUT,'(/,3(A20,6X),2(A26,13X))') 'LOCATION','ETHETA','EPHI','directivity in dB', &
                                               'POLARISATION'
    WRITE (FILEOUT,'(12(A13))') 'THETA','PHI','magn.','phase','magn.','phase','vert.','horiz.', &
                                'total','axial r.','angle','direction'
  CASE (3) ! gain
    WRITE (FILEOUT,'(/,3(A20,6X),2(A26,13X))') 'LOCATION','ETHETA','EPHI','gain in dB','POLARISATION'
    WRITE (FILEOUT,'(12(A13))') 'THETA','PHI','magn.','phase','magn.','phase','vert.','horiz.', &
                                'total','axial r.','angle','direction'
  CASE (4) ! nothing
    WRITE (FILEOUT,'(/,3(A20,6X),2(A26,13X))') 'LOCATION','ETHETA','EPHI',' ','POLARISATION'
    WRITE (FILEOUT,'(6(A13),A39,3(A13))') 'THETA','PHI','magn.','phase','magn.','phase','', &
                                          'axial r.','angle','direction'
  END SELECT  

  ! Calculate the far field:
  DO cardk = 0,ntheta-1
    theta_ob = theta0 + REAL(cardk)*dtheta
    
    DO cardl = 0,nphi-1
      phi_ob   = phi0 + REAL(cardl)*dphi

      ! Calculate observation point in Cartesian co-ordinates:
      x_ob = ff_radius*SIN(theta_ob)*COS(phi_ob)
      y_ob = ff_radius*SIN(theta_ob)*SIN(phi_ob)
      z_ob = ff_radius*COS(theta_ob)

      ! Calculate the field values:
      IF (CBAA_ANALYSIS) THEN 
        CALL CBAA_FIELDCALC(x_ob,y_ob,z_ob,E_value,H_value)
      ELSE
        CALL ERROR_FEMFEKO(0,4900)
      END IF

      ! Multiply by the e^jkr/r factor, and convert to the phi and theta basis:
      unit_vec = (/ -SIN(phi_ob) , COS(phi_ob) , 0.0 /) ! phi unit vector
      E_value_ptr(1) = SUM(E_value*unit_vec)
      unit_vec = (/ COS(theta_ob)*COS(phi_ob) , COS(theta_ob)*SIN(phi_ob) , -SIN(theta_ob) /) ! theta unit vector
      E_value_ptr(2) = SUM(E_value*unit_vec)
      E_value_ptr(3) = (0.0,0.0) ! assume that the field is a plane wave, thus no r-component in the FF
      E_value_ptr = (ff_radius/(EXP(-cj*k0*ff_radius))) * E_value_ptr

      ! Calculate RCS/directivity/gain: (see Stutzman & Thiele)
      SELECT CASE (data_type)
      CASE (1) ! RCS:
        RCS = SUM(ABS(E_value_ptr(1:2))**2)/(Einc_abs**2)
      CASE (2) ! directivity:
        intensity_p_t_tot(1:2) = (1.0/Z_zero) * (ABS(E_value_ptr(1:2))**2)
      CASE (3) ! gain:
        intensity_p_t_tot(1:2) = (1.0/Z_zero) * (ABS(E_value_ptr(1:2))**2)
      END SELECT
      intensity_p_t_tot(3) = SUM(intensity_p_t_tot(1:2)) 

      ! Convert to dB in the directivity/gain case:
      DO count1 = 1,3
        IF ((intensity_p_t_tot(count1)).LE.0.0) THEN
          log10_intensity(count1) = -999.0
        ELSE
          log10_intensity(count1) = 10.0*LOG10(intensity_p_t_tot(count1))
        END IF
      END DO

      ! Write to the output file:
      WRITE (FILEOUT,'(6(1X,E12.5),A13)',ADVANCE='NO') (180.0/PI)*theta_ob, (180.0/PI)*phi_ob, &
                                                   ABS(E_value_ptr(2)), PHASE(E_value_ptr(2)), &
                                                   ABS(E_value_ptr(1)), PHASE(E_value_ptr(1))
      SELECT CASE (data_type)
      CASE (1)   ! RCS
        WRITE (FILEOUT,'(13X,1X,E12.5,13X)',ADVANCE='NO') RCS
      CASE (2:3) ! directivity/gain
        WRITE (FILEOUT,'(3(1X,E12.5))',ADVANCE='NO') log10_intensity(2), log10_intensity(1), &
                                                     log10_intensity(3)
      CASE (4)   ! nothing                                               
        WRITE (FILEOUT,'(39X)',ADVANCE='NO')
      END SELECT
      WRITE (FILEOUT,'(2(1X,E12.5),A13)') 0.0, 0.0, 'RIGHT'  

    END DO
  END DO

END SUBROUTINE OUTPUT_FF_RESULTS
!*******************************************************************************


SUBROUTINE FEM_FIELDCALC(xob,yob,zob,E_xyz,H_xyz,eigen_mode)
  USE basis_function, ONLY: EVALUATE_ELEMENTAL_FUNCTIONS
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TET, LOCAL_TO_GLOBAL_INDEX_PRE_TET
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE scattering_analysis_data
  IMPLICIT NONE
  !*******************************************************************************
  ! Calculates E and H field in the GW_ANALYSIS and in the TD_ANALYSIS case, 
  ! at the location specified. 
  ! Also used by CBAA_FIELDCALC if the observation point is not located in the
  ! Green's function region. Also used for calculation of the field in the 
  ! eigen analysis case, but then the optional parameter <eigen_mode> must also 
  ! be specified.
  !*******************************************************************************
  ! Last changed:
  !
  ! 28 March 2002: Extended to also support LT/LN and QT/QN elements. DBD.
  ! 02 April 2003: Extended to also support TD analysis. DBD.
  ! 07 Aug   2004: Extended to also support FD scattering analysis with prescribed dof's. DBD.
  !  
  !*******************************************************************************
  REAL(SP), INTENT(IN) :: xob,yob,zob
  COMPLEX(SPC), DIMENSION(3), INTENT(OUT) :: E_xyz,H_xyz
  INTEGER(I4B), INTENT(IN), OPTIONAL :: eigen_mode

  INTEGER(I4B) :: inside_el,ifunc,row,row_pre,kk
  REAL(SP), DIMENSION(3) :: r_vec              ! evaluation co-ordinates
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: func_values     ! Basis function values at the specified location
  INTEGER(I4B) :: f_dimension                  ! number of basis functions within an element
  COMPLEX(SPC) :: dof_value
  integer(i4b) :: tmpi
  ! Initialize:
  H_xyz = (0.0,0.0) 
  E_xyz = (0.0,0.0)

  inside_el = XYZ_TO_ELNUM(xob,yob,zob)
  IF (inside_el.EQ.0) RETURN ! thus the field is zero and not located within an element

  ! Choose the dimension of the elemental contributions:
  SELECT CASE (MAX_ORDER(inside_el))
  CASE (1) ! CT/LN and LT/QN
     IF (MIXED_ORDER(inside_el)) THEN ! CT/LN
        f_dimension = 6  ! 6 E1 in the element
     ELSE 
        f_dimension = 12  ! and additional 6 E2 in the element
     END IF
  CASE (2) ! LT/QN and QT/QN
     IF (MIXED_ORDER(inside_el)) THEN ! CT/LN
        f_dimension = 20 ! 6 E1, 6 E2, 4 F1, 4 F2 in the element
     ELSE 
        f_dimension = 30  ! and additional 6 E3 and 4 F3 in the element
     END IF
  CASE DEFAULT
     STOP 'IE: invalid hierarchal order in FEM_FIELDCALC.'
  END SELECT

  !write(fileout,'(//)') 
  !write(fileout,*) 'element #:', inside_el, 'r_vec:',  r_vec  

  r_vec = (/ xob , yob , zob /)
  CALL EVALUATE_ELEMENTAL_FUNCTIONS(inside_el,r_vec,0,func_values)
  !write(fileout,'(A)') 'Basis function values   x     y     z:'
  !do kk=1,ELEM_TET_MATRIX_SIZE
  !  write(fileout,'(3(F8.4,2X))') func_values(kk,:)
  !end do
  element_function_loop: DO ifunc = 1,f_dimension ! Correct code
     !element_function_loop: DO ifunc = 1,6 ! testcode
     row = LOCAL_TO_GLOBAL_INDEX_TET(inside_el,ifunc)
     row_pre = 0 ! Re-initialize on each loop.
     IF(.NOT.(FD_SCAT_ANALYSIS.AND.SCAT_FIELD)) THEN 
        IF (row.EQ.0) CYCLE element_function_loop ! not a free dof
     ELSE ! This should work for other types of analyis with prescribed dof's, but is presently
        ! only implemented for FD_SCAT_ANALYSIS.
        IF (row.EQ.0) THEN
           row_pre = LOCAL_TO_GLOBAL_INDEX_PRE_TET(inside_el,ifunc) ! See if it is a prescribed dof.
           IF (row_pre.EQ.0) THEN
              ! For scattering analysis, this may be an error.
              WRITE(FILEOUT,*)'Possible error in FEM_FIELDCALC with a dof which is neither free nor prescribed: row,row_pre',& 
                   row,row_pre
              CYCLE element_function_loop ! not a prescribed dof either
           END IF
        END IF
     END IF
     ! Assign the dof value according to the analysis type:
     IF (CBAA_ANALYSIS)        dof_value = x_vec_c(row)
     IF (GW_ANALYSIS)          dof_value = x_vec_c(row)
     ! Added DBD 18 Dec 2003	and 07 Aug 04
     IF (FD_SCAT_ANALYSIS) THEN
        IF(row.NE.0) THEN
           IF(TEST_INTERPOLATE_FIELD) WRITE (FILEOUT,*) 'Error in FEM_FIELDCALC when DEBUG_INTERPOLATE activated'
           !dof_value = 0 ! test code 
           dof_value = x_vec_c(row) ! correct code... put back in
        ELSE IF(row_pre.NE.0)  THEN
           dof_value = x_pre_vec_c(row_pre)
           !write(fileout,*) 'Prescribed dof ',row_pre, 'values:',dof_value
        ELSE
           STOP 'IE: Error in FEM_FIELDCALC'
        END IF
     END IF
     ! End added DBD 18 Dec 2003	and 07 Aug 04
     IF (REAL_EIGEN_ANALYSIS)  dof_value = CMPLX(eigenvectors(row,eigen_mode)) ! Make complex for convenience; imag. part remains zero
     ! Added DBD 02 Apr 2003	
     IF (TD_ANALYSIS)          dof_value = CMPLX(u_nplus1(row))                ! ditto.
     IF (COUPLED_TD_ANALYSIS) dof_value = CMPLX(coupled_td_x_vec(row))
     ! End added DBD 02 Apr 2003	
     E_xyz = E_xyz + dof_value*func_values(ifunc,1:3) ! Add every function's contribution
     !           H_xyz = H_xyz + tempvec2*edgelen*dof_value/                          &
     !                    (-cj*mu_0*mu_r(elements(inside_el)%material)*c_0*k0)
  END DO element_function_loop

END SUBROUTINE FEM_FIELDCALC
!*******************************************************************************


FUNCTION EVALUATE_ELECTRIC_FIELD(elem,s_vec,curl_power,elem_dof_values)
  USE basis_function
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TET
  USE geometry, ONLY: SIMPLEX_COEFFICIENTS
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine return the curl(curl FEM electric field) at simplex coordinates
! <s_vec> in element <elem>.
!
! 2002-03-07: Created. MMB
! 2002-05-09: Added optional parameter specification of elemental dofs. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem
  REAL(SP), DIMENSION(4), INTENT(IN) :: s_vec
  INTEGER(I4B), INTENT(IN) :: curl_power
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), OPTIONAL, INTENT(IN) :: elem_dof_values
  COMPLEX(SPC), DIMENSION(3) :: EVALUATE_ELECTRIC_FIELD
 
  INTEGER(I4B) :: ifunc,row ! counters
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: func_values ! basis function values at the location <s_vec>
  REAL(SP), DIMENSION(4,4) :: vertmat
  REAL(SP), DIMENSION(4) :: r_vec
  COMPLEX(SPC) :: dof_value

  ! Initialise:
  EVALUATE_ELECTRIC_FIELD = 0.0

  ! Calculate the xyz coords of the evaluation point:
  CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vertmat)
  r_vec = MATMUL(vertmat,s_vec)

  ! Calculate the basis function values:
  IF (PRESENT(elem_dof_values)) THEN
    CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),curl_power,func_values, &
	                                  MAX_PRESENT_ORDER,MAX_PRESENT_MIXED)
  ELSE
    CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),curl_power,func_values)
  END IF

  ! Add the contributions of all the basis functins to the total:
  FUNCTION_LOOP: DO ifunc = 1,ELEM_TET_MATRIX_SIZE

    ! Establish the corresponding dof number and dof value:
    IF (PRESENT(elem_dof_values)) THEN
      dof_value = elem_dof_values(ifunc)
    ELSE
	  row = LOCAL_TO_GLOBAL_INDEX_TET(elem,ifunc)
      IF (row.EQ.0) CYCLE FUNCTION_LOOP ! not a free dof
      IF (CBAA_ANALYSIS.OR.GW_ANALYSIS) dof_value = x_vec_c(row)
      IF (REAL_EIGEN_ANALYSIS) THEN
	    STOP 'IE: EIGEN_ANALYSIS not supported in EVALUATE_ELECTRIC_FIELD.'
	    ! dof_value = CMPLX(eigenvectors(row,eigen_mode))
      END IF
    END IF

    ! Add every function's contribution to the total:
    EVALUATE_ELECTRIC_FIELD = EVALUATE_ELECTRIC_FIELD + dof_value*func_values(ifunc,1:3)

  END DO FUNCTION_LOOP

END FUNCTION EVALUATE_ELECTRIC_FIELD
!*******************************************************************************


FUNCTION L_SQ_NORM_ELEMENTAL(elem,elem_dof_vec)
  USE basis_function
  USE feminterface, ONLY: EVALUATE_ELECTRIC_FIELD
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! This routine calculates the L^2 norm of the electric field on element <elem>.
! The last FEM solution is used, but optionally, the degrees of freedom can be
! supplied as an argument.
!
! 2002-05-09: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), OPTIONAL, INTENT(IN) :: elem_dof_vec
  REAL(SP) :: L_SQ_NORM_ELEMENTAL 

  INTEGER(I4B) :: iquad,num_qpts
  COMPLEX(SPC), DIMENSION(3) :: field_value
  REAL(SP) :: temp_real
  REAL(SP), DIMENSION(4) :: s_vec

  L_SQ_NORM_ELEMENTAL = 0.0 ! init
  num_qpts            = 11  ! init

  QUAD_LOOP: DO iquad = 1,num_qpts
    s_vec = cube_tet_rules(num_qpts)%rule(iquad,1:4)
    IF (PRESENT(elem_dof_vec)) THEN
	  field_value = EVALUATE_ELECTRIC_FIELD(elem,s_vec,0,elem_dof_vec)
    ELSE
	  field_value = EVALUATE_ELECTRIC_FIELD(elem,s_vec,0)
    END IF
	temp_real = DOT_PRODUCT(field_value,field_value)
    L_SQ_NORM_ELEMENTAL = L_SQ_NORM_ELEMENTAL + &
                          cube_tet_rules(num_qpts)%rule(iquad,5) * temp_real
  END DO QUAD_LOOP

  ! Weight with volume (final step in quadrature):
  L_SQ_NORM_ELEMENTAL = VOLUME(elem) * L_SQ_NORM_ELEMENTAL

  ! Take square root of the integral, to obtain the norm:
  L_SQ_NORM_ELEMENTAL = SQRT(L_SQ_NORM_ELEMENTAL)

END FUNCTION L_SQ_NORM_ELEMENTAL
!*******************************************************************************


FUNCTION RESIDUAL_ELEMENT_HOMOGENEOUS(elem)
  USE feminterface, ONLY: EVALUATE_ELECTRIC_FIELD
  USE frequency_data
  USE geometry
  USE nrtype
  USE quad_tables
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! This routine calculates the homogeneous wave equation contribution to the
! (elemental volume residual)^2 of element <elem>. It
! returns the integral of (homogeneous volume residual)^2.
!
! 2002-03-07: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem
  REAL(SP) :: RESIDUAL_ELEMENT_HOMOGENEOUS

  INTEGER(I4B) :: iquad,num_qpts
  COMPLEX(SPC), DIMENSION(3) :: residual_vec
  REAL(SP), DIMENSION(4) :: s_vec

  ! Initialise:
  RESIDUAL_ELEMENT_HOMOGENEOUS = 0.0
  
  ! Choose the number of elemental quadrature points:
  num_qpts = 11 ! 4th order rule

  ! Integrate:
  DO iquad = 1,num_qpts
    
	! Establish the simplex coords of the integration point:
	s_vec = cube_tet_rules(num_qpts)%rule(iquad,1:4)
    
	! Caculate the homogeneous resudual at the point:
	residual_vec = - (1.0/mu_r(elements(elem)%material))*EVALUATE_ELECTRIC_FIELD(elem,s_vec,2)
    residual_vec = residual_vec &
	               + (k0**(2.0))*eps_r(elements(elem)%material)*EVALUATE_ELECTRIC_FIELD(elem,s_vec,0)
    
	! Add to the total, weighted with the quadrature rule weight:
	RESIDUAL_ELEMENT_HOMOGENEOUS = RESIDUAL_ELEMENT_HOMOGENEOUS +              &
	                               cube_tet_rules(num_qpts)%rule(iquad,5)*     &
								   REAL(DOT_PRODUCT(residual_vec,residual_vec))
  END DO

  ! Final step of quadrature: scale by tet volume:
  RESIDUAL_ELEMENT_HOMOGENEOUS = VOLUME(elem)*RESIDUAL_ELEMENT_HOMOGENEOUS

END FUNCTION RESIDUAL_ELEMENT_HOMOGENEOUS
!*******************************************************************************


FUNCTION RESIDUAL_FACE(face_num)
  USE boundary_conditions
  USE coax_feed
  USE feminterface, ONLY: EVALUATE_ELECTRIC_FIELD
  USE gw_sys, ONLY: E_MN, K_Z_MN
  USE cbaa_sys, ONLY: CBAA_FIELDCALC_BI, COAX_MODE
  USE frequency_data
  USE geometry
  USE gw_data
  USE math_tools, ONLY: FIND_IN_LIST,CROSS_PRODUCT
  USE nrtype
  USE quad_tables
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the square of the L^2 norm of the face residual, as
! defined in MMB-PhD.
!
! 2002-03-09: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: face_num
  REAL(SP) :: RESIDUAL_FACE

  INTEGER(I4B) :: elem1,elem2                        ! elements sharing the face
  INTEGER(I4B) :: local_face1,local_face2            ! local face number to element 1
  REAL(SP), DIMENSION(3) :: normaltemp               ! normal vector pointing from elem1 to elem2
  COMPLEX(SPC), DIMENSION(3) :: normal12             ! normal vector pointing from elem1 to elem2
  INTEGER(I4B) :: iquad,num_qpts                     ! counters
  REAL(SP),DIMENSION(4) :: s_vec1,s_vec2             ! simplex coords in elem1 and elem2
  INTEGER(I4B),DIMENSION(3) :: facenodes1,facenodes2 ! local face nodes in elem1 and elem2
  COMPLEX(SPC), DIMENSION(3) :: residual_vec,tempvec
  INTEGER(I4B) :: port_num                           ! GW port number of the face
  COMPLEX(SPC) :: E0                                 ! complex size of the GW excitation
  REAL(SP) :: kzmn,k_coax                            ! GW propagation constant
  REAL(SP), DIMENSION(3) :: r_vec
  COMPLEX(SPC), DIMENSION(3) :: E_field,H_field

  ! Initialise:
  num_qpts = 7
  RESIDUAL_FACE = 0.0
  s_vec1 = 0.0
  s_vec2 = 0.0

  ! If the face is a PEC boundary, then it makes no contribution:
  IF (faces(face_num)%PEC) THEN
	RETURN
  END IF
 
  ! Establish the elements sharing the face:
  elem1 = faceconnectelem(face_num,1)
  elem2 = faceconnectelem(face_num,2)
  IF (elem1.EQ.0) STOP 'IE: invalid first element in RESIDUAL_FACE.'

  ! Establish the local face numbers:
  CALL FIND_IN_LIST(elements(elem1)%faces(1:4),4,face_num,local_face1)
  IF (local_face1.EQ.0) STOP 'IE: invalid local face in RESIDUAL_FACE.'
  IF (elem2.NE.0) THEN ! not a boundary face
    CALL FIND_IN_LIST(elements(elem2)%faces(1:4),4,face_num,local_face2)
    IF (local_face2.EQ.0) STOP 'IE: invalid local face in RESIDUAL_FACE.'
  END IF

  ! Calculate the face, normal vector, pointing away from the first element
  ! sharing the face:
  CALL FIND_IN_LIST(elements(elem1)%faces(1:4),4,face_num,local_face1)
  IF (local_face1.EQ.0) STOP 'IE: invalid local face in RESIDUAL_FACE.'
  normaltemp = ELEMENT_NORMAL_VECTOR(elem1,local_face1)
  normal12 = normaltemp ! convert to complex values

  ! Establish local face nodes of elem1 and elem2:
  facenodes1 = LOCAL_FACENODES(local_face1)
  IF (elem2.NE.0) THEN ! not a boundary face
    facenodes2 = LOCAL_FACENODES(local_face2)
  END IF

  
!!test::
IF (elem2.NE.0) THEN ! not a boundary face
IF (elements(elem1)%nodes(facenodes1(1)).NE.elements(elem2)%nodes(facenodes2(1))) THEN
STOP 'Groot kak!'
END IF
IF (elements(elem1)%nodes(facenodes1(2)).NE.elements(elem2)%nodes(facenodes2(2))) THEN
STOP 'Groot kak!'
END IF
IF (elements(elem1)%nodes(facenodes1(3)).NE.elements(elem2)%nodes(facenodes2(3))) THEN
STOP 'Groot kak!'
END IF
END IF


  ! Integrate over the face using quadrature:
  DO iquad = 1,num_qpts
    
	! Init at this quad point:
	residual_vec = ZERO_C

    ! Coordinates:
	s_vec1(facenodes1) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
    r_vec              = XYZ_COORDINATES(elem1,s_vec1)

	! Contribution from elem1:
	residual_vec = residual_vec + (1.0/mu_r(elements(elem1)%material)) * &
	                              CROSS_PRODUCT(normal12,EVALUATE_ELECTRIC_FIELD(elem1,s_vec1,1))
    ! Boundary contribution:
	IF (elem2.EQ.0) THEN
      
	  ! Waveguide port:
      IF (faces(face_num)%port) THEN
	    port_num = faces(face_num)%portnumber
        kzmn = K_Z_MN(k0,port_num,mode_m,mode_n)
        E0 = ports(port_num)%excitation
        tempvec = E0*E_MN(r_vec(1),r_vec(2),r_vec(3),port_num,mode_m,mode_n)
        tempvec = 2.0*tempvec - EVALUATE_ELECTRIC_FIELD(elem1,s_vec1,0)
	    tempvec = CROSS_PRODUCT(normal12,tempvec)
	    tempvec = CROSS_PRODUCT(normal12,tempvec)
        tempvec = -j*kzmn*tempvec
        residual_vec = residual_vec + tempvec
	  END IF

	  ! Coax port:
      IF (faces(face_num)%coax_aperture) THEN
	    port_num     = faces(face_num)%coaxnumber
        k_coax       = k0*SQRT(AXdata(port_num)%coax_eps*AXdata(port_num)%coax_mu)
        E0           = AXdata(port_num)%coax_Iabs *                               &
                       EXP(j*AXdata(port_num)%coax_Iphase*PI/180.0) *             &
			           (Z_zero/(2.0*PI*AXdata(port_num)%coax_a)) *                &
				        SQRT(AXdata(port_num)%coax_mu/AXdata(port_num)%coax_eps)
        tempvec      = E0*COAX_MODE(AXdata(port_num)%coaxcentre,AXdata(port_num)%coax_a,r_vec)
        tempvec      = 2.0*tempvec - EVALUATE_ELECTRIC_FIELD(elem1,s_vec1,0)
	    tempvec      = CROSS_PRODUCT(normal12,tempvec)
	    tempvec      = CROSS_PRODUCT(normal12,tempvec)
        tempvec      = -j*(k_coax/AXdata(port_num)%coax_mu)*tempvec
        residual_vec = residual_vec + tempvec
	  END IF

      ! CBAA face:
	  CBAA_FACE_TEST: IF (faces(face_num)%CBAA_aperture) THEN

        ! Avoiding the singularity:
        r_vec(3) = r_vec(3) + 0.01*lambda0
		CALL CBAA_FIELDCALC_BI(r_vec,E_field,H_field,elem1,25)
        r_vec(3) = r_vec(3) - 0.01*lambda0

        tempvec = -j*k0*Z_zero*H_field
        tempvec = CROSS_PRODUCT(normal12,tempvec)

        ! Minus, because it subtracts from the residual (see MMB PhD):
		residual_vec = residual_vec - tempvec
      END IF CBAA_FACE_TEST
	
	END IF

	! Contribution from elem2:
    IF (elem2.NE.0) THEN ! not a boundary face
      s_vec2(facenodes2) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
	  residual_vec = residual_vec - (1.0/mu_r(elements(elem2)%material)) * &
	                                CROSS_PRODUCT(normal12,EVALUATE_ELECTRIC_FIELD(elem2,s_vec2,1))
    END IF

    ! Add to total:
	RESIDUAL_FACE = RESIDUAL_FACE + &
	                quad_tri_rules(num_qpts)%rule(iquad,4) * &
					REAL(DOT_PRODUCT(residual_vec,residual_vec))
  END DO

  ! Final step in quadrature: scale by triangle area:
  RESIDUAL_FACE = FACE_AREA(elem1,local_face1)*RESIDUAL_FACE

END FUNCTION RESIDUAL_FACE
!*******************************************************************************


SUBROUTINE EVALUATE_RESIDUAL_INDICATORS
  USE adaptive_fem
  USE feminterface, ONLY: RESIDUAL_ELEMENT_HOMOGENEOUS, RESIDUAL_FACE
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine calculates all the elemental and facial residual. They are
! combined to obtain total, elemental residual error indicators. See MMB-PhD.
! The datastructures to store these quatities are also allocated.
! 
! 2002-03-10: Created. MMB.
!*******************************************************************************
  INTEGER(I4B) :: ielem,iface     ! counters
  REAL(SP) :: residual,diameter
  REAL(SP) :: alpha               ! weighting factor for relative face and element contributions

  ! Assign the relative contributions of the faces and elements:
  alpha = ADdata%relative_res ! (alpha) element and (1-alpha) face
 

  ! Allocate storage for residuals and indicators:
  ALLOCATE(residuals_elements(num_elements))
  ALLOCATE(residuals_faces(num_faces))
  ALLOCATE(error_indicators(num_elements))
  ALLOCATE(diameters_elements(num_elements))
  ALLOCATE(diameters_faces(num_faces))

  ! Evaluate element diameters:
  DO ielem = 1,num_elements
    diameters_elements(ielem) = MAXVAL(T_LENGTHS(ielem))
  END DO

  ! Evaluate face diameters:
  DO iface = 1,num_faces
    IF (faceconnectelem(iface,2).NE.0) THEN ! a shared face
	  diameters_faces(iface) = MAX(diameters_elements(faceconnectelem(iface,1)), &
	                               diameters_elements(faceconnectelem(iface,2)))
    ELSE ! a boundary face
	  diameters_faces(iface) = diameters_elements(faceconnectelem(iface,1))
	END IF
  END DO

  ! Evaluate elemental residuals squared:
  DO ielem = 1,num_elements
	PRINT *, 'Evaluating element residual:',ielem,'/',num_elements
	residuals_elements(ielem) = RESIDUAL_ELEMENT_HOMOGENEOUS(ielem)
  END DO

  ! Evaluate facial residuals squared:
  DO iface = 1,num_faces
	PRINT *, 'Evaluating face residual:',iface,'/',num_faces
    residuals_faces(iface) = RESIDUAL_FACE(iface)
  END DO

  ! Calculate elemental indicators:
  error_indicators = 0.0 ! initisalise
  DO ielem = 1,num_elements

	! Add the volume contribution:
	residual = residuals_elements(ielem)
	diameter = diameters_elements(ielem)
	error_indicators(ielem) = error_indicators(ielem) + alpha*(diameter**(2.0))*residual

	! Add the contributions of the 4 faces:
    DO iface = 1,4
      residual = residuals_faces(elements(ielem)%faces(iface))
	  diameter = diameters_faces(elements(ielem)%faces(iface))
	  error_indicators(ielem) = error_indicators(ielem) + 0.5*(1.0-alpha)*diameter*residual
	END DO

  END DO

END SUBROUTINE EVALUATE_RESIDUAL_INDICATORS
!*******************************************************************************


SUBROUTINE ASSIGN_ELEMENT_ORDERS
  USE adaptive_fem
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine assigns new hierarchal orders to every element according to its  
! error indicator value (if it falls within the frac_upgrade fraction of largest 
! values). The new order is specified by the flag <upgrade_type>.
!
! 2002-03-11: Created. MMB.
! 2002-03-21: Added frac_upgrade variable. MMB.
! 2002-04-15: Added upgrade_type flag. MMB.
!*******************************************************************************

  INTEGER(I4B) :: ielem
  REAL(SP) :: list_fraction
  INTEGER(I4B) :: assign_count

  ! random assignment
  !  assign_count = 0
  !  DO ielem = 1,num_elements
  !    IF ((REAL(ielem)*ADdata%frac_upgrade).GT.REAL(assign_count)) THEN
  !      elements(ielem)%order = elements(ielem)%order + 1
  !      assign_count = assign_count + 1
  !    END IF
  !  END DO

  DO ielem = 1,num_elements
    list_fraction = REAL(ielem)/REAL(num_elements)
	IF (list_fraction.GE.(1-ADdata%frac_upgrade)) THEN

      SELECT CASE (ADdata%upgrade_type)
      CASE (0) ! None (CT/LN)
        elements(sorted_order(ielem))%order = 1
        elements(sorted_order(ielem))%mixed = .TRUE.
      CASE (1) ! LT/LN
        elements(sorted_order(ielem))%order = 1
        elements(sorted_order(ielem))%mixed = .FALSE.
      CASE (2) ! LT/QN
        elements(sorted_order(ielem))%order = 2
        elements(sorted_order(ielem))%mixed = .TRUE.
      CASE (3) ! QT/QN
        elements(sorted_order(ielem))%order = 2
        elements(sorted_order(ielem))%mixed = .FALSE.
      CASE (4) ! Graded
		IF ((1.0-list_fraction).LT.(ADdata%frac_upgrade/3.0)) THEN ! 1st third becomes LT/LN
          elements(sorted_order(ielem))%order = 1
          elements(sorted_order(ielem))%mixed = .FALSE.
		ELSE IF ((1.0-list_fraction).LT.(2.0*ADdata%frac_upgrade/3.0)) THEN ! 2nd third becomes LT/QN
          elements(sorted_order(ielem))%order = 2
          elements(sorted_order(ielem))%mixed = .TRUE.
		ELSE ! 3rd third becomes QT/QN
          elements(sorted_order(ielem))%order = 2
          elements(sorted_order(ielem))%mixed = .FALSE.
		END IF
      END SELECT

	END IF
  END DO

END SUBROUTINE ASSIGN_ELEMENT_ORDERS
!*******************************************************************************


SUBROUTINE ASSIGN_ERROR_LEVEL_FLAGS
  USE adaptive_fem
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine assigns inter flags to every element according to its error 
! indicator. These labels are used to output the error distribution via
! a *.fek-file. This routine used a slatec-library routine. The library is available
! www.netlib.org.
!
! 2002-03-11: Created. MMB.
!*******************************************************************************

  INTEGER(I4B) :: ielem
  REAL(SP) :: list_fraction

  ALLOCATE(error_level_flags(num_elements))

  DO ielem = 1,num_elements
  
    list_fraction = REAL(ielem)/REAL(num_elements)
    
	IF      (list_fraction.LE.0.80) THEN
      error_level_flags(sorted_order(ielem)) = 1
    ELSE IF (list_fraction.LE.0.825) THEN
      error_level_flags(sorted_order(ielem)) = 1
    ELSE IF (list_fraction.LE.0.85) THEN
      error_level_flags(sorted_order(ielem)) = 1
    ELSE IF (list_fraction.LE.0.875) THEN
      error_level_flags(sorted_order(ielem)) = 1
    ELSE IF (list_fraction.LE.0.9) THEN
      error_level_flags(sorted_order(ielem)) = 1
    ELSE IF (list_fraction.LE.0.925) THEN
      error_level_flags(sorted_order(ielem)) = 2
    ELSE IF (list_fraction.LE.0.95) THEN!
      error_level_flags(sorted_order(ielem)) = 3
    ELSE IF (list_fraction.LE.0.975) THEN
      error_level_flags(sorted_order(ielem)) = 4
    ELSE 
      error_level_flags(sorted_order(ielem)) = 5
	END IF

! Temp code to get nice WinFEKO colours:
!	IF      (list_fraction.LE.0.8) THEN
!      error_level_flags(sorted_order(ielem)) = 1
!    ELSE 
!      error_level_flags(sorted_order(ielem)) = 8
!	END IF

  END DO

END SUBROUTINE ASSIGN_ERROR_LEVEL_FLAGS
!*******************************************************************************


SUBROUTINE INPUT_ELEMENT_ORDERS
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine reads a text file storing the element orders <element.ord>. 
! The format of a line is [order  0/1(=not mixed/mixed)], for every element.
!
! 2002-03-18: Created. MMB.
!*******************************************************************************

  INTEGER(I4B) :: ord_unit_num            ! unit number of the output file
  INTEGER(I4B) :: ielem,mixed_val,ios_error ! counters

  ! Initialise the unit number:
  ord_unit_num = 16

  ! Open the file:
  OPEN(UNIT=ord_unit_num,STATUS='UNKNOWN',IOSTAT=ios_error,FILE='element.ord')
  IF (ios_error.NE.0) STOP 'IE: Error opening file in INPUT_ELEMENT_ORDERS.'

  DO ielem = 1,num_elements
    READ(ord_unit_num,*) elements(ielem)%order, mixed_val
    elements(ielem)%mixed = .TRUE.
	IF (mixed_val.EQ.0) elements(ielem)%mixed = .FALSE.
  END DO

  ! Close the file:
  CLOSE(UNIT=ord_unit_num)

END SUBROUTINE INPUT_ELEMENT_ORDERS
!*******************************************************************************


SUBROUTINE OUTPUT_ELEMENT_ORDERS
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine creates a text file storing the newly established element
! orders <element.ord>. The format of a line is [order  0/1(=not mixed/mixed)].
!
! 2002-03-17: Created. MMB.
!*******************************************************************************

  INTEGER(I4B) :: ord_unit_num    ! unit number of the output file
  INTEGER(I4B) :: ielem,mixed_val ! counters

  ! Initialise the unit number:
  ord_unit_num = 16

  ! Create the new *.fek file:
  OPEN(UNIT=ord_unit_num,STATUS='REPLACE',FILE='element.ord')
  
  DO ielem = 1,num_elements
    mixed_val = 0
	IF (elements(ielem)%mixed) mixed_val = 1
	WRITE (ord_unit_num,'(2(I4))') elements(ielem)%order, mixed_val
  END DO

  ! Close the *.fek file:
  CLOSE(UNIT=ord_unit_num)

END SUBROUTINE OUTPUT_ELEMENT_ORDERS
!*******************************************************************************


SUBROUTINE OUTPUT_ERROR_DISTRIBUTION
  USE adaptive_fem
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine creates a *.fek file with the elements labels indicating 
! different error levels, as specified in <error_level_flags>.
!
! 2002-03-11: Created. MMB.
!*******************************************************************************

  INTEGER(I4B) :: fek_unit_num  ! unit number of the output file
  INTEGER(I4B) :: ielem         ! counters
  REAL(SP),DIMENSION(4,4) :: vertmat

  ! Initialise the unit number:
  fek_unit_num = 15

  ! Create the new *.fek file:
  OPEN(UNIT=fek_unit_num,STATUS='REPLACE',FILE='err_dist.fek')

  ! Start the FEM block within the fek file:
  WRITE (fek_unit_num,'(A)') 'BEGIN_FEM'

  ! Write out the data of every element, with the label equal to
  ! the values in <error_level_flags>: 
  DO ielem = 1,num_elements

    ! Find the vertices coordinates:
	CALL SIMPLEX_COEFFICIENTS(ielem,coord_mat=vertmat)

    ! Write out this element's data in the *.fek file format (see TE card):
    WRITE (fek_unit_num,'(A2,3X,5(I5),6(2X,E18.4))') 'TE', error_level_flags(ielem), 0, 0, 0, 0, &
	                                              vertmat(1,1), vertmat(2,1), vertmat(3,1), 0.0, 0.0, 0.0
    WRITE (fek_unit_num,'(5X,5(I5),6(2X,E18.4))') 0, 0, 0, 0, 0, &
	                                              vertmat(1,2), vertmat(2,2), vertmat(3,2), 0.0, 0.0, 0.0
    WRITE (fek_unit_num,'(5X,5(I5),6(2X,E18.4))') 0, 0, 0, 0, 0, &
	                                              vertmat(1,3), vertmat(2,3), vertmat(3,3), 0.0, 0.0, 0.0
    WRITE (fek_unit_num,'(5X,5(I5),6(2X,E18.4))') 0, 0, 0, 0, 0, &
	                                              vertmat(1,4), vertmat(2,4), vertmat(3,4), 0.0, 0.0, 0.0
  END DO

  ! End the FEM block within the fek file:
  WRITE (fek_unit_num,'(A)') 'END_FEM'

  ! Close the *.fek file:
  CLOSE(UNIT=fek_unit_num)

END SUBROUTINE OUTPUT_ERROR_DISTRIBUTION
!*******************************************************************************


SUBROUTINE ERM_INTER_ELEMENT_BC_MAKE(elem1,local_face1,Nvec)
  USE basis_function
  USE feminterface, ONLY: EVALUATE_ELECTRIC_FIELD
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  USE material_properties
  IMPLICIT NONE
!*******************************************************************************
! Constructs the facial contribution of an averaged, inter-element Neumann boundary
! condition in the Element Residual Method.
!
! 2002-05-09: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem1,local_face1
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE), INTENT(OUT) :: Nvec

  INTEGER(I4B) :: num_qpts                             ! number of quadrature points
  INTEGER(I4B) :: global_face                          ! global face number
  INTEGER(I4B) :: iquad,ifunc,mm,row                   ! counters
  INTEGER(I4B) :: elem2,local_face2                    ! neighbouring element
  INTEGER(I4B), DIMENSION(3) :: facenodes1,facenodes2  ! Local face nodes
  REAL(SP),DIMENSION(4,4) :: vertmat
  REAL(SP),DIMENSION(4) :: r_vec,s_vec1,s_vec2
  REAL(SP), DIMENSION(3) :: normal
  COMPLEX(SPC), DIMENSION(3) :: normal_complex,flux_av ! Error corrected DBD 13 Aug 04.
  REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: vbf_values
  COMPLEX(SPC),DIMENSION(ELEM_TET_MATRIX_SIZE) :: Nvec_large

  ! Set the number of quadrature points:
  num_qpts = 7 ! 5th order complete

  ! Establish info about the neighbour that share the face:
  global_face = elements(elem1)%faces(local_face1)
  elem2 = faceconnectelem(global_face,2)
  IF (elem2.EQ.0) STOP 'IE: invalid second element in ERM_INTER_ELEMENT_BC_MAKE.'
  CALL FIND_IN_LIST(elements(elem2)%faces(1:4),4,global_face,local_face2)
  IF (local_face2.EQ.0) STOP 'IE: invalid second local face in ERM_INTER_ELEMENT_BC_MAKE.'

  ! Establish local face nodes of elem1 and elem2:
  facenodes1 = LOCAL_FACENODES(local_face1)
  facenodes2 = LOCAL_FACENODES(local_face2)

  ! Calculate the face, normal vector, pointing away from <elem1> at <local_face1>:
  normal         = ELEMENT_NORMAL_VECTOR(elem1,local_face1)
  normal_complex = normal ! convert to complex values

  ! Establish vertices matrix for repeated use:
  CALL SIMPLEX_COEFFICIENTS(elem1,coord_mat=vertmat)
  s_vec1 = 0.0     ! vector assignment
  s_vec2 = 0.0     ! vector assignment
  Nvec_large = 0.0 ! vector assignment

  ! Integrate over the face using quadrature:
  QUAD_LOOP: DO iquad = 1,num_qpts

    ! Calculate the coords of current quad point:
    s_vec1(facenodes1) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
    s_vec2(facenodes2) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
    r_vec = MATMUL(vertmat,s_vec1)

    ! Evaluate basis functions:
    CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem1,r_vec(1:3),0,vbf_values, &
	                                  MAX_PRESENT_ORDER,MAX_PRESENT_MIXED)

    ! Cross basis functions with the normal:        
    DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	  vbf_values(ifunc,1:3) = CROSS_PRODUCT(normal,vbf_values(ifunc,1:3))
    END DO

	! Calculate the averaged $\frac{1}{\mu_r} ^n \times ^n \times \nabla \times E$
	! at the current quadrature point (<flux_av>):
	flux_av = ZERO_C ! init
	flux_av = & ! Contribution from elem1
	  flux_av + 0.5*(1.0/mu_r(elements(elem1)%material))*EVALUATE_ELECTRIC_FIELD(elem1,s_vec1,1)
	flux_av = & ! Contribution from elem2
	  flux_av + 0.5*(1.0/mu_r(elements(elem2)%material))*EVALUATE_ELECTRIC_FIELD(elem2,s_vec2,1)
    flux_av = CROSS_PRODUCT(normal_complex,flux_av) ! first  ^n \times
    flux_av = CROSS_PRODUCT(normal_complex,flux_av) ! second ^n \times
    
	! Add to total:
	DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	  Nvec_large(ifunc) = Nvec_large(ifunc) + &
	                      quad_tri_rules(num_qpts)%rule(iquad,4) * &
						  SUM(flux_av*vbf_values(ifunc,1:3))
    END DO
  END DO QUAD_LOOP

  ! Scale by face area (final step in quadrature):
  Nvec_large = FACE_AREA(elem1,local_face1) * Nvec_large

  ! Extract <Nvec> fron <Nvec_large>:
  DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row loop
    row = LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face1,mm)
    Nvec(mm) = Nvec_large(row)
  END DO

END SUBROUTINE ERM_INTER_ELEMENT_BC_MAKE
!*******************************************************************************


SUBROUTINE ERM_LOCAL_SOLVE(elem,elem_error_dof_vec)
  USE basis_function
  USE B_matrix, ONLY: B_MAKE_HIERARCHAL
  USE S_and_T_matrix, ONLY: S_AND_T_MAKE_HIERARCHAL
  USE boundary_conditions
  USE coax_feed
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TET, &
                          LOCAL_TO_GLOBAL_INDEX_TRI,         &
                          U_MAKE_HIERARCHAL ! Added DBD 31 March 2003.
  USE cbaa_sys, ONLY:CBAA_MAKE_BI_MATRIX_ELEMENTAL,     &
       MAKE_ERM_CBAA_BI_LINEAR, &
       CBAA_MAKE_COAX_MODE_INTEGRAL
  USE gw_sys, ONLY: K_Z_MN
  USE frequency_data
  USE gw_data
  USE matrix
  USE nrtype
  USE problem_info
  USE geometry
  IMPLICIT NONE
!*******************************************************************************
! This routine sets up the Element Residual Method (ERM) problem for the element 
! <elem> and solves it.
!
! 2002-05-07: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem 
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE), INTENT(OUT) :: elem_error_dof_vec

  INTEGER(I4B) :: ifunc,iface                  ! Counters
  INTEGER(I4B) :: row,row_approx,column,mm,nn, &
                  global_face,port_num         ! Indices
  INTEGER(I4B) :: elem_dof                     ! Number of dofs of the local problem
  INTEGER(I4B) :: temp_order
  LOGICAL(LGT) :: temp_mixed
  REAL(SP) :: wave_number                      ! for port contributions
  COMPLEX(SPC) :: E0                           ! complex size of the excitation
  COMPLEX(SPC) :: gamma                        ! scaling factor
  REAL(SP), DIMENSION(3) :: normal             ! Outward directed unit normal on port.
  COMPLEX(SPC) :: dof_approx                   ! dof value from last global solution
  INTEGER(I4B) :: info                         ! Info from LAPACK factorization.
  INTEGER(I4B), DIMENSION(3) :: tempedges
  LOGICAL(LGT) :: test_free  
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ipiv_error
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: A_mat_error
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: x_vec_error,b_vec_error
  INTEGER(I4B), DIMENSION(ELEM_TET_MATRIX_SIZE) :: vbf2dof_index
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE,ELEM_TET_MATRIX_SIZE) :: Se,Te ! Elemental FE matrices
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) :: Bs,CBAA_BI_mat
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: Us
  REAL(SP),     DIMENSION(ELEM_TRI_MATRIX_SIZE) :: mode_vbf_integrals
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: Nvec,BI_vec


  
  ! Create index from local elemental matrix contributions to the 
  ! actual elemental error problem (actual dofs of the elemental problem)
  ! and count the number of local element problem dofs in the process: 
  CALL ELEMENTAL_DOF_FLAGS(elem,vbf2dof_index)
  elem_dof = 0
  DO ifunc = 1,ELEM_TET_MATRIX_SIZE
    IF (vbf2dof_index(ifunc).EQ.1) THEN
	  elem_dof = elem_dof + 1
      vbf2dof_index(ifunc) = elem_dof
	END IF
  END DO   


  
  ! Allocate the matrix equation storage and initialise: 
  ALLOCATE(A_mat_error(elem_dof,elem_dof))
  ALLOCATE(x_vec_error(elem_dof))
  ALLOCATE(b_vec_error(elem_dof))
  ALLOCATE(ipiv_error(elem_dof))
  A_mat_error = (0.0,0.0)
  b_vec_error = (0.0,0.0)

  

  ! Add the straight FEM contributions ([S] & [T]) to the A-matrix:
  ! **AND**
  ! Add the straight FEM contributions ([S] & [T]) to the b-vector (I.e
  ! the bilinear form of the last approximate solution and the weighting
  ! functions of the current problem.):
  CALL S_AND_T_MAKE_HIERARCHAL(elem,MAX_PRESENT_ORDER,MAX_PRESENT_MIXED,Se,Te) ! this incorporates epsilon and mu
  COL_LOOP1: DO nn = 1,ELEM_TET_MATRIX_SIZE ! Cycle through all VBF components of the last FEM solution

	column = vbf2dof_index(nn)

    ! Store the dof value of this VBF in the last FEM solution (if it exists):
	row_approx = LOCAL_TO_GLOBAL_INDEX_TET(elem,nn)
	IF (row_approx.GT.0) dof_approx = x_vec_c(row_approx)

	ROW_LOOP1: DO mm = 1,ELEM_TET_MATRIX_SIZE

  	  row = vbf2dof_index(mm)

      IF ((row.GT.0).AND.(column.GT.0)) THEN
		! A-matrix contribution:
        A_mat_error(row,column) = A_mat_error(row,column) + Se(mm,nn) - (k0**2)*Te(mm,nn)

		IF (row_approx.GT.0) THEN
          ! B-vector contribution:
		  ! Note the minus: the bilinear form is subtracted at the linear side of the VBVP:
		  b_vec_error(row) = b_vec_error(row) - dof_approx*(Se(mm,nn)-(k0**2)*Te(mm,nn))
        END IF
      END IF
    END DO ROW_LOOP1
  END DO COL_LOOP1

  

  ! Port (coax or square-waveguide) contributions through the bilinear 
  ! forms to the A-matrix and b-vector:
  FACE_LOOP1: DO iface = 1,4

    global_face = elements(elem)%faces(iface)
	
	! Cycle if this is not a port face:
	IF (.NOT.(faces(global_face)%coax_aperture.OR.faces(global_face)%port)) CYCLE FACE_LOOP1

    ! Assign port properties:    
    IF (faces(global_face)%port) THEN ! square waveguide
      port_num    = faces(global_face)%portnumber
      wave_number = K_Z_MN(k0,port_num,mode_m,mode_n)
      gamma       = j*wave_number
      normal      = ports(port_num)%normal
	ELSE ! coaxial port
      port_num    = faces(global_face)%coaxnumber
      wave_number = k0*SQRT(AXdata(port_num)%coax_mu)*SQRT(AXdata(port_num)%coax_eps)
      gamma       = j*wave_number/AXdata(port_num)%coax_mu
	  normal      = AXdata(port_num)%normal
    END IF

    ! Calculate the elemental port face matrix:
    ! (To make B_MAKE_HIERARCHAL calculate the maximum-full matrix,
	! the elemental order is temporarily upped.)
	temp_order = faces(global_face)%order
	temp_mixed = faces(global_face)%mixed
	faces(global_face)%order = MAX_PRESENT_ORDER
	faces(global_face)%mixed = MAX_PRESENT_MIXED
	CALL B_MAKE_HIERARCHAL(elem,iface,normal,Bs)
	faces(global_face)%order = temp_order
	faces(global_face)%mixed = temp_mixed

    ! Now add <Bs> to the system matrix:
    COL_LOOP2: DO nn = 1,ELEM_TRI_MATRIX_SIZE ! column counter

      column = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,nn))

      ! Store the dof value of this VBF in the last FEM solution (if it exists):
	  row_approx = LOCAL_TO_GLOBAL_INDEX_TRI(elem,iface,nn)
	  IF (row_approx.GT.0) dof_approx = x_vec_c(row_approx)

      ROW_LOOP2: DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row counter

        row = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,mm))

        IF ((row.GT.0).AND.(column.GT.0)) THEN
		  ! A-matrix contribution:
          A_mat_error(row,column) = A_mat_error(row,column) + gamma*Bs(mm,nn)

		  IF (row_approx.GT.0) THEN
            ! B-vector contribution:
		    ! Note the minus: the bilinear form is subtracted at the linear side of the VBVP:
		    b_vec_error(row) = b_vec_error(row) - dof_approx*gamma*Bs(mm,nn)
          END IF
        END IF
      END DO ROW_LOOP2
    END DO COL_LOOP2
  END DO FACE_LOOP1



  ! CBAA contributions through the bilinear 
  ! forms to the A-matrix and b-vector:
  FACE_LOOP_CBAA: DO iface = 1,4

    global_face = elements(elem)%faces(iface)
	
	! Cycle if this is not a port face:
	IF (.NOT.faces(global_face)%CBAA_aperture) CYCLE FACE_LOOP_CBAA

    ! Calculate the elemental face's CBAA BI contribution:
    CALL CBAA_MAKE_BI_MATRIX_ELEMENTAL(elem,iface,ELEM_TRI_MATRIX_SIZE,13, &
                                       elem,iface,ELEM_TRI_MATRIX_SIZE,13, &
									   CBAA_BI_mat)

    ! Now add <CBAA_BI_mat> to the system matrix:
    COL_LOOP_CBAA: DO nn = 1,ELEM_TRI_MATRIX_SIZE ! column counter

      column = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,nn))

      ! Store the dof value of this VBF in the last FEM solution (if it exists):
	  row_approx = LOCAL_TO_GLOBAL_INDEX_TRI(elem,iface,nn)
	  IF (row_approx.GT.0) dof_approx = x_vec_c(row_approx)

      ROW_LOOP_CBAA: DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row counter

        row = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,mm))

        IF ((row.GT.0).AND.(column.GT.0)) THEN
		  ! A-matrix contribution:
          A_mat_error(row,column) = A_mat_error(row,column) + CBAA_BI_mat(mm,nn)

		  IF (row_approx.GT.0) THEN
            ! B-vector contribution:
		    ! Note the minus: the bilinear form is subtracted at the linear side of the VBVP:
		    b_vec_error(row) = b_vec_error(row) - dof_approx*CBAA_BI_mat(mm,nn)
          END IF
        END IF
      END DO ROW_LOOP_CBAA
    END DO COL_LOOP_CBAA
  END DO FACE_LOOP_CBAA



  ! CBAA BI contributions through the linear 
  ! form to the b-vector:
  FACE_LOOP_CBAA_LINEAR: DO iface = 1,4


! TEMPPPP!!!!
! CYCLE FACE_LOOP_CBAA_LINEAR


    global_face = elements(elem)%faces(iface)
	
	! Cycle if this is not a CBAA face:
	IF (.NOT.faces(global_face)%CBAA_aperture) CYCLE FACE_LOOP_CBAA_LINEAR

    ! Test whether this aperture facet has any dofs:
    tempedges = LOCAL_FACEEDGES(iface)
    test_free =                                          &
      faces(elements(elem)%faces(iface))%free.OR.        &
      edges(elements(elem)%edges(tempedges(1)))%free.OR. &
      edges(elements(elem)%edges(tempedges(2)))%free.OR. &
      edges(elements(elem)%edges(tempedges(3)))%free
       
    ! If no dofs, then no contribution to [P] possible from this element:
    IF (.NOT.test_free) CYCLE FACE_LOOP_CBAA_LINEAR

    ! Calculate the elemental face contribution:
    CALL MAKE_ERM_CBAA_BI_LINEAR(elem,iface,ELEM_TRI_MATRIX_SIZE,13,BI_vec,elem,0)

    ROW_LOOP_CBAA_LINEAR: DO mm = 1,ELEM_TRI_MATRIX_SIZE

      row = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,mm))
      IF (row.GT.0) THEN ! this is a dof
		b_vec_error(row) = b_vec_error(row) - BI_vec(mm)
      END IF
    END DO ROW_LOOP_CBAA_LINEAR
  END DO FACE_LOOP_CBAA_LINEAR





  ! Port (coax or square-waveguide) contributions through the linear 
  ! form to the b-vector:
  FACE_LOOP2: DO iface = 1,4

    global_face = elements(elem)%faces(iface)
	
	! Cycle if this is not a port face:
	IF (.NOT.(faces(global_face)%coax_aperture.OR.faces(global_face)%port)) CYCLE FACE_LOOP2

    ! Assign port properties:    
    IF (faces(global_face)%port) THEN ! square waveguide
      port_num    = faces(global_face)%portnumber
      wave_number = K_Z_MN(k0,port_num,mode_m,mode_n)
      normal      = ports(port_num)%normal
      E0          = ports(port_num)%excitation
      gamma       = (-2.0)*j*wave_number*E0
	ELSE ! coaxial port
      port_num    = faces(global_face)%coaxnumber
      wave_number = k0*SQRT(AXdata(port_num)%coax_mu)*SQRT(AXdata(port_num)%coax_eps)
	  normal      = AXdata(port_num)%normal
      E0          = AXdata(port_num)%coax_Iabs *                               &
                    EXP(j*AXdata(port_num)%coax_Iphase*PI/180.0) *             &
					(Z_zero/(2.0*PI*AXdata(port_num)%coax_a)) *                &
					SQRT(AXdata(port_num)%coax_mu/AXdata(port_num)%coax_eps)
     gamma        = E0*2.0*j*(wave_number/AXdata(port_num)%coax_mu)
    END IF

    ! Calculate the elemental port face incident b-vector:
    ! (To make B_MAKE_HIERARCHAL calculate the maximum-full matrix,
	! the elemental order is temporarily upped.)
	temp_order = faces(global_face)%order
	temp_mixed = faces(global_face)%mixed
	faces(global_face)%order = MAX_PRESENT_ORDER
	faces(global_face)%mixed = MAX_PRESENT_MIXED
    IF (faces(global_face)%port) THEN ! square waveguide
! Changed DBD 31 March 2003
	  CALL U_MAKE_HIERARCHAL(elem,iface,2,port_num,normal,Us)
! End changed DBD 31 March 2003

	ELSE ! coaxial port
	  CALL CBAA_MAKE_COAX_MODE_INTEGRAL(elem,iface,normal,           &
	                                    AXdata(port_num)%coaxcentre, &
										AXdata(port_num)%coax_a,     &
										mode_vbf_integrals)
    END IF
	faces(global_face)%order = temp_order
	faces(global_face)%mixed = temp_mixed

    ROW_LOOP3: DO mm = 1,ELEM_TRI_MATRIX_SIZE

      row = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,mm))
      IF (row.GT.0) THEN ! this is a dof
        IF (faces(global_face)%port) THEN ! square waveguide
		  b_vec_error(row) = b_vec_error(row) + gamma*Us(mm)
	    ELSE ! coaxial port
		  b_vec_error(row) = b_vec_error(row) + gamma*mode_vbf_integrals(mm)
        END IF
      END IF
    END DO ROW_LOOP3
  END DO FACE_LOOP2



  ! Add the contributions of interelement, approximate Neumann
  ! boundary conditions to the b-vector:
  FACE_LOOP_N: DO iface = 1,4

    global_face = elements(elem)%faces(iface)
	
	! Cycle if this face is not free or not internal:
	IF ((.NOT.faces(global_face)%free).OR.(faceconnectelem(global_face,2).EQ.0)) CYCLE FACE_LOOP_N

    ! Calculate the elemental interelement face matrix:
    CALL ERM_INTER_ELEMENT_BC_MAKE(elem,iface,Nvec)

    ! Now add <Nvec> to the system matrix equation:

    ROW_LOOP_N: DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row counter

      row = vbf2dof_index(LOCAL_TRI_TO_LOCAL_TET_INDEX(iface,mm))

      IF ((row.GT.0)) THEN
        ! B-vector contribution:
		b_vec_error(row) = b_vec_error(row) - Nvec(mm)
      END IF
    END DO ROW_LOOP_N
  END DO FACE_LOOP_N



  ! Solve the local problem:
  ! Uses LAPACK routine SGETRF to LU factorize (single precision)
  CALL CGETRF(elem_dof,elem_dof,A_mat_error,elem_dof,ipiv_error,info) ! LU decomposition of A_mat_c
  IF(Info.NE.0) CALL ERROR_FEMFEKO(1,4302,int1=info) ! Check error condition
  CALL CGETRS('N',elem_dof,1,A_mat_error,elem_dof,ipiv_error,b_vec_error,elem_dof,info) ! Back substitution
  IF(Info.NE.0) CALL ERROR_FEMFEKO(1,4303,int1=info) ! Check error condition
  x_vec_error = b_vec_error


  
  ! Extract the local dofs to a maximum length vector:
  elem_error_dof_vec = ZERO_C ! init
  DO ifunc = 1,ELEM_TET_MATRIX_SIZE
    IF (vbf2dof_index(ifunc).GT.0) THEN
	  row = vbf2dof_index(ifunc)
      elem_error_dof_vec(ifunc) = x_vec_error(row)
	END IF
  END DO   



  ! Deallocate the matrix equation storage: 
  DEALLOCATE(A_mat_error)
  DEALLOCATE(x_vec_error)
  DEALLOCATE(b_vec_error)
  DEALLOCATE(ipiv_error)

END SUBROUTINE ERM_LOCAL_SOLVE
!*******************************************************************************


SUBROUTINE EVALUATE_ERM_ESTIMATORS
  USE adaptive_fem
  USE feminterface, ONLY: L_SQ_NORM_ELEMENTAL
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine evaluates the L^2 norm of the estimated error field on every 
! element. The error field is estimated with the Element Residual Method.
! 
! 2002-05-09: Created. MMB.
!*******************************************************************************
  INTEGER(I4B) :: ielem,iquad,num_qpts     ! counters
  COMPLEX(SPC), DIMENSION(ELEM_TET_MATRIX_SIZE) :: elem_error_dof_vec

  ! Allocate storage for residuals and indicators:
  ALLOCATE(error_indicators(num_elements))

  ! Cycle through all elements, calculating the estimated error norm on each:
  ELEMENT_LOOP: DO ielem = 1,num_elements



! TEMP debug code:
print *, 'ERM of element:',ielem,'/',num_elements


    
	! Caculate the dofs of the approximate error on the current element:
	CALL ERM_LOCAL_SOLVE(ielem,elem_error_dof_vec)

    ! Calculate the L^2 norm of the estimated error as the indicator value:
	error_indicators(ielem) = L_SQ_NORM_ELEMENTAL(ielem,elem_error_dof_vec)

  END DO ELEMENT_LOOP

END SUBROUTINE EVALUATE_ERM_ESTIMATORS
!*******************************************************************************


SUBROUTINE DEALLOCATE_ADAPTIVE
  USE adaptive_fem
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Deallocates all datastructures associated with the module <adaptive_fem>.
!
! 2002-03-11: Created. MMB.
!*******************************************************************************

  IF (ALLOCATED(residuals_elements)) DEALLOCATE(residuals_elements)
  IF (ALLOCATED(residuals_faces))    DEALLOCATE(residuals_faces)
  IF (ALLOCATED(error_indicators))   DEALLOCATE(error_indicators)
  IF (ALLOCATED(diameters_elements)) DEALLOCATE(diameters_elements)
  IF (ALLOCATED(diameters_faces))    DEALLOCATE(diameters_faces)
  IF (ALLOCATED(error_level_flags))  DEALLOCATE(error_level_flags)
  IF (ALLOCATED(sorted_order))       DEALLOCATE(sorted_order)

END SUBROUTINE DEALLOCATE_ADAPTIVE
!*******************************************************************************


SUBROUTINE CONTROL_ADAPTIVE
  USE adaptive_fem
  USE feminterface, ONLY: ASSIGN_ERROR_LEVEL_FLAGS, DEALLOCATE_ADAPTIVE, &
                          EVALUATE_RESIDUAL_INDICATORS, OUTPUT_ERROR_DISTRIBUTION
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine calls all the other routines and cantrols all adaptive related
! tasks such as: error estinmation/indication, output of error analysis results,
! assigning of new element orders.
!
! 2002-03-11: Created. MMB.
! 2002-04-15: Updated. MMB.
!*******************************************************************************
  INTEGER(I4B) :: error_indic

  ! Exit of no error analysis is requested:
  IF (ADdata%analysis_type.EQ.0) RETURN ! no error estimator requested

  ! Evaluate the elemental error indicators/estimators:
  SELECT CASE (ADdata%analysis_type)
  CASE (0) ! no error analysis
    CONTINUE
  CASE (1) ! explicit residual analysis
    CALL EVALUATE_RESIDUAL_INDICATORS
  CASE (2) ! random assignment
    ! Place subroutine call here
  CASE (3) ! Element residual method, 0.5:0.5 weighing, L^2 norm of error.
    CALL EVALUATE_ERM_ESTIMATORS
  END SELECT

  ! Normalise the error indicator values:
  error_indicators = error_indicators - MINVAL(error_indicators)
  error_indicators = (1.0/MAXVAL(error_indicators)) * error_indicators

  ! Create an index list to the sorted order, from smallest indicated
  ! to largest: (slatec library routine)
  ALLOCATE(sorted_order(num_elements))
  STOP 'This routine must be supplied. Removed during SGI port,DBD 8 Feb 05'
!  CALL SPSORT(error_indicators,num_elements,sorted_order,1,error_indic)

  ! Assign flags to set error levels (for *.fek file visualization):
  CALL ASSIGN_ERROR_LEVEL_FLAGS

  ! Output the error levels via label assignments
  ! in a new *.fek file (for visualization):
  CALL OUTPUT_ERROR_DISTRIBUTION

  ! Assign new element orders as requested:
  CALL ASSIGN_ELEMENT_ORDERS

  ! Output the element orders if required:
  SELECT CASE (ADdata%file_orders)
  CASE (1,3)
    CALL OUTPUT_ELEMENT_ORDERS
  END SELECT

  ! Deallocate all allocatable variables in <MODULE adaptive_fem>:
  CALL DEALLOCATE_ADAPTIVE

END SUBROUTINE CONTROL_ADAPTIVE
!*******************************************************************************


