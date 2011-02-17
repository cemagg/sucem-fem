MODULE output_error
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Internal routines:
! (A routine is either listed as an interface, or in the private statement.)


PRIVATE
PUBLIC :: error_femfeko

CONTAINS

!*******************************************************************************

! NB - documentation in SUBROUTINE DESCRIPTION needs to be updated, 
! usage has changed since Matthys' mods.
! Also - errors have changed, input errors now extend beyond 4049. 
! 
! Last changed:
! 27 February 2002 DBD - error 4106 added.
! 25 April 2001 DBD - error 4082 added.
!  3 Apr 2001 DBD - errors 4500-4502 description slightly.
! 27 March 2001 DBD - error 4509 removed (after adding). 
! 20 March 2001 DBD - errors 4080 and 4081 added. 
! 19 March 2001 DBD - error 4902 added. 
! 16 March 2001 DBD - error 4037 changed. 
! 17 Feb 2001 DBD - additional errors 4058-60 added.


SUBROUTINE ERROR_FEMFEKO(err_type, err_number, int1, int2, int3)
  USE unit_numbers
  IMPLICIT NONE   
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine writes out error messages arising from within FEMFEKO.
! Note that due to the optional arguments, routines calling it must
! USE the relevant INTERFACE block in feminterface:
! eg:   USE feminterface, ONLY: ERROR_FEMFEKO, .... other routines....
! 
!*******************************************************************************
! AUTHORS
!*******************************************************************************
! DB Davidson,  RH Geschke, MM Botha
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! Original version Jan 2000 by DBD
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Error type, error number and explanatory message.
!
! Note that this routine uses the OPTIONAL feature of FORTRAN 90. Only the 
! first 2 arguments are required, viz. err_type, err_number. 
! As noted above, to use this routine, a 
!   USE feminterface, ONLY: ERROR_FEMFEKO, .... other routines....
! statement is REQUIRED in routines calling ERROR_FEMFEKO. 
! 
! The other optional arguments may be included, so that one can print an integer(s) 
! and/ or additional lines of error reporting, and there is provision for possible
! German translations at some future stage. 
!
! USAGE: Basic usage is CALL ERRROR_FEMFEKO(1,NNNN,
!        These two arguments are REQUIRED. Use 2 for an error that occured before
!        the output file was opened and has to be echoed to screen, then stopping,
!        1 for an error (this stops the code), and 0 for a warning (the code 
!        continues running). NNNN is a four-digit error number - see the table below.
!        Please include an ERROR DIRECTORY in your routines - see example in 
!        meshin.f90
!
!        If the error is an internal error, then precede text with:
!        'IE: Internal error: ...further explanation here...'
!        and add a German line, with the text 'IF'.
!
! ADVANCED USAGE: As an example, here is an error message for an internal error
!        with an integer  that is incorrect (wrong_int),
!        2nd and 3rd line of English text and a German line:
!
!         CALL ERRROR_FEMFEKO(1,4999,&
!         'IE:Internal error, this integer caused a failure: ',int1=wrong_int,&
!         english_line2='in routine SOME_ROUTINE, please',&
!         english_line3='report fault to authors.',deutsch_line1='IF'). 
!
! FEMFEKO error numbers are in the range 4000...4999 (reserved
! within FEKO for FEMFEKO errors.)
! Within FEMFEKO, error numbers have the following ranges:
!   4000...4049: Data inputting section errors.
!   4050...4059: Elemental matrix computation errors.
!   4060...4069: Connection routine errors. 
!   4070...4079: Renumbering routine errors.
!   4080...4099: Reserved for future pre-processing errors.
!
! 4000...4099: Inputting and pre-processing errors.
! 4100...4199: General errors in analysis sections, utility routines &
!              functions.
! 4200...4299: Eigenanalysis section errors (real & complex)
! 4300...4399: Cavity-backed aperture analysis section errors.
! 4400...4499: FSS (periodic) analysis section errors.
! 4500...4599: Wavguide analysis section errors.
! 4600...4699: FEM/BEM analysis section errors.
! 4700...4799: Time domain analysis section errors.
! 4800...4899: Reserved for future analysis section. 
! 4900...4999: Post-processing errors.
!   4900...4949: Main post-processing routines.
!   4950...4959: Simplex coordinate routines.
!   4960...4969: Element searching routines.
!   4970...4979: Field computation routines/functions.
! 5000...5099: More input file errors (NM 2005.03.23)
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is an explanatory message w.r.t. the error, printed in the 
! standard FEMFEKO output file. 
!
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! None
!
!*******************************************************************************

  INTEGER(I4B), INTENT(IN) :: err_type       ! 1: Error (code stops)
                                             ! 0: Warning (code continues)
  INTEGER(I4B), INTENT(IN) :: err_number     ! See above.
  

  INTEGER(I4B), INTENT(IN), OPTIONAL  :: int1,int2,int3      
                                             ! Optional integers for error info.
  

  ! Write the error type out:
  SELECT CASE(err_type)
  CASE(0) ! Warning
     WRITE(FILEOUT,'(/,1X,A)',ADVANCE='NO')   'WARNING:'
     WRITE(FILEOUT,'(1X,I5,A)',ADVANCE='NO')  err_number,':'
  CASE(1) ! Error
     WRITE(FILEOUT,'(/,1X,A)',ADVANCE='NO')   'ERROR:'
     WRITE(FILEOUT,'(1X,I5,A)',ADVANCE='NO')  err_number,':'
  CASE(2) ! Output file not open yet
     WRITE(*,'(/,1X,A)',ADVANCE='NO')   'ERROR:'
     WRITE(*,'(1X,I5,A)',ADVANCE='NO')  err_number,':'
  CASE DEFAULT ! Shouldn't happen, but for safety:
    WRITE(FILEOUT,'(/,1X,A)',ADVANCE='NO')   'EXCEPTION:'
    WRITE(FILEOUT,'(1X,I5,A)',ADVANCE='NO')  err_number,':'
  END SELECT
  
  ! Write the massage out:
  SELECT CASE(err_number)


  ! Inputting and pre-processing errors: *****************************************
  CASE (4000)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Invalid excitation type (',int1,') in A0 card number ',int2,'.'
  CASE (4001)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'A0 card number ',int1,' was forced to be a new exitation.'
  CASE (4002)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'A0 card number ',int1,' cannot have multiple angles of incidence as well.'
  CASE (4003)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Invalid excitation type (',int1,') in A8 card number ',int2,'.'
  CASE (4004)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Invalid probe direction (',int1,') in A8 card number ',int2,'.'
  CASE (4005)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'A8 card number ',int1,' was forced to be a new exitation.'
  CASE (4007)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid type of boundary in BC card number ',int1,'.'
  CASE (4008)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid type of ground plane (',int1,') in BO card.'
  CASE (4009)
    WRITE (FILEOUT,'(1X,A)') 'Invalid ARPACK spectrum control option. Possibilities are SM, LM, BE, SA or LA.'
  CASE (4010)
    WRITE (FILEOUT,'(1X,A)') 'Upper limit on number of ARPACK iterations too small. Increased to default of 500.'
  CASE (4011)
    WRITE (FILEOUT,'(1X,A)') 'ARPACK modes from 2 to 5 available.'
  CASE (4012)
    WRITE (FILEOUT,'(1X,A)') 'Eigenvalue spectrum shift must be zero of positive.'
  CASE (4013)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Unavailable field averaging option (',int1,') in FA card.'
  CASE (4014)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Unavailable hierarchal order (',int1,') in FC card.'
  CASE (4015)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid analysis type (',int1,') in FC card.'
  CASE (4016)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Unimplemented FE type (',int1,') in FE card number ',int2,'.'
  CASE (4017)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Unimplemented FE co-ordinate system (',int1,') in FE card number ',int2,'.'
  CASE (4018)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'No corresponding excitation for FE card number ',int1,'.'
  CASE (4019)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid parameter FFREQ specified in FF card number ',int1,'.'
  CASE (4020)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid parameter RIGE specified in FF card number ',int1,'.'
  CASE (4021)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'No corresponding excitation for FF card number ',int1,'.'
  CASE (4022)
    WRITE (FILEOUT,'(1X,A)') 'Invalid matrix storage format requested in FM card.'
  CASE (4023)
    WRITE (FILEOUT,'(1X,A)') 'Invalid type of matrix solver requested in FM card.'
  CASE (4024)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid type of frequency stepping requested by FR card number ',int1,'.'
  CASE (4025)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Negative material label not allowed in MA card number ',int1,'.'
  CASE (4026)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Invalid material type in MA card number ',int1,'.'
  CASE (4027)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Maximum material referenced in geometry (',int1,') not defined.'
  CASE (4028)
    WRITE (FILEOUT,'(1X,A)') 'Sparse eigenanalysis no longer supported.'
  CASE (4029)
    WRITE (FILEOUT,'(1X,A)') 'Sparse, complex eigenanalysis not yet implemented.'
  CASE (4030)
    WRITE (FILEOUT,'(1X,A)') 'At least 1 positive frequency is required for frequency domain analysis.'
  CASE (4031)
    WRITE (FILEOUT,'(1X,A)') 'No ground plane (BO card) defined for CBAA.'
  CASE (4032)
    WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRF in subroutine FACE_IN_QUADRILATERAL.'
  CASE (4033)
    WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRI in subroutine FACE_IN_QUADRILATERAL.'
  CASE (4034)
    WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRF in subroutine SIMPLEX_COEFFICIENTS.'
  CASE (4035)
    WRITE (FILEOUT,'(1X,A)') 'IE: Error with SGETRI in subroutine SIMPLEX_COEFFICIENTS.'
  CASE (4036)
    WRITE (FILEOUT,'(1X,A)') 'Guided wave analysis requires at least one port.'
  CASE (4037)
    WRITE (FILEOUT,'(1X,A,I2)') 'Unimplemented GW analysis type requested in GW card: ',int1
  CASE (4038)
    WRITE (FILEOUT,'(1X,A)') 'IE: Invalid hierarchal order in routine FIND_FREEDOM'
  CASE (4039)
    WRITE (FILEOUT,'(1X,A)') 'No degrees of freedom, likely cause of error: invalid mesh.'
  CASE (4040)
    WRITE (FILEOUT,'(1X,A)') 'IE: Internal error in NUMBER_DOF. Error during renumbering. This should not happen! Contact authors.'
  CASE (4041)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Invalid material type (',int1,') for material no. ',int2,'.'
  CASE (4042)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'Real eigenanalysis code cannot handle lossy material no. ',int1,'.'
  CASE (4043)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Eigenanalysis code cannot handle material no. ',int1,' of type ',int2,'.'
  CASE (4044)
    WRITE (FILEOUT,'(1X,A)') 'FF card(s) not allowed with eigenanalysis.'
  CASE (4045)
    WRITE (FILEOUT,'(1X,A)') 'A8 card(s) not allowed with eigenanalysis.'
  CASE (4046)
    WRITE (FILEOUT,'(1X,A)') 'A0 card(s) not allowed with eigenanalysis.'
!!!
!!! FR cards have been removed, and the new input code won't try to add frequency
!!! info for eigen analyses NM
!!!
!!!  CASE (4047)
!!!    WRITE (FILEOUT,'(1X,A)') 'FR card(s) not allowed with eigenanalysis.'
!!! 
  CASE (4048)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Invalid excitation type (',int1,') in A9 card number ',int2,'.'
  CASE (4049)
    WRITE (FILEOUT,'(1X,A,I4,A)') 'A9 card number ',int1,' was forced to be a new exitation.'
  CASE (4050)
    WRITE (FILEOUT,'(1X,A)') 'Number of ports must be equal to number of A9 cards'
  CASE (4051)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Port label (',int2,'), defined in BC card no. ',int1,' is out of range.'
  CASE (4052)
    WRITE (FILEOUT,'(1X,A)') 'Multiple ports are assigned the same label(s).'
  CASE (4053)
    WRITE (FILEOUT,'(1X,A,I4,A,I4,A)') 'Port label (',int2,'), defined in A9 card no. ',int1,' is out of range.'
  CASE (4054)
    WRITE (FILEOUT,'(1X,A,I4)',ADVANCE='NO') 'Excitation for port no. ',int1
    WRITE (FILEOUT,'(A,I4,A)') ' is assigned a socond time in the same analysis by A9 card no. ',int2,'.'
  CASE (4055)
    WRITE (FILEOUT,'(1X,A)') 'A0 card(s) not allowed in GW analysis.'
  CASE (4056)
    WRITE (FILEOUT,'(1X,A)') 'A8 card(s) not allowed in GW analysis.'
  CASE (4057)
    WRITE (FILEOUT,'(1X,A)') 'FF card(s) not allowed in GW analysis.'
  CASE (4058)
    WRITE (FILEOUT,'(1X,A)') 'Invalid elemental matrix computation parameter in FX card'
  CASE (4059)
    WRITE (FILEOUT,'(1X,A)') 'Edge scaling may not be used with A&V-type elements'
  CASE (4060)
    WRITE (FILEOUT,'(1X,A)') 'Only cubature implemented for higher-order A&V elements'
  CASE (4061)  
    WRITE (FILEOUT,'(1X,A,I4)') 'Center direction must be = 1/2/3 at coax port: ',int1
  CASE (4062)
    WRITE (FILEOUT,'(1X,A,I4)') 'Center length cannot be equal to 0 at coax port: ',int1
  ! Additional pre-processing errors
  CASE (4080)
    WRITE (FILEOUT,'(1X,A,I4)') 'Normal incorrectly specified at port: ',int1
  CASE (4081)
    WRITE (FILEOUT,'(1X,A,I4)') 'Tangent incorrectly specified at port: ',int1
  CASE (4082)
    WRITE (FILEOUT,'(1X,A,I4)') 'Preconditioner not implemented for this solver: ',int1
! Added DBD March-June 03
  CASE (4083)
    WRITE (FILEOUT,'(1X,A)') 'Only one TD card permitted.'
  CASE (4084)
    WRITE (FILEOUT,'(1X,A)') 'Unimplemented time domain option (TD card).'
  CASE (4085)
    WRITE (FILEOUT,'(1X,A,I4)') 'Normal defined is not a unit normal on ABC: ',int1
  CASE (4086)
    WRITE (FILEOUT,'(1X,A,I4)') 'More than one boundary condition has been defined as an ABC with this number: ',int1
  CASE (4087)
    WRITE (FILEOUT,'(1X,A,I4)') 'No boundary condition has been defined as an ABC with this number: ',int1
  CASE (4088)
    WRITE (FILEOUT,'(1X,A)') 'Only one QU card permitted.'
  CASE (4089)
    WRITE (FILEOUT,'(1X,A)') 'A QU card is required for TD analysis.'
  CASE (4090)
    WRITE (FILEOUT,'(1X,A)') 'Six ABC''s are required, numbered 1 through 6.'
  CASE (4091)
    WRITE (FILEOUT,'(1X,A)') 'ABC bounding box point S1 or S2 not defined.'
  CASE (4092)
    WRITE (FILEOUT,'(1X,A)') 'For ABC bounding box, S1 must be closer to origin than S2.'
  CASE (4093)
    WRITE (FILEOUT,'(1X,A)') 'One (and only one) A0 card required for time-domain analysis.'
  CASE (4094)
    WRITE (FILEOUT,'(1X,A)') 'Non-zero phase specified in A0 card for time-domain analysis.'
  CASE (4095)
    WRITE (FILEOUT,'(1X,A,I6)') 'Unsupported type of ABC, ABC number ',int1
  CASE (4096)
    WRITE (FILEOUT,'(1X,A)') 'Only one KU card permitted.'
  CASE (4097)
    WRITE (FILEOUT,'(1X,A)') 'Only one PL card permitted.'
  CASE (4098)
    WRITE (FILEOUT,'(1X,A)') 'Error in definition of spherical boundary for scattered field treatment.'
  CASE (4099)
    WRITE (FILEOUT,'(1X,A)') 'PML can only be used with scattered field formulation.'
! See 5000...5099 for more input errors.
! End DBD additions March-June 03
  ! General errors in analysis sections, utility routines functions: *************
  CASE (4100)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'IE: In function SIMPLEX_COORDINATES following SGETRF, info=',int1,'.'
  CASE (4101)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'IE: In function SIMPLEX_COORDINATES following SGETRI, info=',int1,'.'
  CASE (4102)
    WRITE (FILEOUT,'(1X,A)') 'IE: In function SIMPLEX_COORDINATES, simplex coordinates do not sum to 1.'
  CASE (4103)
    WRITE (FILEOUT,'(1X,A)') 'IE: Invalid simplex sum in XYZ_COORDINATES.'
  CASE (4104)
    WRITE (FILEOUT,'(1X,A)') 'IE: Sparse matrix entries not symmetrical.'
  CASE (4105)
    WRITE (FILEOUT,'(1X,A)') 'IE: Necessary variables not allocated in ITER_SOLVE.'
  CASE (4106)
    WRITE (FILEOUT,'(1X,A)') 'IE: Edge length scaling not supported with quadrature'
  CASE (4107)
    WRITE (FILEOUT,'(1X,A)') 'IE: Necessary variables not allocated in USER_DITSOL_PCG.'
  CASE (4108)
    WRITE (FILEOUT,'(1X,A,I8)') 'Iterative solver did not converge on call ',int1
  ! Eigenanalysis section errors (real & complex): *******************************
  CASE (4200)
    WRITE (FILEOUT,'(1X,A)') 'Banded storage not implemented in EIGEN_SYSMAT for complex eigenanalysis.'
  CASE (4201)
    WRITE (FILEOUT,'(1X,A)') 'Error in routine EIGEN_SYSMAT, complex eigenmodes not yet implemented.'
  CASE (4202)
    WRITE (FILEOUT,'(1X,A)') 'Number of Arnoldi vectors must be more than number of eigenvalues requested.'
  CASE (4203)
    WRITE (FILEOUT,'(1X,A,I6)',ADVANCE='NO') 'SSYGV failed to converge.',int1
    WRITE (FILEOUT,'(A)') ' off-diagonal elements of an intermediate tridiagonal form did not converge to zero.'
  CASE (4204)
    WRITE (FILEOUT,'(1X,A)') 'SSBDR: N is greater than MAXN.'
  CASE (4205)
    WRITE (FILEOUT,'(1X,A)') 'SSBDR: NEV is greater than MAXNEV.'
  CASE (4206)
    WRITE (FILEOUT,'(1X,A)') 'SSBDR: NCV is greater than MAXNCV.'
  CASE (4207)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'IE: _sband, info=',int1,'. Check the documentation of _sband.'
  CASE (4208)
    WRITE (FILEOUT,'(1X,A,I6)',ADVANCE='NO') 'SSYGV: Leading minor of order',int1
    WRITE (FILEOUT,'(A)',ADVANCE='NO') ' is not positive definite. Factorization of B could not be completed'
    WRITE (FILEOUT,'(A)') ' and no eigenvalues or eigenvectors were computed.'


  ! Cavity-backed aperture analysis section errors: ******************************
  CASE (4300)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'IE: info=',int1,', SGETRF in CBAA_SELF_TERM routine.'
  CASE (4301)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'IE: info=',int1,', SGETRI in CBAA_SELF_TERM routine.'
  CASE (4302)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'info=',int1,', of CGETRF in DIRECT_SOLVE routine.'
  CASE (4303)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'info=',int1,', of CGETRS in DIRECT_SOLVE routine.'


  ! Waveguide analysis section errors: ********************************************
  CASE (4500)
    WRITE (FILEOUT,'(1X,A,A)') 'Waveguide a and b dimensions identical,',&
                               'degenerate modes may exist'
  CASE (4501)
    WRITE (FILEOUT,'(1X,A)') 'Waveguide dimensions at different ports not identical.'
  CASE (4502)
    WRITE (FILEOUT,'(1X,A)') 'Dominant waveguide mode is under cut-off.'
! DBD change 25 March 2003 and 17 Dec 03
  CASE (4503)
    WRITE (FILEOUT,'(1X,A)') 'Banded storage not implemented for GW_SYSMAT or TD analysis or FD scattering analysis.'
! End DBD change 25 March 2003 and 17 Dec 03
  CASE (4504)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'More than one face of element no. ',int1,' on a port.'
  CASE (4505)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'Argument number ',int1,' of call to CSYTRF was invalid in GW_SYS.'
  CASE (4506)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'On exit of CSYTRF. Diagonal element ',int1,' is exactly zero and matrix is singular in GW_SYS.'
  CASE (4507)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'Argument number ',int1,' of call to CSYTRS was invalid in GW_SYS.'
  CASE (4508)
    WRITE (FILEOUT,'(1X,A,I6)',ADVANCE='NO') 'During assembly by elements at element no. ',int1
    WRITE (FILEOUT,'(A)') ' more than one port was defined for this element.'
  CASE (4509)
    WRITE (FILEOUT,'(1X,A)') 'Waveguide can support multi-mode propagation.'
  CASE (4510)
    WRITE (FILEOUT,'(1X,A,I2,A)') 'IE: port corner sorting failed in WG_DIMS, at port:',&
    int1, ' Possible workaround: change order of corners.'
  CASE(4511) 
    WRITE (FILEOUT,'(1X,A,I2,A)') 'Error in definition of port: ',int1,&
    '. No element was found with a face aligned with this port.'
  CASE(4512)
    WRITE (FILEOUT,'(1X,A,I2,A)') 'Port ',int1, &
    ' shares an edge with another port.'
  ! Time domain analysis section errors: ********************************************
  CASE(4700)  
    WRITE (FILEOUT,'(1X,A,I2,A)') 'Error in definition of ABC, label number: ',int1,&
    '. No element was found with an external face aligned with this ABC.'
  CASE(4701)  
    WRITE (FILEOUT,'(1X,A,I6,A)') 'Argument number ',int1,' of call to SGETRF was invalid in TD_SYSMAT.'
  CASE(4702)  
    WRITE (FILEOUT,'(1X,A,I6,A)') 'Factor U(i,i) exactly singular, i=',int1,&
	  'solution cannot be computed. (Call to SGETRF in TD_SYSMAT).'
  CASE(4703)  
    WRITE (FILEOUT,'(1X,A,I6,A)') 'Argument number ',int1,' of call to SGETRS was invalid in TD_TIMESTEPPING.'
  CASE (4704)
    WRITE (FILEOUT,'(1X,A,I6,A)') 'More than one face of element no. ',int1,' on the same ABC.'
  CASE(4705)  
    WRITE (FILEOUT,'(1X,A,I6)') 'Sparse factorization routine DSSKYF exited abnormally in TD_SYSMAT, error =',int1
  CASE(4706)  
    WRITE (FILEOUT,'(1X,A,I6)') 'Sparse solver routine DSSKYS exited abnormally in TD_TIMESTEP, error =',int1
  ! Post-processing errors: ******************************************************
  CASE (4900)
    WRITE (FILEOUT,'(1X,A)') 'No FF request possible for this analysis type.'
  CASE (4901)
    WRITE (FILEOUT,'(1X,A)') 'FF 2 card with no incident direction specified.'    
  CASE (4902)
    WRITE (FILEOUT,'(1X,A)') 'Near fields not computed for S-parameter analysis.'
  CASE (4903)
    WRITE (FILEOUT,'(1X,A)') 'Scattered fields not available with this type of analysis.'
 CASE (5000)
    WRITE (FILEOUT,'(1X,A)') 'ARPACK sparse eigenanalysis no longer supported.'
 CASE (5001)
    WRITE (FILEOUT,'(1X,A)') 'Error reading input.dat. Possibly incorrect variable specification'
    WRITE (*,'(1X,A)') 'Error reading input.dat. Possibly incorrect variable specification'
  CASE DEFAULT
    WRITE(FILEOUT,'(1X,A)') 'IE: Not a valid FEMFEKO error number.'
  END SELECT

  ! STOP in the case of an error:
  SELECT CASE(err_type)
  CASE(1)
    STOP ! Stop on error
  CASE(2)
    STOP ! Stop on error
  CASE DEFAULT 
    CONTINUE 
  END SELECT

  RETURN 
END SUBROUTINE ERROR_FEMFEKO
!*******************************************************************************


END MODULE output_error

