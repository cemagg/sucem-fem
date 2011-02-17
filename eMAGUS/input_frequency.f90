!*******************************************************************************
! MODULE input_frequency_data : Subroutines and related data structures for 
!                      reading the frequency control elements from
!                      the input file namelist
!*******************************************************************************

MODULE input_frequency
  USE nrtype
  USE frequency_data
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_frequency_namelist
  
  INTEGER(I4B) :: nfreq           ! Number of frequencies
  INTEGER(I4B) :: freqf           ! type of stepping: 0=linear; 1=logarithmic
  REAL(SP) :: freq0               ! Starting frequency
  REAL(SP) :: dfreq               ! Frequency delta

  NAMELIST/FREQUENCIES/nfreq, freqf, freq0, dfreq
  
CONTAINS

!*******************************************************************************
!   SUBROUTINE read_frequency_namelist(infileunit)
!
! Reads the frequency namelist from unit infileunit, and then
! puts the data in the frequency_data  module structures
!*******************************************************************************
  SUBROUTINE read_frequency_namelist(infileunit)
    USE namelist_tools
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios
    
    CALL frequency_data_defaults ! set frequency_data module details
    CALL frequency_input_defaults ! set input routine defaults

    READ(unit=infileunit, nml=FREQUENCIES,  iostat=ios) 
    CALL check_namelist_read(ios,infileunit) 
    ! Set the module data entries to the entries as read
    FRdata%nfreq = nfreq 
    FRdata%freqf = freqf 
    FRdata%freq0 = freq0 
    FRdata%dfreq = dfreq 
    CALL frequency_error_check
  END SUBROUTINE read_frequency_namelist


  SUBROUTINE frequency_input_defaults
    nfreq = 0
    freqf = 0
    freq0 = 0
    dfreq = 0
  END SUBROUTINE frequency_input_defaults

END MODULE input_frequency
