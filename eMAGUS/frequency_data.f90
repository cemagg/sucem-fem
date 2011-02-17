MODULE frequency_data ! This is the FEMFEKO frequency module.
   USE nrtype
   SAVE

   REAL(SP) frequency, k0, lambda0 ! Frequency, free space wavenumber and
                                   ! wavelength respectively of incident field.
   TYPE FR_card_record
     INTEGER(I4B) :: nfreq           ! Number of frequencies
     INTEGER(I4B) :: freqf           ! type of stepping: 0=linear; 1=logarithmic
     REAL(SP) :: freq0               ! Starting frequency
     REAL(SP) :: dfreq               ! Frequency delta
   END TYPE FR_card_record

   TYPE(FR_card_record)  :: FRdata

 CONTAINS

   SUBROUTINE frequency_data_defaults
     IMPLICIT NONE
     FRdata%nfreq = 0
     FRdata%freqf = 0
     FRdata%freq0 = 0
     FRdata%dfreq = 0
   END SUBROUTINE frequency_data_defaults

   SUBROUTINE frequency_error_check
     USE output_error
     IF ((FRdata%nfreq < 1).OR.(FRdata%freq0 <= 0)) CALL ERROR_FEMFEKO(1,4030)
   END SUBROUTINE frequency_error_check
   
END MODULE frequency_data
