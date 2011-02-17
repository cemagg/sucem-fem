MODULE namelist_tools
  USE nrtype
  
CONTAINS

  SUBROUTINE check_namelist_read(ios,fileunit)
    USE output_error
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: ios, fileunit

    IF (ios > 0) THEN           ! This means some real file error has occured
       CALL error_femfeko(2, 5001)

    ELSEIF (ios < 0) then       ! Probably the namedlist is not contained in the
       REWIND(fileunit)         ! input file. Rewind so that next namelist read
    END IF                      ! doesn't fail.
    REWIND(fileunit)
  END SUBROUTINE check_namelist_read

END MODULE namelist_tools
