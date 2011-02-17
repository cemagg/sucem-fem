!*******************************************************************************
! MODULE input_aux_inputfile: 
!
! Parses the auxilary input file used for multi-issue control cards.
!                             
!*******************************************************************************

MODULE input_aux_inputfile
  USE nrtype
  USE problem_info
  USE blockinput_tools
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_aux_inputfile

CONTAINS
  SUBROUTINE read_aux_inputfile(auxinfile)
    USE input_aux_materials
    IMPLICIT NONE
    INTEGER(i4b), INTENT(in) :: auxinfile
    CHARACTER(optionlength) :: block_type
    INTEGER(i4b) :: ios
    
    DO 
       block_type =  blockinput_nextblock(auxinfile, ios)
       IF (ios < 0) EXIT        ! End of file, presumably
       SELECT CASE(block_type)
       CASE('materials')        ! Material properties
          CALL input_materials(auxinfile)
       CASE default
          PRINT*, 'Unknown BLOCK type: ', TRIM(block_type), ' encountered'
          STOP
       END SELECT
    END DO
  END SUBROUTINE read_aux_inputfile

END MODULE input_aux_inputfile
