!!!
!!! MODULE blockinput_tools:
!!!
!!! Contains some tools to easy the parsing of femfeko block-format input files
!!!

MODULE blockinput_tools
  USE nrtype
  USE problem_info
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: blockinput_nextblock

CONTAINS

!!!
!!! FUNCTION blockinput_nextblock(blockfile,ios)
!!!
!!! Finds the beginning of the next block in an femfeko block-format input
!!! file. 
!!!
!!! input: Unit number of open block-format input file
!!!
!!! output: 
!!!
!!! Puts the current file position at the line after the BLOCK statement
!!! and returns a string containing the entry type (i.e. the string after the
!!! BLOCK command in the input file.
!!! 
!!! Returns the file read status in ios. If EOF was reached, ios will be < 0
!!!
  FUNCTION blockinput_nextblock(blockfile, ios)
    IMPLICIT NONE
!!!
!!! Interface variables
!!!
    INTEGER(i4b), INTENT(in) :: blockfile
    INTEGER(i4b), INTENT(out) :: ios
    CHARACTER(optionlength) :: blockinput_nextblock
!!!
!!! Internal variables
!!!
    CHARACTER(6) :: temp

    DO 
       READ(blockfile, '(A6, A)', IOSTAT=ios) temp, blockinput_nextblock 
       SELECT CASE(ios) 
       CASE(:-1) ! < 0, End of file, presumably
          EXIT         
       CASE(1:) ! > 0, Serious error
          PRINT*, 'Error reading from block file with unit: ', blockfile
          STOP
       END SELECT
       
       IF (temp .NE. 'BLOCK ') THEN
          CYCLE                 ! Does not seem to be the beginning of a BLOCK
       ELSE
          EXIT                  ! Block found!
       END IF
    END DO
    
  END FUNCTION blockinput_nextblock
  
END MODULE blockinput_tools

  
