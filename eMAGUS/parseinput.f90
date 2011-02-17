MODULE parseinput
  USE problem_info
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_namelists
  PUBLIC :: meshfilename, outfilename, fekoout_filename, run_label, meshout_filename
  PUBLIC :: mesh_type, extra_meshfilename
  SAVE

  ! NAMELIST/INPUTFILE/ 
  CHARACTER(FILENAMELENGTH)  :: meshfilename, aux_inputfilename, extra_meshfilename
  CHARACTER(RUN_LABEL_LINE_LENGTH), DIMENSION(RUN_LABEL_NUM_LINES)  :: run_label
  CHARACTER(optionlength) :: mesh_type
  NAMELIST/INPUTFILE/ meshfilename, aux_inputfilename, run_label, mesh_type, &
       extra_meshfilename

  ! NAMELIST/OUTPUTFILE/ 
  CHARACTER(FILENAMELENGTH)  :: outfilename
  CHARACTER(FILENAMELENGTH)  :: fekoout_filename
  CHARACTER(FILENAMELENGTH)  :: meshout_filename
  NAMELIST/OUTPUTFILE/ outfilename, fekoout_filename, meshout_filename

    

CONTAINS

  SUBROUTINE read_namelists(infileunit)
    USE input_problem_info
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios
    
    CALL general_defaults       ! Set up default values for namelist variables.
    CALL read_general_namelists(infileunit)
  END SUBROUTINE read_namelists
    

  SUBROUTINE read_general_namelists(infileunit)
    USE namelist_tools
    USE input_problem_info
    USE unit_numbers
    USE input_aux_inputfile
    USE problem_info, ONLY: WRITE_MESH
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: infileunit
    INTEGER(I4B) :: ios

    READ(unit=infileunit, nml=outputfile, iostat=ios)
    CALL check_namelist_read(ios,infileunit)
    ! Open the output file first, so that errors can be reported.
    OPEN (unit=fileout,status='UNKNOWN',file=outfilename)
    
    READ(unit=infileunit, nml=inputfile,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)

    READ(unit=infileunit, nml=debugging,  iostat=ios)
    CALL check_namelist_read(ios,infileunit)
    CALL read_problem_namelists(infileunit) ! Reads all other namelists.
    rewind(infileunit)

    IF (.NOT. WRITE_MESH) THEN
       ! Open aux inputfile for reading
       OPEN(unit=auxinfile, status='OLD', file=aux_inputfilename)
       CALL read_aux_inputfile(auxinfile) ! Parses the aux inputfile
    END IF
    
  END SUBROUTINE read_general_namelists
  

  SUBROUTINE general_defaults
!*********************************************************************
! Setup defaults for the namelist entries.
!*********************************************************************
    IMPLICIT NONE
    !   NAMELIST/INPUTFILE/
    meshfilename = 'mesh.fek'
    aux_inputfilename = 'auxfile.dat'
    run_label = ' '  
    mesh_type = 'oldfek'
    !   NAMELIST/DEBUGGING/
    DEBUG_AREA                    = .FALSE.
    DEBUG_EIGENVECTORS            = .FALSE.
    OUTPUT_ELEMENT_DATA           = .TRUE.
    DEBUG_ELEMENT_E_FIELD         = .FALSE.
    OUTPUT_ELEMENT_SHAPE          = .FALSE.
    DEBUG_EDGEMAKE_VECTORS        = .FALSE. 
    DEBUG_EDGEMAKE_EDGES          = .FALSE. 
    DEBUG_FACEMAKE_FACES          = .FALSE. 
    DEBUG_GW_PRE                  = .FALSE.
    DEBUG_VOLUME                  = .FALSE. 
    DEBUG_FACEAREA                = .FALSE. 
    DEBUG_EDGE_CONNECT            = .FALSE.  
    DEBUG_FACE_CONNECT            = .FALSE.  
    DEBUG_NEWTON_RAPHSON          = .FALSE.    
    DEBUG_MESHIN                  = .FALSE.  
    DEBUG_NUMBER_DOF              = .FALSE.
    DEBUG_NODAL_E_FIELDS          = .FALSE.
    DEBUG_SKYLINE                 = .FALSE.
    DEBUG_SPARSE_MATVEC           = .FALSE.
    DEBUG_SYSTEM_ELEMENTS         = .FALSE.
    DEBUG_SYSTEM_MATRIX           = .FALSE.
    DEBUG_SOLVE_SPARSE            = .FALSE.
    DEBUG_TD                      = .FALSE.
    DEBUG_TMAKE                   = .FALSE. 
    DEBUG_WIJ                     = .FALSE.  
    DEBUG_WRITE_MATRIX            = .FALSE.

    !  NAMELIST/OUTPUTFILE/ outfilename, fekoout_filename, meshout_filename
    outfilename = 'femfeko.out'
    fekoout_filename = 'output.efe'
    meshout_filename = 'output.femmesh'

  END SUBROUTINE general_defaults

END MODULE parseinput
