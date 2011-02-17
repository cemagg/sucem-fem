MODULE input_aux_materials
  USE nrtype
  USE problem_info
  USE material_properties
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: input_materials

  INTEGER(i4b) :: no_materials  ! Number of materials to be allocated


CONTAINS
!!!
!!! SUBROUTINE input_aux_materials(auxinfile) : Read material properties from
!!! the auxilary input file, unit auxinfile. Stores results in the data-
!!! structures in module material_properties
!!!
  SUBROUTINE input_materials(auxinfile) 
    IMPLICIT NONE
    INTEGER(i4b), INTENT(in) :: auxinfile
    INTEGER(i4b) :: ios
    CHARACTER(optionlength) :: tmp_line

    READ(auxinfile, *) no_materials 
    PRINT*, 'no_materials: ',no_materials
    CALL material_properties_allocate(no_materials) ! Allocate materials memory 
    PRINT*, 'materials allocated'
    CALL material_properties_init ! Initialise the material values to free-space

    DO                          ! Find next subblock
       READ(auxinfile, '(A)', IOstat=ios) tmp_line
       IF (ios < 0) THEN
          PRINT*, "Unexpected EOF while reading auxilary input file for materials"
          STOP
       END IF
       SELECT CASE(tmp_line) 
       CASE('BLOCK')
          PRINT*, "Unexpected BLOCK while reading auxilary input file for materials"
          STOP
       CASE('SUBBLOCK')
          CALL input_materials_subblock(auxinfile) ! Parse the sub block
       CASE('ENDBLOCK')
          EXIT
       CASE default
          CYCLE
       END SELECT
    END DO
    
  END SUBROUTINE input_materials

!!!
!!! SUBROUTINE input_aux_materials_subblock(auxinfile) : Parses a materials
!!! sub block in the auxilary input file. Assumes that the fileposition of 
!!! auxinfile is at the line following the SUBBLOCK statement. It will parse
!!! until the it finds an ENDSUBBLOCK statement.
!!!

  SUBROUTINE input_materials_subblock(auxinfile)
    USE material_properties
    IMPLICIT NONE

    INTEGER(i4b), INTENT(in) :: auxinfile
    CHARACTER(optionlength) :: tmp_line
    CHARACTER(optionlength) :: tmp_type = 'isotropic'
    REAL(sp) :: rtemp1, rtemp2
    COMPLEX(SPC) :: tmp_eps_r    ! Relative permittivity
    COMPLEX(SPC) :: tmp_eps_r_xx ! Diagonally 
    COMPLEX(SPC) :: tmp_eps_r_yy ! anisotropic 
    COMPLEX(SPC) :: tmp_eps_r_zz ! relative permittivity
    COMPLEX(SPC) :: tmp_mu_r     ! Relative permeability
    COMPLEX(SPC) :: tmp_mu_r_xx  ! Diagonally 
    COMPLEX(SPC) :: tmp_mu_r_yy  ! anisotropic 
    COMPLEX(SPC) :: tmp_mu_r_zz  ! relative permeability
    INTEGER(i4b) :: label        ! Label of the current material

    ! Initialise temp variables.
    tmp_eps_r    = 1
    tmp_eps_r_xx = 1
    tmp_eps_r_yy = 1
    tmp_eps_r_zz = 1
    tmp_mu_r     = 1
    tmp_mu_r_xx  = 1
    tmp_mu_r_yy  = 1
    tmp_mu_r_zz  = 1
    label = -1      

    DO
       READ(auxinfile, '(A)') tmp_line
       SELECT CASE(tmp_line)
          CASE('ENDSUBBLOCK')
             EXIT
          CASE('material_type')
             READ(auxinfile, '(A)') tmp_type
          CASE('material_label')
             READ(auxinfile, *) label
          CASE('eps_r')
             READ(auxinfile, *) rtemp1, rtemp2
             tmp_eps_r = CMPLX(rtemp1, rtemp2)
          CASE('eps_r_xx')
             READ(auxinfile, *) rtemp1, rtemp2
             tmp_eps_r_xx = CMPLX(rtemp1, rtemp2)
          CASE('eps_r_yy')
             READ(auxinfile, *) rtemp1, rtemp2
             tmp_eps_r_yy = CMPLX(rtemp1, rtemp2)
          CASE('eps_r_zz')
             READ(auxinfile, *) rtemp1, rtemp2
             tmp_eps_r_zz = CMPLX(rtemp1, rtemp2)
          CASE('mu_r')
             READ(auxinfile, *) rtemp1, rtemp2
             mu_r = CMPLX(rtemp1, rtemp2)
          CASE('mu_r_xx')
             READ(auxinfile, *) rtemp1, rtemp2
             mu_r_xx = CMPLX(rtemp1, rtemp2)
          CASE('mu_r_yy')
             READ(auxinfile, *) rtemp1, rtemp2
             mu_r_yy = CMPLX(rtemp1, rtemp2)
          CASE('mu_r_zz')
             READ(auxinfile, *) rtemp1, rtemp2
             mu_r_zz = CMPLX(rtemp1, rtemp2)
          CASE default
             PRINT*, "Unknown member in materials subblock found while reading axiliary input file."
             STOP
          END SELECT
       END DO

       IF ( (label < 0) .OR. (label > max_materials)) THEN
          PRINT*, "Invalid material label specified"
          STOP
       END IF
       
       SELECT CASE(tmp_type)    ! Assign material type
       CASE('PEC')
          material_type(label) = 0
          PRINT*, "PEC material specification not currently supported"
          STOP
       CASE('isotropic')
          material_type(label) = 1
       CASE('diag_anisotropic')
          material_type(label) = 2
       CASE default
          PRINT*, "Invalid material type ", TRIM(tmp_type), " specified."
       END SELECT
       
       ! Assign read values to appropriate member of the material structures
       ! from module material_properties
       eps_r(label)     = tmp_eps_r   
       eps_r_xx(label)  = tmp_eps_r_xx
       eps_r_yy(label)  = tmp_eps_r_yy
       eps_r_zz(label)  = tmp_eps_r_zz
       mu_r(label)      = tmp_mu_r    
       mu_r_xx(label)   = tmp_mu_r_xx 
       mu_r_yy(label)   = tmp_mu_r_yy 
       mu_r_zz(label)   = tmp_mu_r_zz 
       material_defined(label) = .TRUE.
  
  END SUBROUTINE input_materials_subblock
END MODULE input_aux_materials



