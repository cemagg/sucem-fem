MODULE write_matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: write_sparse_c, write_sparse_r, write_vector_c

CONTAINS
  
  SUBROUTINE write_sparse_c(cval, col_ind, row_ind)
    COMPLEX(SPC), DIMENSION(:), INTENT(in) :: cval
    INTEGER(I4B), DIMENSION(:), INTENT(in) :: col_ind, row_ind
    
    INTEGER :: i
    INTEGER :: funit

    funit = 500
    OPEN(unit = funit, file = 'output.cval', status = 'replace')
    
    DO i=1,SIZE(cval)
       WRITE (funit, '(2G20.10)') REAL(cval(i)), AIMAG(cval(i))
    END DO

    CLOSE(funit)
  END SUBROUTINE write_sparse_c

  SUBROUTINE write_sparse_r(cval, col_ind, row_ind)
    REAL(SP), DIMENSION(:), INTENT(in) :: cval
    INTEGER(I4B), DIMENSION(:), INTENT(in) :: col_ind, row_ind
    
    INTEGER :: i
    INTEGER :: funit

    funit = 500
    OPEN(unit = funit, file = 'output.cval', status = 'replace')
    
    DO i=1,SIZE(cval)
       WRITE (funit, '(1G20.10)') cval(i)
    END DO

    CLOSE(funit)
  END SUBROUTINE write_sparse_r

  SUBROUTINE write_vector_c(cvec, filename)
    COMPLEX(SPC), DIMENSION(:), INTENT(in) :: cvec
    CHARACTER(FILENAMELENGTH), INTENT(in) :: filename
    INTEGER :: i, funit
    
    funit = 500
    OPEN(unit = funit, file = filename, status = 'replace')
    
    DO i=1,SIZE(cvec)
       WRITE (funit, '(2G20.10)') REAL(cvec(i)), AIMAG(cvec(i))
    END DO

    CLOSE(funit)
  END SUBROUTINE write_vector_c

END MODULE write_matrix
  
