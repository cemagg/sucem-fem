MODULE output_mesh
  USE nrtype
  USE unit_numbers
  use geometry
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: output_mesh_nodes

CONTAINS

  SUBROUTINE output_mesh_nodes(mesh_unit)
    IMPLICIT NONE
!!!
!!! interface variables
!!!
    INTEGER(i4b), INTENT(in) :: mesh_unit
!!!
!!! Internal variables
!!!
    INTEGER(i4b) :: i

    WRITE(mesh_unit, '(A)') 'BLOCK nodes'
    WRITE(mesh_unit, '(I0)') num_nodes
    DO i = 1, num_nodes
       WRITE(mesh_unit, '(I0, G0, G0, G0)') i, vertices(i)%coord
    END DO
    WRITE(mesh_unit, '(A)') 'ENDBLOCK'

    WRITE(mesh_unit, '(A)') 'BLOCK tets'
    WRITE(mesh_unit, '(I0)') num_elements
    DO i = 1, num_elements
       WRITE(mesh_unit, '(6I0)') i, elements(i)%material, elements(i)%nodes
    END DO
    WRITE(mesh_unit, '(A)') 'ENDBLOCK'

  END SUBROUTINE output_mesh_nodes
  

END MODULE output_mesh
