MODULE output_data
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: output_materials_info, output_mesh_info


CONTAINS 
  
!!!
!!! SUBROUTINE output_materials_info:
!!!
!!! Writes information about the material parameters to the output file.
!!!
  SUBROUTINE output_materials_info
    USE material_properties
    USE output_error
    
    INTEGER(i4b) :: i_mat
    CHARACTER(12) :: status  ! string to write the material definition status out

    ! Write material information to the output file:
    WRITE(FILEOUT,'(//,20X,A,/)') 'CONSTITUTIVE PARAMETERS'
    WRITE(FILEOUT,'(3(1X,A12),2(1X,A25))') 'Material no.','Status','Type','eps_r       ','mu_r       '
    PRINT*, 'max_materials: ', max_materials
    DO i_mat = 0,MAX_MATERIALS 
       ! Assign the status string:
       IF (material_defined(i_mat)) THEN
          status = '     defined'
       ELSE
          status = '   undefined'
       END IF
       ! write this material's information:
       SELECT CASE (material_type(i_mat))
       CASE (0)  ! PEC
          CALL ERROR_FEMFEKO(1,4041,int1=material_type(i_mat),int2=i_mat)
       CASE (1)  ! isotropic
          WRITE(FILEOUT,'(1X,I12,2(1X,A12),2(1X,F8.3,A7,F8.3,A2))') &
               i_mat, status, 'isotropic', &
               REAL(eps_r(i_mat)), ' + i*( ', AIMAG(eps_r(i_mat)), ' )', &
               REAL(mu_r(i_mat)),  ' + i*( ', AIMAG(mu_r(i_mat)),  ' )'
       CASE (2)
          WRITE(FILEOUT,'( A,I3,A,/,2(A,3(2(G10.4,1X),2X),/)/)') &
               ' Material number ', i_mat, & 
               ' is anisotropic with diagonal tensors:', &
               ' Permittivity (xx,yy,zz): ', & 
               eps_r_xx(i_mat), eps_r_yy(i_mat),  eps_r_zz(i_mat), &
               ' Permeability (xx,yy,zz): ', &
               mu_r_xx(i_mat), mu_r_yy(i_mat),  mu_r_zz(i_mat) 
          IF ( ( ABS(eps_r_xx(i_mat)-1.0_SP).LT.EPS ) .OR. &
               ( ABS(eps_r_yy(i_mat)-1.0_SP).LT.EPS ) .OR. &
               ( ABS(eps_r_zz(i_mat)-1.0_SP).LT.EPS ) ) THEN
             WRITE(FILEOUT,*) 'WARNING. Material number ',i_mat,' has some free-space ',&
                  'dielectric material characteristics.', ' If not intended, ',&
                  'check that all components of permittivity tensor are set.'
          END IF
          IF ( ( ABS(mu_r_xx(i_mat)-1.0_SP).LT.EPS ) .OR. &
               ( ABS(mu_r_yy(i_mat)-1.0_SP).LT.EPS ) .OR. &
               ( ABS(mu_r_zz(i_mat)-1.0_SP).LT.EPS ) ) THEN
             WRITE(FILEOUT,*) 'WARNING. Material number ',i_mat,' has some free-space ',&
                  'magnetic  material characteristics.', ' If not intended, ',&
                  'check that all components of permeability tensor are set.'
          END IF
       CASE DEFAULT
          CALL ERROR_FEMFEKO(1,4041,int1=material_type(i_mat),int2=i_mat)
       END SELECT
    END DO
  END SUBROUTINE output_materials_info

  SUBROUTINE output_mesh_info
    USE geometry
    USE problem_info
    USE quad_tables
    USE matrix
    USE near_field_data
    USE boundary_conditions
    IMPLICIT NONE

    INTEGER(i4b) :: i_BC, i_vert, i_node, i_elem

    ! ANALYSIS INFORMATION:
    WRITE (FILEOUT,'(//,20X,A,/)') 'ANALYSIS INFORMATION'
    IF (CBAA_ANALYSIS) &
         WRITE (FILEOUT,'(1X,A)') 'Analysis type:          Cavity backed aperture analysis'
    IF (GW_ANALYSIS) &
         WRITE (FILEOUT,'(1X,A)') 'Analysis type:          Guided wave analysis'
    IF (REAL_EIGEN_ANALYSIS) &
         WRITE (FILEOUT,'(1X,A)') 'Analysis type:          Real eigen analysis of a cavity'
    SELECT CASE(HIERARCHAL_ORDER)
       ! Changes DBD 28 Feb 02
    CASE(1) ! CT/LN 
       IF(MIXED_ORDER_FLAG) THEN
          WRITE (FILEOUT,'(1X,A)') 'Elements used:          Constant Tangential/Linear Normal'
       ELSE 
          WRITE (FILEOUT,'(1X,A)') 'Elements used:          Linear Tangential/Linear Normal'
       END IF
    CASE(2) ! LT/QN 
       IF(MIXED_ORDER_FLAG) THEN
          WRITE (FILEOUT,'(1X,A)') 'Elements used:          Linear Tangential/Quadratic Normal'
       ELSE
          WRITE (FILEOUT,'(1X,A)') 'Elements used:          Quadratic Tangential/Quadratic Normal'
       END IF
       SELECT CASE(ELEMENT_TYPE)
       CASE(1)
          WRITE (FILEOUT,'(25X,A)') 'Savage and Peterson type'
       CASE(2)
          WRITE (FILEOUT,'(25X,A)') 'Andersen and Volakis type'
       CASE(3)
          WRITE (FILEOUT,'(25X,A)') 'Webb99 type'
       END SELECT
    END SELECT
    ! End DBD changes 28 Feb 02.
! Changes DBD 4 Aug 05

    IF (OUTPUT_ELEMENT_SHAPE) THEN
   IF (CURVILINEAR) THEN
     WRITE (FILEOUT,'(1X,A)')   'Tet shape:              2nd order curvilinear (where selected).'
   ELSE IF (ALL_CURVILINEAR) THEN
     WRITE (FILEOUT,'(1X,A)')   'Tet shape:              2nd order curvilinear (all).'
   ELSE
     WRITE (FILEOUT,'(1X,A)')   'Tet shape:              Rectilinear.'
  END IF
END IF
! End DBD changes 4 Aug 05
    
    IF (CUBATURE) THEN
       WRITE (FILEOUT,'(1X,A)') 'Elemental matrices:     Cubature used.'
    ELSE
       WRITE (FILEOUT,'(1X,A)') 'Elemental matrices:     Closed form expressions used.'
    END IF


    IF (GW_ANALYSIS) THEN
       WRITE (FILEOUT,'(1X,A,I2)') 'Number of quadrature points used for boundary integrals:',& 
            gauss_points
    END IF

    IF (SPARSE) THEN
       WRITE (FILEOUT,'(1X,A)') 'Sparse matrix storage:  Yes'
    ELSE
       WRITE (FILEOUT,'(1X,A)') 'Sparse matrix storage:  No'
    END IF
    SELECT CASE(SOLVER_TYPE)
    CASE(0)
       WRITE (FILEOUT,'(1X,A)') 'Matrix solution scheme: LU factorization (direct)'
    CASE(1)
       WRITE (FILEOUT,'(1X,A)') 'Matrix solution scheme: Bi-conjugate gradient (iterative)'
       WRITE (FILEOUT,'(1X,A,G10.4)') 'Residual norm target:   ',residual_norm
    CASE(2)
       WRITE (FILEOUT,'(1X,A)') 'Matrix solution scheme: Conjugate gradient (iterative)'
       WRITE (FILEOUT,'(1X,A,G10.4)') 'Residual norm target:   ',residual_norm
    END SELECT
    IF (FIELD_AVERAGING) THEN
       WRITE (FILEOUT,'(1X,A)') 'FE field averaging:     Yes'
    ELSE
       WRITE (FILEOUT,'(1X,A)') 'FE field averaging:     No'
    END IF

    ! BOUNDARY CONDITIONS:
    WRITE(FILEOUT,'(//,20X,A,/)') 'BOUNDARY CONDITIONS'
    WRITE(FILEOUT,'(2(1X,A14),1X,A30)') 'BC number','BC type','co-ordinates (x,y,z)'
    DO i_BC = 1,num_BCs
       SELECT CASE (BCs(i_BC)%type)
       CASE(0)
          WRITE(FILEOUT,'(1X,I14,1X,A14)',ADVANCE='NO') i_BC,'Free'
       CASE(1)
          WRITE(FILEOUT,'(1X,I14,1X,A14)',ADVANCE='NO') i_BC,'Port'
       CASE(2)
          WRITE(FILEOUT,'(1X,I14,1X,A14)',ADVANCE='NO') i_BC,'PEC'
       CASE(5)
          IF (ABCs(BCs(i_BC)%ABC_num)%type.EQ.1) THEN
             WRITE(FILEOUT,'(1X,I14,1X,A14)',ADVANCE='NO') i_BC,'ABC (1st)'
          ELSE
	     STOP 'MESHIN_FEK: Unimplemented ABC type'
          END IF
       CASE DEFAULT 
          WRITE(FILEOUT,'(I5,10X,A)') i_BC,'Unimplemented'
       END SELECT
       DO i_vert = 1,4
          WRITE(FILEOUT,'(T31,2X,A,I2,A,3(1X,F12.4))') 'corner',i_vert,':',BCs(i_BC)%corner(i_vert,1), &
               BCs(i_BC)%corner(i_vert,2),BCs(i_BC)%corner(i_vert,3)
       END DO
       IF(BCs(i_BC)%type.EQ.5) THEN 
          WRITE(FILEOUT,'(T31,2X,A,G12.5)') 'Y_c',ABCs(BCs(i_BC)%ABC_num)%Yc
       END IF
    END DO

    IF (PML_PRESENT) THEN
       WRITE(FILEOUT,'(//,20X,A,/)') 'PERFECTLY MATCHED LAYER USED'
       IF (PML%full) THEN
          WRITE(FILEOUT,'(20X,A,/)') 'FULL IMPLEMENTATION'
       ELSE
          WRITE(FILEOUT,'(20X,A,/)') 'APPROXIMATE IMPLEMENTATION'
       END IF
       WRITE(FILEOUT,'(2(1X,A20,G10.4,/),1X,A20,I2)') 'Thickness: ',  &
            PML%thickness,'Maximum sigma: ',PML%sigma_max,'Polynomial order: '&
            ,PML%m
       WRITE(FILEOUT,'(3(A11,1X))') 'Absorb x','Absorb y','Absorb z'
       WRITE(FILEOUT,'(3(5X,L1,6X))') PML%absorb_x,PML%absorb_y,PML%absorb_z
    END IF



    ! LIST OF NODES
    WRITE(FILEOUT,'(//,20X,A)') 'LIST OF NODES'
    WRITE(FILEOUT,'(/,32X,A)') 'LOCATION'
    WRITE(FILEOUT,'(4(A14))') 'Node no.', 'X', 'Y', 'Z'
    IF (OUTPUT_ELEMENT_DATA) THEN
       DO i_node = 1,NUM_nodes
          WRITE(FILEOUT,'(1X,I10,3X,3(2X,E12.5))') i_node, vertices(i_node)%coord(1:3)
       END DO
    ELSE 
       WRITE(FILEOUT,'(1X,A)') 'Suppressed by option OUTPUT_ELEMENT_DATA'
    END IF



    !    WRITE (FILEOUT,'(10X,A,I8)') 'Degrees of freedom used in FEA: ',dof
    !    WRITE (FILEOUT,'(10X,A,G10.4)') 'Average edge length in mesh: ',avg_edge_length
    IF(DEBUG_MESHIN) THEN
       WRITE(FILEOUT,*) 'Raw, unsorted, elements and associated nodes'
       DO i_elem=1,num_elements
          WRITE(FILEOUT,*) i_elem,elements(i_elem)%nodes(:)
       END DO
    END IF
  END SUBROUTINE output_mesh_info
  
  
END MODULE output_data
