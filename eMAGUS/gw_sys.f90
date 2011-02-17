! Last changed 03 May 2002 DBD. 
! NB! Routine TRANS_MATRIX corrected: line 497-499 and 602-609

MODULE gw_sys

CONTAINS

SUBROUTINE GW_SYSMAT
   USE boundary_conditions
   USE frequency_data
   USE geometry
   USE gw_data
   USE matrix
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE problem_info
   USE unit_numbers
   IMPLICIT NONE   
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine constructs and solves the system matrix for guided wave
! problems. It is based on the theory in Section 8.5, J-M Jin, "The Finite
! Element Method in Electromagnetics", Wiley 1993. 
!
! Mixed-order CT/LN and LT/QN elements have been implemented. Work presently in 
! progress on polynomial complete LT/LN and QN/QN. 
! Element definitions are now entirely contained with MODULE basisfun. 
! Presently, several different types of elements are avaiable.
!
! Changed (updated) DBD 30 March 2001
! NB: Assumptions: 
! o The waveguide may be arbitrarily orientated, and may have an 
!   arbitrary number of ports. The phase reference for each port is set
!   at the port surface, as would be the case for a measurement on a 
!   vector network analyzer. 
! o The waveguide should support only dominant TE10 mode propagation 
!   at the ports. Although the code permits
!   the dimensions to be different, the normalizing impedances are not
!   presently thus adjusted. At present, the code does not test for this. 
! o The code tests whether the w/g is operating above cut-off. However, 
!   it does NOT test whether higher-order modes, or degenerate modes, could
!   exists. 
!
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! DB Davidson
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! 24 Mar 2000: First release. DBD.
! 28 Jun 2000: First sparse matrix release. DBD.
! 13 Jul 2000: Extended to include LT/QN elements. DBD.
! 16 Jan 2001: Extended to support FEKO-like analysis - MMB. 
! 30 Mar 2001: Extended to handle multi-port, bent waveguide analysis. 
!              DBD. 
! 28 Feb 2002: Extended to support polynomail complete elements. 
! 
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data is the material composition, geometry and boundary conditions.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Output data is the edge-based degrees of freedom for subsequent 
! post-processing.
!
!*******************************************************************************
! EXTERNAL ROUTINES REQUIRED
!*******************************************************************************
! LAPACK routines CSYTRF and CSYTRS (for full matrix solution only). 
!
!*******************************************************************************

  ! Check that waveguide is operating above cut-off  
  IF (k0**2.LT.(PI/ports(1)%a)**2) THEN 
    CALL ERROR_FEMFEKO(0,4502)
  END IF
 
  ! Check if additional TE/TM modes could propagate (not a complete check, 
  ! just the most probable)
  IF (k0**2.GE.((PI/ports(1)%a)**2+(PI/ports(1)%b)**2)) THEN  ! TE/TM_11
    CALL ERROR_FEMFEKO(0,4509)
  ELSE IF (k0**2.GE.(PI/ports(1)%b)**2) THEN    ! TE/TM_01
    CALL ERROR_FEMFEKO(0,4509)
  ELSE IF (k0**2.GE.(2*PI/ports(1)%a)**2) THEN  ! TE/TM_20
    CALL ERROR_FEMFEKO(0,4509)
  ELSE IF (k0**2.GE.(2*PI/ports(1)%b)**2) THEN  ! TE/TM_02
    CALL ERROR_FEMFEKO(0,4509)
  END IF

  CALL GW_MATRIX_FILL ! Fill FEM system matrix.

CONTAINS

!*******************************************************************************

SUBROUTINE GW_MATRIX_FILL
  USE feminterface, ONLY: MAKE_FE_AMATRIX 
  USE frequency_data
  IMPLICIT NONE
!*******************************************************************************
! This routine fills the FEM system matrices, by calling a number of 
! routines that first compute, and then assemble, the elemental matrices.
!*******************************************************************************
! Note - matrix cannot be naively symmetrized, since entries are in fairly
! randomly entered in U and L parts. Entries in the corresponding
! symmetrical submatrices are entered, term by term,  as computed for  
! terms such as Se1f1 etc. (i.e. Sf1e1 in this case). 
! This is not required for sub-matrices such as Se1e1 etc (note that the 
! whole sub-matrix is presently built - for simplicity - 
! for these entries at the cost of a little computational efficiency). 
! 
! Last changed: DBD, 30 March 2001; normal vector extended. 
!*******************************************************************************
  INTEGER(I4B) ielem,jface,row,col,line,linelength
  INTEGER(I4B) m,n,itemp
  LOGICAL(LGT) matrix_symmetry

  ! (Re-)initialize arrays and matrix. (Re-initialization needed for 
  ! frequency loop).
  IF (SPARSE) THEN
    Asparse_c = ZERO_C ! Array initializations.
  ELSE 
    A_mat_c = ZERO_C 
  END IF
  x_vec_c = ZERO_C 

  ! Add the FE (S&T matrices) contribution to the system matrix:
  CALL MAKE_FE_AMATRIX(k0)

  ! Add contribution of ports to the system matrix:
  CALL GW_MAKE_PORT_AMATRIX(k0)

  DEBUG_MATRICES: IF (DEBUG_SYSTEM_MATRIX.AND..NOT.SPARSE) THEN
    ! Test matrix symmetry - asymmetry is an indication of error
    ! Only implemented for full matrices, and only for testing.
    matrix_symmetry = .TRUE.
    DO row = 1,dof
      DO col =1,dof
    IF (ABS(A_mat_c(row,col)-A_mat_c(col,row)).GT.EPS) THEN
          WRITE (FILEOUT,'(A,I4,I4,A,G10.4)') 'Matrix entry (',row,col,&
            ') asymmetry exceeds',EPS
          matrix_symmetry = .FALSE.
    END IF
      END DO
    END DO
    IF (matrix_symmetry) THEN 
      WRITE (FILEOUT,'(/A,G10.4)') 'Matrix was symmetrical within ', EPS
    END IF
    WRITE(FILEOUT,'(//A,A)') '****** [A] matrix ', & 
                           '**************'
    DO row=1,dof
      WRITE (FILEOUT,'(A,I4)') 'Row = ', row
      linelength = 10
      DO line = 1,CEILING(dof/REAL(linelength))
        DO col=(line-1)*linelength+1,line*linelength
          IF (col.GT.dof) THEN
            EXIT
          END IF
          WRITE(FILEOUT,'(A,G10.4,1X,A,G10.4,A,1X)',ADVANCE='NO') & 
              '(',REAL(A_mat_c(row,col)),'+i*',AIMAG(A_mat_c(row,col)),')'
        END DO
        WRITE (FILEOUT,'()') ! New line
      END DO
      WRITE (FILEOUT,'(/)') ! Skip a line
    END DO 
  
    WRITE(FILEOUT,'(//A,A)') '****** [b] vector ', & 
                           '**************'
    DO row=1,dof      
      WRITE(FILEOUT,'(E12.5,1X,E12.5,3X)') b_vec_c(row)       
    END DO
  END IF DEBUG_MATRICES

END SUBROUTINE GW_MATRIX_FILL
!*******************************************************************************


SUBROUTINE GW_FULL_MATRIX_LU_SOLVE(dim)
  IMPLICIT NONE
!*******************************************************************************
! This routine solves the linear system. At present, full matrix storage only
! is supported. 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: dim                 ! Dimension info
  INTEGER(I4B) lda,ldb                            ! Dimension info for 
                                                  ! CSYTRF and CSYTRS
  INTEGER(I4B) :: row                             ! Counter
  INTEGER(I4B) lwork                              ! LAPACK workspace size
  INTEGER(I4B) info                               ! Flag from 
                                                  ! CSYTRF and CSYTRS

  !********************************************************************    
  ! The LAPACK LU solver and back-substitution for full, complex-valued 
  ! matrices is used
  !********************************************************************
  lda = dim ! Set leading dimensions for LAPACK routine for 'A' 
  CALL CSYTRF('L',(dim),A_mat_c,(lda),ipiv,work_c,lwork,info)
  ! This is the LAPACK routine for symmetric complex  single
  ! precision matrices.
  ! It factors the matrix [A] 
  ! See notes in subroutine description re. obtaining and compiling.

  SELECT CASE(info) ! Check error status on exit.
    CASE(:-1) ! < 0
      CALL ERROR_FEMFEKO(1,4505,int1=ABS(info))
    CASE(0)
      WRITE(FILEOUT,'(A)') 'Routine CSYTRF exited successfully'
      WRITE(FILEOUT,'(A,I8/)') 'Optimal value of lwork :',INT(work(1))
    CASE(1:)  ! > 0
      CALL ERROR_FEMFEKO(1,4506,int1=info)
  END SELECT  

  ! Above settings for factorize are correct for solve.
  ldb=dim
  CALL CSYTRS('L',(dim),(1),A_mat_c,(lda),ipiv,b_vec_c,(ldb),info)
  x_vec_c= b_vec_c
  ! This is the LAPACK routine for symmetric complex  single
  ! precision matrices.
  ! It  solves the linear sysmtem [A][x]=[b] using the already-computed
  ! factorization of [A].
  ! Note that [b] is overwritten by [x] on exit.
  ! See notes in subroutine description re. obtaining and compiling.

  SELECT CASE(info) ! Check error status on exit.
   CASE(:-1) ! < 0
     CALL ERROR_FEMFEKO(1,4507,int1=ABS(info))
   CASE(0)
     WRITE(FILEOUT,'(A)') 'Routine CSYTRS exited successfully'
   CASE DEFAULT
     CONTINUE ! No associated error. 
  END SELECT  

  
  DEBUG_SOLUTION: IF (DEBUG_SYSTEM_MATRIX) THEN
    WRITE(FILEOUT,'(A,I8)') 'Element order: ',HIERARCHAL_ORDER
    WRITE(FILEOUT,'(//A,A)') '****** [x] vector following LAPACK solution ', & 
                           '**************'
    DO row=1,dof
      WRITE(FILEOUT,'(10(E12.5,1X))') x_vec_c(row)
    END DO 
  END IF DEBUG_SOLUTION
  ! NB -  Do NOT deallocate workspace and pivots here - needed for next frequency
  ! point.

END SUBROUTINE GW_FULL_MATRIX_LU_SOLVE
!*******************************************************************************

END SUBROUTINE GW_SYSMAT

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

SUBROUTINE GW_PREPROCESSING
  USE geometry
  USE problem_info
  USE math_tools, ONLY: CROSS_PRODUCT, VECTOR_LENGTH
  USE boundary_conditions
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Calculates the relevant data in <ports> to be used through out the code.
! Note that it is possible that the ports are not defined in the sequence
! 1,2,...etc in the input file; hence the reference via the A9data datastructure. 
! Usually, port_count  and A9data(i_port)%port_num will be the same.
!
! 2002-03-12: Correction to temporary normal calculation. MMB.
!*******************************************************************************
  INTEGER(I4B) :: i_port, port_count,BC_count,BC_num ! counters
  REAL(SP), DIMENSION(3) :: temp_normal           ! temporary storage of normal. 
  REAL(SP), DIMENSION(3) :: temp_tangent          ! temporary storage of tangent.
  REAL(SP), DIMENSION(3) :: temp_vec1,temp_vec2,temp_vec3   ! temporary storage.
  REAL(SP) a_temp,b_temp
  
  ! Assign BC_number:
  DO BC_count = 1,NUM_BCs
    IF (BCs(BC_count)%type.EQ.1) THEN ! this is a port geometry specification
      ports(BCs(BC_count)%port_num)%BC_number = BC_count
    END IF
  END DO

! Start DBD change 20 Mar 2001.
  DO i_port = 1,NUM_A9_cards ! = num_ports, checked in MESHIN. 

   port_count = A9data(i_port)%port_num

! Next lines ARE NOW INCORRECT!
!    ! Assign normal:
!    SELECT CASE(port_count) ! Assume two ports normal to ^z
!    CASE (1) 
!      ports(port_count)%normal = (/0.0_SP,0.0_SP,-1.0_SP/)
!    CASE (2) 
!      ports(port_count)%normal = (/0.0_SP,0.0_SP,+1.0_SP/)
!    CASE DEFAULT 
!      STOP 'IE: In GW_MAKE_PORT_BVECTOR, Invalid port number, currrently only 2 ports allowed.'
!    END SELECT

    ports(port_count)%normal(1:3)= A9data(i_port)%normal(1:3)
    ports(port_count)%tangent(1:3)= A9data(i_port)%tangent(1:3)

    ! Create unit normal and tangential vectors:
    temp_normal(1:3) = ports(port_count)%normal(1:3)/ &
      VECTOR_LENGTH(ports(port_count)%normal(1:3))
    temp_tangent(1:3) = ports(port_count)%tangent(1:3)/ &
      VECTOR_LENGTH(ports(port_count)%tangent(1:3))

    IF (DEBUG_GW_PRE) THEN
      print *,''
      print *,'Normal at port',port_count,'=',ports(port_count)%normal
      print *,'Unit normal at port',port_count,'=',temp_normal  
      print *,'Tangent at port',port_count,'=',ports(port_count)%tangent
      print *,'Unit tangent at port',port_count,'=',temp_tangent  
      print *,''
    END IF

    ports(port_count)%normal(1:3) = temp_normal(1:3) ! Unit normal vector
    ports(port_count)%tangent(1:3) = temp_tangent(1:3) ! Unit tangential vector

    ! Check that this user-specified normal really IS normal to the 
    ! port surface. Get two vectors lying the plane of the port, and 
    ! cross-product them with the unit vector (the answer should be zero).
    ! These vectors are arbitrarily chosen as those from corners 1 to 2 
    ! and 1 to 4. The temporary normal vector is normalized to make the test size-independent.
    ! The reason that the vector must be user-specified and not computed automatically 
    ! thus is the difficulty of making the vector always point OUTWARDS.
    BC_num = ports(port_count)%BC_number ! Find appropriate BC number (often
                                         ! the same as the port and A9 card
                                         ! counter, but not guaranteed).
    temp_vec1(1:3) = BCs(BC_num)%corner(2,1:3)-BCs(BC_num)%corner(1,1:3)
    temp_vec2(1:3) = BCs(BC_num)%corner(4,1:3)-BCs(BC_num)%corner(1,1:3)
!print *,'temp_vec1',temp_vec1
!print *,'temp_vec2',temp_vec2
    temp_normal(1:3) = CROSS_PRODUCT(temp_vec1,temp_vec2) ! Vector perp. to port
!print *,'Raw temp_normal',temp_normal
    temp_normal(1:3) = temp_normal/VECTOR_LENGTH(temp_normal) ! Temp. unit vector
    temp_vec3(1:3) = CROSS_PRODUCT(temp_normal(1:3),ports(port_count)%normal(1:3))
!print *,'Normalized temp_normal',temp_normal
!print *,'temp_vec3',temp_vec3
!print *,'length of Xproduct:',VECTOR_LENGTH(temp_vec3)
!print *,'eps',EPS

    IF (ABS(VECTOR_LENGTH(temp_vec3)).GT.EPS) THEN 
      CALL ERROR_FEMFEKO(1,4080,port_count)
    END IF

    ! Ditto the user-specified tangential vector. This is done by
    ! constructing a temporary normal, using the specified tangential vector 
    ! and a vector lying in the plane of the port  (this vector is arbitarily 
    ! chosen - again - as the vector from corners 1 to 2). 
    ! This temporary normal vector is normalized to make the test size-independent.
    ! It is then cross-producted with the normal (already tested above) 
    ! which should yield zero. 
    ! This user-specified tangential must be user-specified  to resolve the 
    ! potential phase ambiguity which arises from the definition of the 
    ! direction of the mode. 
    temp_vec1(1:3) = BCs(BC_num)%corner(2,1:3)-BCs(BC_num)%corner(1,1:3)
    temp_normal(1:3) = CROSS_PRODUCT(ports(port_count)%tangent(1:3),temp_vec1)
    
! Added MMB 2002-03-12 begin
	! Check that <temp_vec1> and <tangent> are not  in line:
	IF (VECTOR_LENGTH(temp_normal).LT.EPS) THEN 
      temp_vec1(1:3) = BCs(BC_num)%corner(4,1:3)-BCs(BC_num)%corner(1,1:3)
      temp_normal(1:3) = CROSS_PRODUCT(ports(port_count)%tangent(1:3),temp_vec1)
    END IF
! Added MMB 2002-03-12 end

    temp_normal(1:3) = temp_normal/VECTOR_LENGTH(temp_normal) ! Temp. unit vector
    temp_vec3(1:3) =    CROSS_PRODUCT(temp_normal(1:3),ports(port_count)%normal(1:3))
    IF (ABS(VECTOR_LENGTH(temp_vec3)).GT.EPS) THEN 
      CALL ERROR_FEMFEKO(1,4081,port_count)
    END IF
! End DBD change 20 Mar 2001.

    ! Assign a,b:
    CALL WG_DIMS(port_count,ports(port_count)%a,ports(port_count)%b)

    ! Find geometrical transformation matrix for port
    CALL TRANS_MATRIX(port_count)

    IF (DEBUG_GW_PRE) THEN
      CALL TRANS_TEST(port_count) ! For debugging only
    END IF 

! Start DBD change 30 Mar 2001 - following calls removed, no longer needed
! for abritrarily-orientated w/g. 
    ! Assign x0,y0:
!    CALL WG_ORIGIN(port_count,ports(port_count)%x0,ports(port_count)%y0)

    ! Assign z_coord:
!    CALL WG_ORIENTATION(port_count,ports(port_count)%z_coord)
! End DBD change 28 Mar 2001. 

  END DO

! Start DBD change 30 Mar 2001 - following lines removed, no longer needed
! for abritrarily-orientated w/g. 
!  ! Check assumption of relative position and size of ports 1 and 2:
!  IF(ports(1)%z_coord.GE.ports(2)%z_coord) THEN
!    CALL ERROR_FEMFEKO(1,4500)
!  END IF
! End DBD change 30 Mar 2001

! Start DBD change 3 Apr 2001
! Test if waveguide port dimensions are the same:
  a_temp = ports(1)%a
  b_temp = ports(1)%b
  DO port_count = 2,num_ports
    IF(ABS(ports(num_ports)%a-ports(num_ports)%b).LT.EPS) THEN
      CALL ERROR_FEMFEKO(0,4500) ! Square waveguide, degenerate modes possible.
    END IF
    IF(ABS(ports(num_ports)%a-a_temp).GE.EPS.OR.ABS(ports(num_ports)%b-b_temp).GE.EPS) THEN
      CALL ERROR_FEMFEKO(0,4501) ! Waveguide dimensions at ports not identical.
    END IF
  END DO
! End DBD change 3 Apr 2001

CONTAINS

! New subroutine
  SUBROUTINE TRANS_MATRIX(port_num)
!*******************************************************************************
! This routine computes the transformation matrix for port number port_num.
! "Local" port computations assume the port lies in the xi-eta plane, 
! centered on the zeta axis. 
! See, for eg. Rogers and Adams, " Mathematical Elements for Computer Graphics",
! McGraw-Hill, 1990. Note that for the special case when the waveguide
! lies on the x-axis, special treatment is required to avoid
! a divide-by-zero error.
!
! The transformation matrix (with no change of perspective and no overall
! scaling)
! has the form:
!
! [xi eta zeta 1] = [x y z 1]][a   b   c   | 0] 
!                             [d   e   f   | 0]
!                             [g   i   j   | 0]
!                             [---------------}
!                             [t_x t_y t_z | 1]
!
! where (t_x, t_y,t_z_) is the translation vector, and the upper left 3x3
! matrix produces linear transformations, ie scaling, shearing, reflection,
! rotation. 
!
! The approach used here is the following: (operations are performed on the 
! normal vector). 
! Firstly, construct the direction cosines of the normal vector. 
! Secondly, rotate the normal (unit) vector about the x-axis and then about
! the y-axis, as described in Rogers and Adams, p.121-123. 
! Finally, a further rotation by 180 degrees about the z axis may be needed
! to preserve the correct sense of the tangent vector. 
!
! Note the following local numbering convention in the zeta=0 plane. 
! (The numbering in the global coordinate system is arranged in the 
! same sequence. This is done in routine WG_DIMS). 
!
!            eta 
!       4 ___|____ 3      __
!        |   |    |        | 
!        |   |----|-->xi   b    normal (zeta) is out of the page
!        |        |        |
!       1|--------|2      --
!
!        <---a --->
!
! eta has the same interpretation as the tangent vector, i.e. defines
! the positive sense of the eigenmode.
!
! At the same time, the inverse transformation matrix is also computed - i.e.
! to take a point from the local port (xi,eta,zeta) coordinate system back 
! to the Cartesian coordinate system. 
! 
! Author: DB Davidson. First version  02 April 2001.
! Corrected DBD 4 June 2001.
!
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: port_num
    REAL(SP), DIMENSION(4,4) :: M,T,Rx,Ry,Rz,temp_mat
    REAL(SP), DIMENSION(4,4) :: M_inv,T_inv,Rx_inv,Ry_inv,Rz_inv
    REAL(SP) c_x, c_y, c_z                ! Direction cosines. 
    REAL(SP), DIMENSION(3) :: n_x,n_y,n_z ! Unit normal vectors
    REAL(SP) x_0, y_0, z_0                ! Translation vector
    REAL(SP) d                            ! Projection of unit vector on yz plane
    REAL(SP) t_x,t_y                      ! Projection of tangent vector onto
                                          ! x and y axes. 
    INTEGER(I4B) BC_num
    REAL(SP), DIMENSION(4) :: x_cnr
    REAL(SP), DIMENSION(4) :: y_cnr
    REAL(SP), DIMENSION(4) :: z_cnr
    REAL(SP), DIMENSION(4) :: coords
 ! Changed DBD 02 May 2002
    REAL(SP), DIMENSION(4) :: tx_coords1,tx_coords2
 ! End DBD changes
    REAL(SP), DIMENSION(3) :: tx_tangent,x_hat,y_hat

    ! Firstly, check that the normal and tangent vectors are normalized:
    IF(ABS(VECTOR_LENGTH(ports(port_num)%normal ))-1.0_SP.GT.EPS.OR.&
       ABS(VECTOR_LENGTH(ports(port_num)%tangent))-1.0_SP.GT.EPS) THEN
      STOP 'IE: Internal error in TRANS_MATRIX, unit vector(s) un-normalized'
    END IF

    BC_num = ports(port_num)%BC_number
    x_cnr(1:4) = BCs(BC_num)%corner(1:4,1)
    y_cnr(1:4) = BCs(BC_num)%corner(1:4,2)  
    z_cnr(1:4) = BCs(BC_num)%corner(1:4,3)
    ! Find centre of port (average of each coordinate)
    ! This defines the translation vector
    x_0 = (x_cnr(1) + x_cnr(2) + x_cnr(3) + x_cnr(4))/4.0_SP
    y_0 = (y_cnr(1) + y_cnr(2) + y_cnr(3) + y_cnr(4))/4.0_SP
    z_0 = (z_cnr(1) + z_cnr(2) + z_cnr(3) + z_cnr(4))/4.0_SP
    IF(DEBUG_GW_PRE) print *,'Translation vector',-x_0,-y_0,-z_0

    ! Find direction cosines.
    n_x = (/1.0_SP,0.0_SP,0.0_SP/)
    n_y = (/0.0_SP,1.0_SP,0.0_SP/)
    n_z = (/0.0_SP,0.0_SP,1.0_SP/)
    c_x = DOT_PRODUCT(ports(port_num)%normal,n_x)
    c_y = DOT_PRODUCT(ports(port_num)%normal,n_y)
    c_z = DOT_PRODUCT(ports(port_num)%normal,n_z)


    ! Find projection of unit vector onto yz-plane
    d = SQRT(c_y**2 + c_z**2)
    IF(DEBUG_GW_PRE) print *,'d: ',d

    ! Construct T, Rx and Ry matrices (and also their inverses)
    T      = 0.0_SP ! Array assignment
    Rx     = 0.0_SP ! Ditto
    Ry     = 0.0_SP ! Ditto
    Rz     = 0.0_SP ! Ditto
    T_inv  = 0.0_SP ! Array assignment
    Rx_inv = 0.0_SP ! Ditto
    Ry_inv = 0.0_SP ! Ditto
    Rz_inv = 0.0_SP ! Ditto
    
    ! Translation matrix: 
    T(1,1) =  1.0_SP
    T(2,2) =  1.0_SP
    T(3,3) =  1.0_SP
    T(4,1) = -x_0
    T(4,2) = -y_0
    T(4,3) = -z_0
    T(4,4) =  1.0_SP

    ! Inverse translation matrix: 
    T_inv(1,1) =  1.0_SP
    T_inv(2,2) =  1.0_SP
    T_inv(3,3) =  1.0_SP
    T_inv(4,1) = +x_0
    T_inv(4,2) = +y_0
    T_inv(4,3) = +z_0
    T_inv(4,4) =  1.0_SP


    ! Compute rotation (and inverse rotation) matrix about x-axis.
    IF (ABS(d).GT.EPS) THEN 
      Rx(1,1) =  1.0_SP
      Rx(2,2) =  c_z/d
      Rx(2,3) =  c_y/d
      Rx(3,2) = -c_y/d
      Rx(3,3) =  c_z/d
      Rx(4,4) =  1.0_SP
      Rx_inv(1,1) =  1.0_SP
      Rx_inv(2,2) =  c_z/d
      Rx_inv(2,3) = -c_y/d
      Rx_inv(3,2) = +c_y/d
      Rx_inv(3,3) =  c_z/d
      Rx_inv(4,4) =  1.0_SP
    ELSE ! Already aligned with x axis, special case (d=0). 
         ! No rotation necessary.
      Rx(1,1) = 1.0_SP
      Rx(2,2) = 1.0_SP
      Rx(3,3) = 1.0_SP
      Rx(4,4) = 1.0_SP
      Rx_inv(1,1) = 1.0_SP
      Rx_inv(2,2) = 1.0_SP
      Rx_inv(3,3) = 1.0_SP
      Rx_inv(4,4) = 1.0_SP
    END IF 

    ! Compute rotation (and inverse rotation) matrix about y-axis.
    Ry(1,1) =  d
    Ry(1,3) =  c_x
    Ry(2,2) =  1.0_SP
    Ry(3,1) = -c_x
    Ry(3,3) =  d
    Ry(4,4) =  1.0_SP
    Ry_inv(1,1) =  d
    Ry_inv(1,3) = -c_x
    Ry_inv(2,2) =  1.0_SP
    Ry_inv(3,1) = +c_x
    Ry_inv(3,3) =  d
    Ry_inv(4,4) =  1.0_SP

    ! Compute rotation (and inverse rotation) about the z-axis in the xy plane. 
    temp_mat = MATMUL(Rx,Ry)
    temp_mat = MATMUL(T,temp_mat)
	! DBD corrections 2 May 2002
    ! Now, transform tangent vector.
    coords = (/ ports(port_num)%tangent(1:3), 1.0_SP /)
    tx_coords1 = MATMUL(coords,temp_mat) ! Transform head of vector
    coords = (/ 0.0_SP, 0.0_SP, 0.0_SP, 1.0_SP /)
    tx_coords2 = MATMUL(coords,temp_mat) ! Transform tail of vector
    tx_tangent = tx_coords1(1:3) - tx_coords2(1:3) ! New transformed tangent vector
	! End DBD corrections 2 May 2002
    IF(DEBUG_GW_PRE) THEN 
      print *,'Temporary transformation matrix'
      print *,temp_mat(1,1:4)
      print *,temp_mat(2,1:4)
      print *,temp_mat(3,1:4)
      print *,temp_mat(4,1:4)
      print *,'Transformed tangent vector',tx_tangent
    END IF
    x_hat = (/1.0_SP,0.0_SP,0.0_SP/) ! Unit vector in x (xi) direction
    y_hat = (/0.0_SP,1.0_SP,0.0_SP/) ! Unit vector in y (eta) direction
    ! Now find projections of unit vector onto x and y axes.
    t_x = DOT_PRODUCT(tx_tangent,x_hat)
    t_y = DOT_PRODUCT(tx_tangent,y_hat)
    ! Note that the angle of rotation is defined by the angle
    ! relative to the Y axis. Thus cos psi = t_y/1, sin psi = t_x/1
    IF(DEBUG_GW_PRE) THEN 
      print *,'Further rotation about z at port',port_num, &
              'by ',ATAN2(t_x,t_y)*180.0_SP/PI,'degrees'
    END IF
    Rz(1,1) = t_y
    Rz(1,2) = t_x
    Rz(2,1) = -t_x
    Rz(2,2) = t_y
    Rz(3,3) = 1.0_SP
    Rz(4,4) = 1.0_SP
    Rz_inv(1,1) = t_y
    Rz_inv(1,2) = -t_x
    Rz_inv(2,1) = t_x
    Rz_inv(2,2) = t_y
    Rz_inv(3,3) = 1.0_SP
    Rz_inv(4,4) = 1.0_SP

    ! Now construct full transformation & inverse transformation matrices 
    temp_mat = MATMUL(Ry,Rz)
    temp_mat = MATMUL(Rx,temp_mat)
    M = MATMUL(T,temp_mat)

    temp_mat = MATMUL(Rx_inv,T_inv)
    temp_mat= MATMUL(Ry_inv,temp_mat)
    M_inv =  MATMUL(Rz_inv,temp_mat)

    ! Save these transformation and inverse transformation matrices for the port.
    ports(port_num)%T = M
    ports(port_num)%T_inv = M_inv
  END SUBROUTINE TRANS_MATRIX

! New subroutine
  SUBROUTINE TRANS_TEST(port_num)
!*******************************************************************************
! This is a test routine to check the geometrical transformation and 
! inverse transformation routines. 
! Only called when debugging option DEBUG_GW_PRE on. 
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: port_num
    REAL(SP), DIMENSION(4,4) :: M
    INTEGER(I4B) BC_num, cnr
    REAL(SP), DIMENSION(4) :: coord, tx_coord ! Coordinates and transformed 
                                              ! coordinates.

    BC_num = ports(port_num)%BC_number
    print *,'***********************************************************'

    DO cnr= 1,4
      coord(1:3) = BCs(BC_num)%corner(cnr,1:3)
      coord(4) = 1.0_SP
      tx_coord(1:4) = MATMUL(coord,ports(port_num)%T)
      print *,'Original coordinates of port number ',port_num, 'corner ',cnr,&
              '(after sorting)'
      print *,coord(1:3)
      print *,'Transformed coordinates' 
      print *,tx_coord(1:3)
      ! Now test inverse transform
      coord(1:4) = MATMUL(tx_coord,ports(port_num)%T_inv)
      print *,'Inverse transformed coordinates' 
      print *,coord(1:3)
    END DO
    print *,'***********************************************************'
  END SUBROUTINE TRANS_TEST


! Routine changed DBD 23 March 2001
  SUBROUTINE WG_DIMS(port_num,a,b)
    USE boundary_conditions
    USE math_tools, ONLY: SORT_VECTOR
    IMPLICIT NONE   
  !*******************************************************************************
  ! This internal routine returns the waveguide width a and height b at port 
  ! number port_num; the conventional assumption of a >= b is made.
  ! The algorithm choses one corner (#1, an abitrary choice) and then computes
  ! the lengths to the other three points. One of these is the diagonal, one 
  ! the width and one the height. Which is which is established by a sort
  ! operation.
  ! Once this is established, the corners are then sorted in accordance with the 
  ! following convention: 
  !
  !            tangent
  !           / \ 
  !       4 ___|____ 3    __
  !        |   |    |      | 
  !        |   0    |      b    normal is out of the page
  !        |        |      |
  !       1|--------|2    --
  !
  !        <---a --->
  !
  ! This is accomplished by swopping corners until the above is obtained. 
  ! This is required for subsequent routines, to ensure for instance 
  ! that the positive sense of modes is correctly retained.  
  !
  !*******************************************************************************
     INTEGER(I4B), INTENT(IN) :: port_num
     REAL(SP),INTENT(OUT) :: a,b

     INTEGER(I4B) BCnum
     REAL(SP), DIMENSION(3) :: ell
     REAL(SP), DIMENSION(4) :: x,y,z
     REAL(SP), DIMENSION(3) :: vec12,vec23,vec24,vec34,vec41

     BCnum = ports(port_num)%BC_number

     IF (DEBUG_GW_PRE) THEN 
       print *,'Port',port_num,' corners before sorting:'
       print *,BCs(BCnum)%corner(1,1:3)
       print *,BCs(BCnum)%corner(2,1:3)
       print *,BCs(BCnum)%corner(3,1:3)
       print *,BCs(BCnum)%corner(4,1:3)
     END IF

     x(1:4) = BCs(BCnum)%corner(1:4,1)
     y(1:4) = BCs(BCnum)%corner(1:4,2)
     z(1:4) = BCs(BCnum)%corner(1:4,3)
     ell(1)  = SQRT((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
     ell(2)  = SQRT((x(3)-x(1))**2+(y(3)-y(1))**2+(z(3)-z(1))**2)
     ell(3)  = SQRT((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)
     CALL SORT_VECTOR(ell,3)! Sorts ell into ascending order. 
     b = ell(1)
     a = ell(2)
     ports(port_num)%b = ell(1)
     ports(port_num)%a = ell(2)

     ! Check on local numbering and re-number as required.  
     ! Starting with corner 1, the possible sequences are:
     ! 12XX
     ! 13XX
     ! 14XX
     ! where XX will be determined subsequently.
     ! Firstly, get numbering 1-2 aligned with width: 
     ! Construct a vector from corner 1 to corner 2
     vec12(1:3) = BCs(BCnum)%corner(2,1:3) - BCs(BCnum)%corner(1,1:3)
     ! Then check its length to see if it is equal to width a: 
     ALIGN12: IF (ABS(VECTOR_LENGTH(vec12)-a).LT.EPS) THEN 
       CONTINUE! No action needed, edge suitably numbered. 
     ELSE 
       ! Swop corners 2 and 3 and try again:
       CALL SWOP_CORNER(BCnum,2,3)
       vec12(1:3) = BCs(BCnum)%corner(2,1:3) - BCs(BCnum)%corner(1,1:3)
       IF (ABS(VECTOR_LENGTH(vec12)-a).LT.EPS) THEN
         CONTINUE! No further action needed here. 
       ELSE
         ! Swop back
         CALL SWOP_CORNER(BCnum,3,2)
         ! Swop corners 2 and 4 and try again:
         CALL SWOP_CORNER(BCnum,2,4)
         vec12(1:3) = BCs(BCnum)%corner(2,1:3) - BCs(BCnum)%corner(1,1:3)
         IF (ABS(VECTOR_LENGTH(vec12)-a).LT.EPS) THEN
           CONTINUE ! Finished. 
         ELSE ! Shouldn't happen.
           STOP 'IE: WG_DIMS: Corner sorting failed: 1st stage' 
         END IF
       END IF
     END IF ALIGN12

     ! Now, get numbering 2-3 aligned with height:  
     vec23(1:3) = BCs(BCnum)%corner(3,1:3) - BCs(BCnum)%corner(2,1:3)
     ! Check its length to see if it is equal to height b. 
     ALIGN23: IF (ABS(VECTOR_LENGTH(vec23)-b).LT.EPS) THEN 
       CONTINUE! No action needed, edge suitably numbered. 
     ELSE 
       CALL SWOP_CORNER(BCnum,3,4)
       vec23(1:3) = BCs(BCnum)%corner(3,1:3) - BCs(BCnum)%corner(2,1:3)
       IF (ABS(VECTOR_LENGTH(vec23)-b).LT.EPS) THEN
         CONTINUE! No further action needed. 
       ELSE ! Shouldn't happen.
         STOP 'IE: WG_DIMS: Corner sorting failed: 2nd stage' 
       END IF
     END IF ALIGN23

     vec23(1:3) = BCs(BCnum)%corner(3,1:3) - BCs(BCnum)%corner(2,1:3)
     ! Align local "sense" with the tangent vector:
     IF (DOT_PRODUCT(vec23,ports(port_num)%tangent).GT.0.0_SP) THEN
       CONTINUE ! Alignment correct
     ELSE ! Interchange corners 1 and 4 and 2 and 3. 
       CALL SWOP_CORNER(BCnum,1,4)
       CALL SWOP_CORNER(BCnum,2,3)
     END IF     

     ! Finally, check that the numbering is counter-clockwise as above:
     ! (right-handed w.r.t. normal)
     vec12(1:3) = BCs(BCnum)%corner(2,1:3) - BCs(BCnum)%corner(1,1:3)
     IF (DOT_PRODUCT(CROSS_PRODUCT(vec12,ports(port_num)%tangent), & 
       ports(port_num)%normal).GT.0.0_SP) THEN
       CONTINUE ! Alignment correct
     ELSE ! Interchange corners 1 and 2 and 3 and 4. 
       CALL SWOP_CORNER(BCnum,1,2)
       CALL SWOP_CORNER(BCnum,3,4)
     END IF     

     IF (DEBUG_GW_PRE) THEN 
       print *,'Port',port_num,' corners after sorting:'
       print *,BCs(BCnum)%corner(1,1:3)
       print *,BCs(BCnum)%corner(2,1:3)
       print *,BCs(BCnum)%corner(3,1:3)
       print *,BCs(BCnum)%corner(4,1:3)
     END IF

     ! Do a final check to ensure that the sorting has worked:
     vec12(1:3) = BCs(BCnum)%corner(2,1:3) - BCs(BCnum)%corner(1,1:3)
     vec23(1:3) = BCs(BCnum)%corner(3,1:3) - BCs(BCnum)%corner(2,1:3)
     vec34(1:3) = BCs(BCnum)%corner(4,1:3) - BCs(BCnum)%corner(3,1:3)
     vec41(1:3) = BCs(BCnum)%corner(1,1:3) - BCs(BCnum)%corner(4,1:3)
     IF (ABS(VECTOR_LENGTH(vec12)-a).GE.EPS.OR. &
         ABS(VECTOR_LENGTH(vec23)-b).GE.EPS.OR. & 
         ABS(VECTOR_LENGTH(vec34)-a).GE.EPS.OR. & 
         ABS(VECTOR_LENGTH(vec41)-b).GE.EPS) THEN
       CALL ERROR_FEMFEKO(1,4510,port_num)
     END IF


  END SUBROUTINE WG_DIMS
! End routine changed DBD 23 March 2001

! New routine added DBD 23 March 2001
  SUBROUTINE SWOP_CORNER(BCnum,corner1,corner2)
    IMPLICIT NONE 
  !*****************************************************************************
  ! This internal subroutine swops two corners around for boundary condition
  ! BC_num
  ! DBD, 23 March 2001.
  !*****************************************************************************
    INTEGER(I4B), INTENT(IN):: BCnum, corner1,corner2
    REAL(SP), DIMENSION(3) :: corner_temp
    corner_temp(1:3) = BCs(BCnum)%corner(corner1,1:3)
    BCs(BCnum)%corner(corner1,1:3) = BCs(BCnum)%corner(corner2,1:3)
    BCs(BCnum)%corner(corner2,1:3) = corner_temp(1:3)
  END SUBROUTINE SWOP_CORNER
!*******************************************************************************

END SUBROUTINE GW_PREPROCESSING
!*******************************************************************************


SUBROUTINE GW_S_PARAMS(excited_port)
  USE feminterface, ONLY: GAUSS_QUAD_CFIELD
  USE boundary_conditions
  USE frequency_data
  USE geometry
  USE gw_data
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! This routine computes the S parameters for this port. 
! It is multi-port, arbitrarily orientated port extension of the equations given
! for reflection and transmission
! coefficients  at port 1 by [Jin, eqns.(8.91-8.92), p. 266]. 
! First version:
! Jun 2000. DBD.
! Revised:
! 13 Jul 2000: LT/QN support added. DBD.
! 29 March 2001: arbitrarily orientated, multi-ports added. DBD. 
! 28 Feb 2002: LT/LN and QT/QN support added. DBD. 
! 
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: excited_port ! excited port.  
  COMPLEX(SPC) R,T
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: int
  REAL(SP), DIMENSION(6) :: ell
  REAL(SP), DIMENSION(3) :: ell_tri
  COMPLEX(SPC), DIMENSION(6):: e1,e2,e3
  COMPLEX(SPC), DIMENSION(4):: f1,f2,f3
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE):: temp_dofs_tri
  COMPLEX(SPC), DIMENSION(3):: e1_tri,e2_tri,e3_tri
  COMPLEX(SPC) :: f1_tri,f2_tri, f3_tri
  INTEGER(I4B) :: i_elem,j_face,k_edge,m_port
  COMPLEX(SPC), DIMENSION(3) :: E_mid
  COMPLEX(SPC) E0,E_y_z1,E_y_z2
  
  ALLOCATE (int(num_ports))


  PORT_LOOP: DO m_port = 1,num_ports
    int(m_port) = ZERO_C
    ELEMENT_LOOP: DO i_elem = 1,num_elements
      FACE_LOOP: DO j_face = 1,4
        ! (Re-)initialize local variables:
		e1  = ZERO_C ! Arrays 
        e2  = ZERO_C
		e3  = ZERO_C
        f1  = ZERO_C ! Scalars
        f2  = ZERO_C
		f3  = ZERO_C
        PORT_TEST: IF (faces(elements(i_elem)%faces(j_face))%portnumber.EQ.m_port) THEN
          ! Find the d.o.f.'s for the tetrahedron. Note that a face may 
          ! be a port but an edge nonetheless may be prescribed. 
          ! Note that a dof value is allocated for ALL hierarchal orders;
          ! if this exceeds the hierarchal order which has been used
		  ! it remains as initialized to zero
		  ! (as for a prescribed dof) and hence one need not 
		  ! explicitly differentiate between the various orders
		  ! when integrating.
	      
          ! Store the three edge lengths of this face:
          ell =  T_LENGTHS(i_elem)
          ell_tri = ell(LOCAL_FACEEDGES(j_face))
          
          ! CT/LN and higher elements:
 
		  ! Store the e1 dof values in <e1>:
          DO k_edge = 1,6
            IF (renumbered_e1(elements(i_elem)%edges(k_edge)).EQ.0) THEN
              e1(k_edge) = ZERO_C !prescribed edge
            ELSE IF (renumbered_e1(elements(i_elem)%edges(k_edge)).GT.dof) THEN
              STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. e1.'
            ELSE 
              e1(k_edge) = x_vec_c(renumbered_e1(elements(i_elem)%edges(k_edge)))
            END IF 
          END DO

          ! LT/LN, LT/QN and higher elements:
          IF (MAX_ORDER(i_elem,j_face).GE.2.OR.&
		    .NOT.MIXED_ORDER(i_elem,j_face).AND.MAX_ORDER(i_elem,j_face).GE.1) THEN
            DO k_edge = 1,6
              IF (renumbered_e2(elements(i_elem)%edges(k_edge)).EQ.0) THEN
                e2(k_edge) = ZERO_C !prescribed edge
              ELSE IF (renumbered_e2(elements(i_elem)%edges(k_edge)).GT.dof) THEN
                STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. e2.'
              ELSE 
                e2(k_edge) = x_vec_c(renumbered_e2(elements(i_elem)%edges(k_edge)))
              END IF 
            END DO
		  END IF
		   
		  ! LT/QN, QT/QN and higher elements
          IF (MAX_ORDER(i_elem,j_face).GE.2) THEN
            IF (renumbered_f1(elements(i_elem)%faces(j_face)).EQ.0) THEN
              f1(j_face) = ZERO_C !prescribed edge
            ELSE IF (renumbered_f1(elements(i_elem)%faces(j_face)).GT.dof) THEN
              STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. f1.'
            ELSE 
              f1(j_face) = x_vec_c(renumbered_f1(elements(i_elem)%faces(j_face)))
            END IF 
            IF (renumbered_f2(elements(i_elem)%faces(j_face)).EQ.0) THEN
              f2(j_face) = ZERO_C !prescribed edge
            ELSE IF (renumbered_f2(elements(i_elem)%faces(j_face)).GT.dof) THEN
              STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. f2.'
            ELSE 
              f2(j_face) = x_vec_c(renumbered_f2(elements(i_elem)%faces(j_face)))
            END IF 
	      END IF

          ! QT/QN, QT/CuN and higher elements:
          IF (MAX_ORDER(i_elem,j_face).GE.3.OR.&
		    .NOT.MIXED_ORDER(i_elem,j_face).AND.MAX_ORDER(i_elem,j_face).GE.2) THEN
            DO k_edge = 1,6
              IF (renumbered_e3(elements(i_elem)%edges(k_edge)).EQ.0) THEN
                e3(k_edge) = ZERO_C !prescribed edge
              ELSE IF (renumbered_e3(elements(i_elem)%edges(k_edge)).GT.dof) THEN
                STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. e3.'
              ELSE 
                e3(k_edge) = x_vec_c(renumbered_e3(elements(i_elem)%edges(k_edge)))
              END IF 
            END DO
			IF (renumbered_f3(elements(i_elem)%faces(j_face)).EQ.0) THEN
              f3(j_face) = ZERO_C !prescribed edge
            ELSE IF (renumbered_f3(elements(i_elem)%faces(j_face)).GT.dof) THEN
              STOP 'IE: In GW_S_PARAMS, ref. to out of range d.o.f. f3.'
            ELSE 
              f3(j_face) = x_vec_c(renumbered_f3(elements(i_elem)%faces(j_face)))
            END IF 
		  END IF
          temp_dofs_tri = TET2TRI_DOF(j_face,e1,e2,e3,f1,f2,f3)
          e1_tri(1:3) = temp_dofs_tri(1:3)
          e2_tri(1:3) = temp_dofs_tri(4:6)
		  f1_tri = temp_dofs_tri(7)
          f2_tri = temp_dofs_tri(8)
          e3_tri(1:3) = temp_dofs_tri(9:11)
          f3_tri = temp_dofs_tri(12)
          int(m_port) = int(m_port) + & 
            GAUSS_QUAD_CFIELD(m_port,i_elem,j_face,e1_tri,&
            e2=e2_tri,e3=e3_tri,f1=f1_tri,f2=f2_tri,f3=f3_tri,&
			ell=ell_tri)
          IF (MAX_ORDER(i_elem,j_face).GE.4.OR.&
		    .NOT.MIXED_ORDER(i_elem,j_face).AND.MAX_ORDER(i_elem,j_face).GE.3) THEN
            STOP 'IE: Unimplemented hierarchal order in GW_S_PARAMS.'
          END IF
        END IF PORT_TEST
      END DO FACE_LOOP
    END DO ELEMENT_LOOP

    ! Now fill S matrix. 
    IF (m_port.EQ.excited_port) THEN ! Compute R (over port m_port):
      S_params(m_port,excited_port) = &
        2.0_SP*int(m_port)/(ports(m_port)%a*ports(m_port)%b*ports(excited_port)%excitation) & 
        - 1.0_SP
    ELSE ! Compute T (over port m_port): 
      S_params(m_port,excited_port) = &
        2.0_SP*int(m_port)/(ports(m_port)%a*ports(m_port)%b*ports(excited_port)%excitation) 
    END IF

  END DO PORT_LOOP

  ! Outputting moved to subroutine OUTPUT_S_PARAMS in post_pro.f90
  DEALLOCATE (int)

CONTAINS

  FUNCTION TET2TRI_DOF(jface,e1,e2,e3,f1,f2,f3)
    USE nrtype
    USE problem_info
    IMPLICIT NONE
  !*******************************************************************************
  ! This internal function returns: 
  !  - For CT/LN elements, the appropriate three degrees of freedom for 
  ! the triangle corresponding to face jface, given the six d.o.f.'s
  ! for the parent tetrahedron-based edge numbering scheme.  
  !  - For LT/LN elements, the appropriate six as above. 
  !  - For LT/QN elements, the appropriate eight degrees of freedom for 
  ! the triangle corresponding to face jface, given the twenty d.o.f.'s
  ! for the parent tetrahedron-based edge numbering scheme. 
  !  - For QT/QN elements, the appropriate 12 dofs from the 30. 
  ! etc. 
  !
  ! The numbering conversion scheme, using the S&P edge and face numbering scheme,
  ! is as for TET2TRI_EDGENUMS.
  ! Extended DBD 2 Mar 2002 to include LT/LN, and QT/QN
  !*******************************************************************************
    INTEGER (I4B), INTENT(IN):: jface
    COMPLEX (SPC), DIMENSION(6), INTENT(IN):: e1
    COMPLEX (SPC), DIMENSION(6), INTENT(IN), OPTIONAL :: e2
	COMPLEX (SPC), DIMENSION(6), INTENT(IN), OPTIONAL :: e3
    COMPLEX (SPC), DIMENSION(4), INTENT(IN), OPTIONAL :: f1,f2,f3
    COMPLEX (SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: TET2TRI_DOF
    INTEGER(I4B), DIMENSION(3) :: tempfaceedges

    ! Initialize:
    TET2TRI_DOF(1:ELEM_TRI_MATRIX_SIZE) = ZERO_C

    ! Find local edges:
    tempfaceedges = LOCAL_FACEEDGES(jface)

    ! e1 related dofs:
    TET2TRI_DOF(1:3) = e1(tempfaceedges)

    ! e2 related dofs:
    IF (PRESENT(e2)) THEN
      TET2TRI_DOF(4:6) = e2(tempfaceedges)
    END IF

    ! f1 related dofs:
    IF (PRESENT(f1)) THEN
      TET2TRI_DOF(7) = f1(jface)
    END IF
    
    ! f2 related dofs:
    IF (PRESENT(f2)) THEN
      TET2TRI_DOF(8) = f2(jface)
    END IF

    ! e3 related dofs:
    IF (PRESENT(e3)) THEN
      TET2TRI_DOF(9:11) = e3(tempfaceedges)
    END IF
    
	! f3 related dofs:
    IF (PRESENT(f3)) THEN
      TET2TRI_DOF(12) = f3(jface)
    END IF

  END FUNCTION TET2TRI_DOF

END SUBROUTINE GW_S_PARAMS
!*******************************************************************************


SUBROUTINE GW_MAKE_PORT_BVECTOR
  USE boundary_conditions
  USE feminterface, ONLY: U_MAKE_HIERARCHAL, &
                          LOCAL_TO_GLOBAL_INDEX_TRI
  USE frequency_data
  USE geometry
  USE gw_data
  USE matrix
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  IMPLICIT NONE
!*******************************************************************************
! Adds the contribution of all ports to the excitation vector.
! Note that for S parameter analysis, only ONE of the ports is excited at 
! a time. Also, the phase references planes per port are now the port surfaces 
! (as for network analyzer measurement).
! It is possible that more than one port may be excited simultaneously 
! when NOT performing S parameter analysis.
!
! This routine calls LOCAL_TO_GLOBAL_INDEX_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating them. 
!
! Changed DBD 29 March 2001 for S parameter analysis.
! Reworked by MMB to decrease size. 2001-10-12.
! Minor changes for QT/QN elements DBD. 28 Feb 2002. 
! Minor changed DBD 31 March 2003.
!*******************************************************************************
  INTEGER(I4B) ielem,iface,row,port_num,idof
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE) ::  Us 
  LOGICAL(LGT) port_face_found
  COMPLEX(SPC) :: E0                  ! complex size of the excitation
  REAL(SP) :: kzmn

  ELEMENT_LOOP: DO ielem = 1,num_elements
    ! Note that an element should NOT have more than one face on a port,
    ! nor have faces on more than one port.
    port_face_found = .FALSE.
    FACE_LOOP: DO iface = 1,4
      PORT_TEST: IF (faces(elements(ielem)%faces(iface))%port) THEN ! This face is part of a port
        
        ! Error: this is the second face of this element in a port
        IF (port_face_found) CALL ERROR_FEMFEKO(1,4504,int1=ielem)

        ! Record the relevant post face information:
        port_face_found = .TRUE.
        port_num = faces(elements(ielem)%faces(iface))%portnumber
        kzmn = K_Z_MN(k0,port_num,mode_m,mode_n)
        E0 = ports(port_num)%excitation

        ! Calculate the elemental contribution:
! Changed DBD 31 March 2003
        CALL U_MAKE_HIERARCHAL(ielem,iface,MAX_ORDER(ielem,iface),port_num,& 
                               ports(port_num)%normal,Us)
! End changed DBD 31 March 2003
        ! Add this face's contribution to <b_vec_c>:
        ROW_LOOP: DO idof = 1,ELEM_TRI_MATRIX_SIZE 

          row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,idof)
          IF (row.GT.0) THEN ! this is a dof
             b_vec_c(row) = b_vec_c(row) + Us(idof)*(-2.0)*j*kzmn*E0
          END IF

        END DO ROW_LOOP
      END IF PORT_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE GW_MAKE_PORT_BVECTOR
!*******************************************************************************

SUBROUTINE GW_MATRIX_ALLOCATE
  USE feminterface, ONLY: MATRIX_SPARSE_ALLOCATE
  USE matrix
  USE output_error, ONLY: ERROR_FEMFEKO
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine allocates storage for the guided wave analysis. Full and sparse 
! storage is now implemented.
! Original DBD. Made external - MMB 7 Feb 2001
!*******************************************************************************
  IF (BANDRENUM_STORE) THEN
    CALL ERROR_FEMFEKO(1,4503)
  END IF

  IF (SPARSE) THEN
    CALL MATRIX_SPARSE_ALLOCATE
  ELSE
    ALLOCATE(A_mat_c(dof,dof))
    ! Allocate some storage for LAPACK pivots and workspace
    ALLOCATE(ipiv(dof)) ! dof = dim in LAPACK notation
    !lwork = 5*dof ! guess block size for initialization
!    READ(INFILE,NML=LAPACKDATA)   
!    ALLOCATE(work_c(lwork) )
  END IF

  ALLOCATE(x_vec_c(dof))
  ALLOCATE(b_vec_c(dof))

END SUBROUTINE GW_MATRIX_ALLOCATE
!*******************************************************************************


SUBROUTINE PORT_BOUNDARY_SEARCH
  USE boundary_conditions
  USE geometry
  USE nrtype
  USE output_error, ONLY: ERROR_FEMFEKO
  IMPLICIT NONE   
!*******************************************************************************
! Cycles through all exterior faces (unconnected faces) and flag them as part 
! of a port if appropriate.
! Original version: MMB 13 Nov 2000. 
! Revised: DBD, 7 Apr 2001. Multi-port capability added. 
!                           Test added to ensure that at least one face 
!                           is found on each port. (This is not a complete
!                           test that the port definition is valid!)
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem,port_counter
  INTEGER(I4B), DIMENSION(3) :: tempedges 
  LOGICAL(LGT) :: in_port_plane
  LOGICAL(LGT), ALLOCATABLE :: port_found(:)    ! A flag to check that a port
                                                ! does indeed lie on the mesh. 

  ALLOCATE(port_found(num_ports))
  port_found = .FALSE. ! Array allocation

  ELEMENT_LOOP: DO ielem = 1,num_elements
    FACE_LOOP: DO iface = 1,4 ! Face-wise search to find edges in the quad.
      FACE_CONNECT_TEST: IF (elements(ielem)%connect2elem(iface).EQ.0) THEN     
        tempedges = LOCAL_FACEEDGES(iface)

        PORT_LOOP: DO port_counter = 1,num_ports

          ! Check if in port (source):
          CALL FACE_IN_QUADRILATERAL( & 
               BCs(ports(port_counter)%BC_number)%corner(1,1:3), &
               BCs(ports(port_counter)%BC_number)%corner(2,1:3), &
               BCs(ports(port_counter)%BC_number)%corner(3,1:3), &
               BCs(ports(port_counter)%BC_number)%corner(4,1:3), &
               elements(ielem)%faces(iface),in_port_plane)
          IF (in_port_plane) THEN

            ! Perform check to see if this port shares an edge with 
            ! another port:
            IF ( &
              ((edges(elements(ielem)%edges(tempedges(1)))%portnumber.NE.port_counter).AND. &
               (edges(elements(ielem)%edges(tempedges(1)))%portnumber.NE.0)).OR.            &
              ((edges(elements(ielem)%edges(tempedges(2)))%portnumber.NE.port_counter).AND. &
               (edges(elements(ielem)%edges(tempedges(2)))%portnumber.NE.0)).OR.            &
              ((edges(elements(ielem)%edges(tempedges(3)))%portnumber.NE.port_counter).AND. &
               (edges(elements(ielem)%edges(tempedges(3)))%portnumber.NE.0))                &
              ) THEN
                CALL ERROR_FEMFEKO(1,4512,port_counter)
            END IF

            ! Flag the three edges as part of port port_counter
            edges(elements(ielem)%edges(tempedges))%port = .TRUE.
            edges(elements(ielem)%edges(tempedges))%portnumber = port_counter
            faces(elements(ielem)%faces(iface))%port = .TRUE.
            faces(elements(ielem)%faces(iface))%portnumber = port_counter
            port_found(port_counter) = .TRUE. ! Will usually be overwritten 
          END IF

        END DO PORT_LOOP
      END IF FACE_CONNECT_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

  ! Check that at least one element (face) was found per port.
  DO port_counter = 1,num_ports
    IF(.NOT.port_found(port_counter)) THEN
      CALL ERROR_FEMFEKO(1,4511,port_counter)
    END IF
  END DO
  DEALLOCATE(port_found)
    
END SUBROUTINE PORT_BOUNDARY_SEARCH
!*******************************************************************************


FUNCTION CONVERT_FROM_LOCAL_PORT_COORD(vector,port_num)
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Convert vector expressed in local coordinates on port port_num
! to global (x,y,z) Cartesian system.
! The algorithm uses the inverse translation matrix for the port, applied 
! to the origin and the point described by "vector". 
! The vector in the Cartesian coordinate system is then the difference
! of these vector components. 
! 
! DBD, 10 April 2001.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) ::   port_num
  REAL(SP), DIMENSION(3), INTENT(IN) :: vector ! Global coordinates  
  REAL(SP), DIMENSION(3) :: CONVERT_FROM_LOCAL_PORT_COORD       
  REAL(SP), DIMENSION(4) :: origin
  REAL(SP), DIMENSION(4) :: tx_origin
  REAL(SP), DIMENSION(4) :: vector_head
  REAL(SP), DIMENSION(4) :: tx_vector_head

  origin = (/0.0_SP,0.0_SP,0.0_SP,1.0_SP/)
  tx_origin = MATMUL(origin,ports(port_num)%T_inv)
  vector_head = (/vector,1.0_SP/)
  tx_vector_head = MATMUL(vector_head,ports(port_num)%T_inv)
  CONVERT_FROM_LOCAL_PORT_COORD(1:3) = tx_vector_head(1:3) - tx_origin(1:3)
END FUNCTION CONVERT_FROM_LOCAL_PORT_COORD
!*******************************************************************************


FUNCTION K_Z_MN(k_0,port_num,m,n)
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE 
!*******************************************************************************
! A function to compute the propagation constant for a rectangular waveguide
! mode, in a waveguide with width (larger dimension) a, and height b
! at free space wavenumber k_o(=2*pi/lambda_0).
! Note that it assumes that the waveguide is operating ABOVE cut-off.
! March 2000, DBD.
! March 2001, DBD: Extended for higher-order cylindrical (rectangular) 
!             waveguide modes
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k_0
  INTEGER(I4B), INTENT(IN):: port_num,m,n
  REAL(SP) K_Z_MN
  REAL(SP) a,b
  a = ports(port_num)%a
  b = ports(port_num)%b
  IF (m.NE.1.AND.n.NE.0) THEN ! Higher order mode
    print *,'This higher-order functionality has not yet been tested in FUNCTION K_Z_MN'
  END IF
  K_Z_MN = SQRT( k_0**2 - (m*PI/a)**2 - (n*PI/b)**2)
END FUNCTION K_Z_MN
!*******************************************************************************


FUNCTION E_MN(x,y,z,port_num,m,n)
  USE boundary_conditions
  USE nrtype
  IMPLICIT NONE 
!*******************************************************************************
! A function to compute the TE_MN waveguide eigenmode analytically at point
! (x,y,z) on port port_num. The point is transformed into 
! the local coordinate system, aligned 
! (as shown in TRANS_MATRIX) as:
!
!            eta 
!       4 ___|____ 3      __
!        |   |    |        | 
!        |   |----|-->xi   b    normal (zeta) is out of the page
!        |        |        |
!       1|--------|2      --
!
!        <---a --->
!
! eta has the same interpretation as the tangent vector, i.e. defines
! the positive sense of the eigenmode. 
! Most textbooks position the origin at corner 1; a simple coordinate shift
! by +a/2 and +b/2 accomplishes this. 
!
! The function then converts the vector from the local (xi,eta,zeta)
! representation back into the global Cartesian (x,y,z) coordinate system
! for further use in the code. 
!
! March 2000, DBD.
! Note that the TE eigenmodes have Ex and Ey out of phase (see, for example
! Ramo, Whinnery van Duzer, "Fields and Waves in Communication Electronics", 
! Wiley  1994, p. 421.
! Revised 28 March, 2 April 10 April 2001, DBD.
!*******************************************************************************
  REAL(SP), INTENT(IN) :: x,y,z ! Coordinates. 
  INTEGER(I4B), INTENT(IN):: port_num,m,n
  REAL(SP), DIMENSION(3):: E_MN
  REAL(SP), DIMENSION(4) :: coords
  REAL(SP), DIMENSION(4) :: tx_coords
  REAL(SP) a,b ! Width and height as above.
  REAL(SP) xi,eta,zeta ! Transformed coordinates
  a = ports(port_num)%a
  b = ports(port_num)%b
  coords = (/x,y,z,1.0_SP/)
  tx_coords = MATMUL(coords,ports(port_num)%T)
  xi = tx_coords(1)
  eta = tx_coords(2)
  IF (m.NE.1.AND.n.NE.0) THEN ! Higher order mode
    print *,'This higher-order functionality has not yet been tested in FUNCTION E_MN'
  END IF
  E_MN(1) =   cos(m*PI*(xi+a/2)/a) * sin(n*PI*(eta+b/2)/b) ! xi component
  E_MN(2) = - sin(m*PI*(xi+a/2)/a) * cos(n*PI*(eta+b/2)/b) ! eta component
  E_MN(3) =   0.0_SP                                       ! zeta component 
  ! Now convert back to Cartesian coordinates. 
  E_MN    = CONVERT_FROM_LOCAL_PORT_COORD(E_MN,port_num)
END FUNCTION E_MN
!*******************************************************************************


SUBROUTINE GW_MAKE_PORT_AMATRIX(k0)
  USE B_matrix, ONLY: B_MAKE_HIERARCHAL
  USE boundary_conditions
  USE feminterface, ONLY: CONVERTCOR, LOCAL_TO_GLOBAL_INDEX_TRI
  USE geometry
  USE gw_data
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This subroutine adds the contribution of the ports to the system matrix.
! The matrix <Bs> is the contribution of a specific port face to the system matrix,
! thus its rows and columns refer to local face quantities (ie. 3 edges and 1 face).
! These local face quantities must be converted to local element quantities (ie. 
! 6 edges and 4 faces) in order that they can be used to obtain the correct dof
! number from the indices renumbered_.. . This conversion is trivial for the face
! quantity, because the local element face number of the port face is know. In the
! case of the edges, the index connecting the local tri edge numbers to local tet
! edge number is obtained simply be calling LOCAL_FACEEDGES and remembering that
! the there local tri edge numbers (1,2 and 3) correspond to the three local tet 
! edge numbers in ascending order (eg. 2, 3, 6).
!
! This routine calls LOCAL_TO_GLOBAL_INDEX_TRI, which returns either 
! the relevant (global) row or column, or zero if the dof is either
! prescribed OR does not exist. This permits the routine to handle
! various element types (both w.r.t. order and mixed/complete order)
! without explicitly differentiating them. 
!
!
! MMB reworked original code by DBD. 2001-10-01
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k0  ! wavenumber at which calculations must be performed

  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) ::  Bs 
  LOGICAL(LGT) :: port_face_found
  INTEGER(I4B) :: ielem,face_num,jface,port_num
  INTEGER(I4B) :: mm,nn,row,column,indpos
  INTEGER(I4B), DIMENSION(3) :: localfaceedges
  REAL(SP), DIMENSION(3) :: normal    ! Outward directed unit normal on port.
  COMPLEX(SPC) :: gamma    ! See [Jin,p264]
  REAL(SP) :: kzmn         ! ditto

  ! Compute the elemental [B] matrix if this element has a face on a port. 
  ! Note that an element should NOT have more than one face on a port,
  ! nor have faces on more than one port.

  ! Cycle through all element faces and add the contributions of those 
  ! that lie on a port:
  ELEMENT_LOOP: DO ielem = 1,num_elements

    ! Initialize for this element:
    port_face_found = .FALSE.
    face_num = 0 ! default (to flag no port for assembly routines)

    FACE_LOOP: DO jface = 1,4
      PORT_TEST: IF (faces(elements(ielem)%faces(jface))%port) THEN

        ! Error check and record the port face's data:
        IF (port_face_found) CALL ERROR_FEMFEKO(1,4504,int1=ielem)
        port_face_found = .TRUE.
        port_num = faces(elements(ielem)%faces(jface))%portnumber
        face_num = jface
        localfaceedges = LOCAL_FACEEDGES(face_num)
        kzmn = K_Z_MN(k0,port_num,mode_m,mode_n)
        gamma = j*kzmn

        ! Set up outward-directed normal:
        normal = ports(port_num)%normal

        ! Calculate the elemental port face matrix:
        CALL B_MAKE_HIERARCHAL(ielem,jface,normal,Bs)

        ! Now add <Bs> to the system matrix:
        DO mm = 1,ELEM_TRI_MATRIX_SIZE ! row counter

          ! Assign the global row:
          row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,face_num,mm)

          DO nn = 1,ELEM_TRI_MATRIX_SIZE ! column counter

            ! Assign the global column:
            column = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,face_num,nn)

            ! Add this elemental matrix element to the system matrix:
            IF ((row.GT.0).AND.(column.GT.0)) THEN ! this is two dofs:
              IF (SPARSE) THEN
                indpos = CONVERTCOR(row,column)
                Asparse_c(indpos) = Asparse_c(indpos) + gamma*Bs(mm,nn)
              ELSE
                A_mat_c(row,column)  = A_mat_c(row,column) + gamma*Bs(mm,nn)
              END IF
            END IF

          END DO
        END DO
      END IF PORT_TEST
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE GW_MAKE_PORT_AMATRIX
!*******************************************************************************

END MODULE gw_sys
