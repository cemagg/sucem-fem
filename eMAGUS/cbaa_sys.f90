! Problem at around:
! Line 478
! Edited out call to SPHBES. Temp. hack during SGI port DBD 8 Feb 05'

! line  1176: 
! Had to edit out call. DBD 8 Feb 05. Still needs fixing.



! 2002-05-09: Changes to CBAA_MAKE_COAX_MODE_INTEGRAL, CBAA_COAX_POSTPRO,
! CBAA_MAKE_COAX_BVECTOR.

! NB!! Matthys - you still need to change CBAA_COUNT_APERTURE_DOFS
! for polynomial complete elements. 

! Last changed: 10 Apr 2001 MMB. Added FMM storage allocation capability 
!               for higher order elements.

MODULE CBAA_SYS

CONTAINS

SUBROUTINE CBAA_SYSMAT
   USE coax_feed
   USE feminterface, ONLY: MAKE_FE_AMATRIX
   USE frequency_data
   USE geometry
   USE matrix
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE problem_info
   USE unit_numbers
   USE CBAA_data
   IMPLICIT NONE   
!*******************************************************************************
! This subroutine implements the cavity in an infinite 
! ground plane formulation presented in 'The Finite Element Method in 
! Electromagnetics', chapter 9, by J. Jin, 1993.
! The system matrix is constructed.
! 
! MM Botha
!*******************************************************************************
! REVISION HISTORY
!*******************************************************************************
! Created: Jan 17, 2000, DBD: Matthys - a template for you to start work with...
! 
! 15 March 2000: Original working version by MMB and DBD.
!
! 19 March 2000: Working program, partly integrated with rest of femfeko. 
!                Validated with results from NASA Technical paper 3544. MMB
!
!*******************************************************************************
! PROGRAM OUTLINE:
!*******************************************************************************
! 2000-01-21:
! 0. Rudimentary error/flag checking.
! 0.5. Identifying edges/faces, comprising the aperture, in terms of their
!      global edge/face numbers.
! 1. Logical flags is set containing edges' states: edges%free
! 2. Set up integer arrays that establish the S/T-matrix position as a function
!    of the global edge(/face in future for higher order elements) number.
! 3. Set up the system matrices with their FEM contributions.
! 4. Add to the system matrices the MoM contributions.
! 5. Sets up the b-vector.
! 6. Solves Ax=b (currently the full matrix is inverted)
! 7. Post processing.
!*******************************************************************************
   INTEGER(I4B) :: row,column          ! more counters
   INTEGER(I4B), DIMENSION(3) :: tempedges          ! temp storage
   INTEGER(I4B), DIMENSION(3) :: tempfacenodes      ! temp storage
   INTEGER(I4B), DIMENSION(2) :: tempnodes          ! temp storage

   ! Constants:
   COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)        ! The complex constant j

   ! Initialise the matrices: (many elements of <A_mat_c> will remain zero)
   IF (.NOT.SPARSE) THEN
     A_mat_c   = (0.0,0.0)
   ELSE
     Asparse_c = (0.0,0.0)
   END IF

   CALL MAKE_FE_AMATRIX(k0)                          ! Calculate the FE contributions to A-matrix 
   CALL CBAA_MAKE_BE_AMATRIX                         ! Calculate the BE contributions to A-matrix
   IF (CBAA_FMM_storage) CALL CBAA_MAKE_FMM_MATRICES ! Fill the FMM matrices
   IF (NUM_AX_cards.GT.0) THEN
     IF (COAX_WHITNEY) THEN
       CALL CBAA_MAKE_COAX_A_AND_B(0) ! Add the coax-feed's contribution to the A-matrix (Whitney approx)
     ELSE     
       CALL CBAA_MAKE_COAX_AMATRIX(k0) ! Add the coax-feed's contribution to the A-matrix (General)
     END IF
   END IF
   IF (CBAA_FMM_debug) CALL CBAA_FMM_DEBUGGER

   RETURN
!*******************************************************************************

   CONTAINS


SUBROUTINE CBAA_MAKE_BE_AMATRIX
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI,    &
                          CONVERTCOR
  USE geometry
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! This routine calculates the CBAA BI contribution to the system matrix and 
! adds it to the system matrix.
!
! 2002-05-13: Created from old CBAA_MAKE_BE_AMATRIX. MMB.
!*******************************************************************************
  LOGICAL(LGT) :: test_free1,test_free2                    ! Flag to test if both MoM elements have any dofs
  LOGICAL(LGT), DIMENSION(8,8) :: near_interact            ! store type of edge/face function interactions for FMM case
  LOGICAL(LGT) :: any_fmm_interactions                     ! true if any exist for (ielem,jelem)
  INTEGER(I4B) :: ielem,jelem,fmm1,fmm2,apgroup1,apgroup2, &
	              jj,kk,val_col_temp                       ! counters
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) :: run_tot  ! Max. dofs per facet = 8
  INTEGER(I4B) :: Pst_dim_row,Pst_dim_col                  ! row/col dimension of a facet-facet BI interaction matrix
  INTEGER(I4B) :: outer_hierarchal_order
  INTEGER(I4B) :: inner_hierarchal_order

  ! Array initialisation:
  IF (SPARSE) THEN
    IF (.NOT.CBAA_FMM_storage) THEN
      CBAA_BE_mat = (0.0,0.0)
    ELSE
      CBAA_BE_val = (0.0,0.0)
      IF (CBAA_FMM_debug) CBAA_BE_mat = (0.0,0.0)
    END IF
  END IF

  ! Now calculate and add the MoM contributions to the A-matrix for every possible 
  ! combination of CBAA elemental faces:
  MoM_contrib1: DO ielem = 1,num_apelements

    PRINT *,'Current CBAA face:', ielem, '/', num_apelements

    ! Test whether this aperture facet has any dofs:
    tempedges = LOCAL_FACEEDGES(which_local_face(ielem))
    test_free1 =                                                                   &
      faces(elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))%free.OR. &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(1)))%free.OR.            &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(2)))%free.OR.            &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(3)))%free
       
    ! If no dofs, then no contribution to [P] possible from this element:
    IF (.NOT.test_free1) CYCLE MoM_contrib1 

    ! Record the hierarchal order of the outer integral element:
    outer_hierarchal_order = MAX_ORDER(ap_elnumbers(ielem),which_local_face(ielem))

    ! Choose the row dimension of the matrix [P^{st}] for repeated later use:
    SELECT CASE (outer_hierarchal_order)
    CASE (1) !CT/LN
      IF (MIXED_ORDER(ap_elnumbers(ielem),which_local_face(ielem))) THEN
	    Pst_dim_row = 3
	  ELSE 
	    Pst_dim_row = 6
	  END IF
    CASE (2) ! LT/QN
      IF (MIXED_ORDER(ap_elnumbers(ielem),which_local_face(ielem))) THEN
	    Pst_dim_row = 8
	  ELSE 
	    Pst_dim_row = 12
	  END IF
    CASE DEFAULT
      STOP 'IE: invalid hierarchal order in CBAA_MAKE_BE_AMATRIX.'
    END SELECT

    ! The A-matrix is symmetric, thus only half of the
    ! MoM contributions have to be calculated.
    MoM_contrib2: DO jelem = ielem,num_apelements
 
      ! Check whether this element does have associated dofs - if not,
      ! then it does not contribute to the system matrix, and the program cycles
      ! to the next pair.
      tempedges = LOCAL_FACEEDGES(which_local_face(jelem))
      test_free2 =                                                                   &
        faces(elements(ap_elnumbers(jelem))%faces(which_local_face(jelem)))%free.OR. &
        edges(elements(ap_elnumbers(jelem))%edges(tempedges(1)))%free.OR.            &
        edges(elements(ap_elnumbers(jelem))%edges(tempedges(2)))%free.OR.            &
        edges(elements(ap_elnumbers(jelem))%edges(tempedges(3)))%free
      IF (.NOT.test_free2) CYCLE MoM_contrib2

      ! Record the hierarchal order of the inner integral element:
      inner_hierarchal_order = MAX_ORDER(ap_elnumbers(jelem),which_local_face(jelem))

      ! Choose the row dimension of the matrix [P^{st}] for repeated later use:
      SELECT CASE (inner_hierarchal_order)
      CASE (1) !CT/LN
        IF (MIXED_ORDER(ap_elnumbers(jelem),which_local_face(jelem))) THEN
          Pst_dim_col = 3 
	    ELSE 
	      Pst_dim_col = 6
	    END IF
      CASE (2) ! LT/QN
        IF (MIXED_ORDER(ap_elnumbers(jelem),which_local_face(jelem))) THEN
          Pst_dim_col = 8
	    ELSE 
	      Pst_dim_col = 12
	    END IF
      CASE DEFAULT
        STOP 'IE: invalid hierarchal order in CBAA_MAKE_BE_AMATRIX.'
      END SELECT

      ! Set up near interaction data and check if any near interactions do occur:
      IF (CBAA_FMM_storage) THEN
        any_fmm_interactions = .FALSE.
        DO fmm1 = 1,Pst_dim_row
          ! Find group that this row belongs to:
          SELECT CASE (fmm1)
          CASE (1:3)
            apgroup1 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
              edgenum=elements(ap_elnumbers(ielem))%edges(which_local_edges(ielem,fmm1)))
          CASE (4:6)
            apgroup1 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
              edgenum=elements(ap_elnumbers(ielem))%edges(which_local_edges(ielem,fmm1-3)))
          CASE (7:8)
            apgroup1 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
              facenum=elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))
          CASE DEFAULT
			STOP 'IE: invalid hierarchal order in FMM interaction in CBAA_MAKE_BE_AMATRIX.'
	      END SELECT
          DO fmm2 = 1,Pst_dim_col
            ! Find group that this column belongs to:
            SELECT CASE (fmm2)
            CASE (1:3)
              apgroup2 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
                edgenum=elements(ap_elnumbers(jelem))%edges(which_local_edges(jelem,fmm2)))
            CASE (4:6)
              apgroup2 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
                edgenum=elements(ap_elnumbers(jelem))%edges(which_local_edges(jelem,fmm2-3)))
            CASE (7:8)
              apgroup2 = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
                facenum=elements(ap_elnumbers(jelem))%faces(which_local_face(jelem)))
            CASE DEFAULT
			  STOP 'IE: invalid hierarchal order in FMM interaction in CBAA_MAKE_BE_AMATRIX.'
            END SELECT
            near_interact(fmm1,fmm2) = NEAR_NEIGHBOURS(apgroup1,apgroup2,gr_dim,gr_numx)
            any_fmm_interactions = any_fmm_interactions.OR.near_interact(fmm1,fmm2)
          END DO
        END DO
        IF (.NOT.CBAA_FMM_debug) THEN
          IF (.NOT.any_fmm_interactions) CYCLE MoM_contrib2
        END IF
      END IF

      ! Calculate the current, elemental face pair's contribution to the system matrix:
      CALL CBAA_MAKE_BI_MATRIX_ELEMENTAL(ap_elnumbers(ielem),which_local_face(ielem),Pst_dim_row,7, &
                                         ap_elnumbers(jelem),which_local_face(jelem),Pst_dim_col,7, &
										 run_tot(1:Pst_dim_row,1:Pst_dim_col))

      ! Now add the elements of run_tot to the system matrix:
      DO jj = 1,Pst_dim_row

        row = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(ielem),which_local_face(ielem),jj)

        DO kk = 1,Pst_dim_col

          column = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(jelem),which_local_face(jelem),kk)

          IF ((row.GT.0).AND.(column.GT.0)) THEN ! Check that both edges/faces are free

            IF (.NOT.SPARSE) THEN
              A_mat_c(row,column) = A_mat_c(row,column) + run_tot(jj,kk)
              IF (ielem.NE.jelem) THEN 
                ! The A-matrix is symmetric, thus only half of the
                ! MoM contributions have to be calculated
                A_mat_c(column,row) = A_mat_c(column,row) + run_tot(jj,kk)
              END IF
               
            ELSE IF ((.NOT.CBAA_FMM_storage).OR.CBAA_FMM_debug) THEN 
              CBAA_BE_mat(row,column) = CBAA_BE_mat(row,column) + run_tot(jj,kk)
              IF (ielem.NE.jelem) THEN 
                ! The A-matrix is symmetric, thus only half of the
                ! MoM contributions have to be calculated
                CBAA_BE_mat(column,row) = CBAA_BE_mat(column,row) + run_tot(jj,kk)
              END IF
            END IF ! CBAA_FMM_storage.AND.SPARSE
               
            IF (CBAA_FMM_storage) THEN ! CBAA_FMM_storage.AND.SPARSE

              IF (near_interact(jj,kk)) THEN
                val_col_temp = CONVERTCOR(row,column,2)
                IF (val_col_temp.LT.1) STOP 'IE: Invalid co-ordinate into CBAA_BE_val.'
                CBAA_BE_val(val_col_temp) = CBAA_BE_val(val_col_temp) + run_tot(jj,kk)
                IF (ielem.NE.jelem) THEN 
                  ! The A-matrix is symmetric, thus only half of the
                  ! MoM contributions have to be calculated
                  val_col_temp = CONVERTCOR(column,row,2)
                  IF (val_col_temp.LT.1) STOP 'IE: Invalid co-ordinate into CBAA_BE_val.'
                  CBAA_BE_val(val_col_temp) = CBAA_BE_val(val_col_temp) + run_tot(jj,kk)
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO MoM_contrib2
  END DO MoM_contrib1

END SUBROUTINE CBAA_MAKE_BE_AMATRIX
!*******************************************************************************


   SUBROUTINE CBAA_MAKE_FMM_MATRICES
     USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI
     USE geometry
     USE quad_tables
     IMPLICIT NONE
!*******************************************************************************
! Calculation of elements of [V][T][Vt] martices. 
!*******************************************************************************
     INTEGER(I4B) :: kcount,igroup,jgroup,grtemp,ordercount,idof,row,ielem ! counters
     REAL(SP) :: divbasistemp                       ! kernel factor in [V] - basis function divergeance
     REAL(SP), DIMENSION(3) :: basistemp            ! kernel factor in [V] - basis function
     REAL(SP), DIMENSION(3) :: diff_val,gr_centre,gr_centre2 ! group centre location
     COMPLEX(SPC) :: exptemp                        ! temp exponential term
     REAL(SP) :: kXmm                               ! argument of Hankel function
     REAL(SP) :: jn,yn,djn,dyn                      ! arguments of SPHBES
     REAL(SP) :: Leg_x,Leg0,Leg1,Leg2               ! temp Legendre value storage
     REAL(SP) :: Leg_sign                           ! sign of P_n(-x)
     COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: hankeltable ! look-up table for Hankel function values
     COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: exptable    ! look-up table for the term i^l(2l+1)
     REAL(SP), DIMENSION(25,3) :: gauss_xyz         ! NB - the first dimension MUST be >= MAX(tri quad points) ever used!
     REAL(SP), DIMENSION(6) :: eval_t1              ! ^z cross (E1 and E2)'s divergeance
     REAL(SP), DIMENSION(25,2) :: gauss_eval_t1     ! ^z cross (F1 and F2)'s divergeance
     REAL(SP), DIMENSION(25,8,3) :: gauss_eval_t2   ! NB - the first dimension MUST be >= MAX(tri quad points) ever used!
     INTEGER(I4B) :: Pst_dimension                  ! dimension of a facet-facet BI interaction matrix
     INTEGER(I4B) :: num_qpoints, qcount            ! number of quadrature points for triangular surface integration
     INTEGER(I4B) :: face_hierarchal_order
     REAL(SP), DIMENSION(25,ELEM_TRI_MATRIX_SIZE)   :: temp_quad_term1
     REAL(SP), DIMENSION(25,ELEM_TRI_MATRIX_SIZE,3) :: temp_quad_term2

     !Temp vars! - maybe permanent ...
     REAL(SP) :: g1x,g1y,g2x,g2y,diffx,diffy
     

     ! Choose the number of quadrature points:
     num_qpoints = 7

     ! Init the value matrices:
     CBAA_groupmat = (0.0,0.0)
     CBAA_elemvec_term1 = (0.0,0.0)
     CBAA_elemvec_term2 = (0.0,0.0)
 
     ! Calculation of [V]:
     ELEMENT_LOOP: DO ielem = 1,num_apelements

       ! Choose the dimension of the matrix [P^{st}] for repeated later use:
       face_hierarchal_order = MAX_ORDER(ap_elnumbers(ielem),which_local_face(ielem))
       SELECT CASE (face_hierarchal_order)
       CASE (1) !CT/LN
         Pst_dimension = 3 ! 3 E1 in the face
       CASE (2) ! LT/QN
         Pst_dimension = 8 ! 3 E1, 3 E2, 1 F1 and 1 F2
       END SELECT

       ! eval_t1 and eval_t2: Evaluate the basis funtions here for every gauss point in the
       ! the current element:
       SELECT CASE (face_hierarchal_order)
       CASE (1)
          CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),&
               num_qpoints,3, gauss_xyz(1:num_qpoints,1:3), &
               temp_quad_term1(1:num_qpoints,1:3),temp_quad_term2(1:num_qpoints,1:3,1:3))
          eval_t1(1:3) = (1.0/quad_tri_rules(num_qpoints)%rule(1,4))*temp_quad_term1(1,1:3)
         gauss_eval_t2(1:num_qpoints,1:3,1:3) = temp_quad_term2(1:num_qpoints,1:3,1:3)
       CASE (2)
         CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,8, &
		                                   gauss_xyz(1:num_qpoints,1:3), &
		                                   temp_quad_term1(1:num_qpoints,1:8),temp_quad_term2(1:num_qpoints,1:8,1:3))
         eval_t1(1:3)                         = (1.0/quad_tri_rules(num_qpoints)%rule(1,4))*temp_quad_term1(1,1:3)
         eval_t1(4:6)                         = (1.0/quad_tri_rules(num_qpoints)%rule(1,4))*temp_quad_term1(1,4:6)
		 gauss_eval_t1(1:num_qpoints,1:2)     = temp_quad_term1(1:num_qpoints,7:8)
         gauss_eval_t2(1:num_qpoints,1:8,1:3) = temp_quad_term2(1:num_qpoints,1:8,1:3)
       CASE DEFAULT
         STOP 'IE: invalid hierarchal order in CBAA_MAKE_FMM_MATRICES.'
       END SELECT

       FUNCTION_LOOP: DO idof = 1,Pst_dimension
         ! Find the dof number of the current basis function:
         row = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(ielem),which_local_face(ielem),idof)
         IF (row.EQ.0) CYCLE FUNCTION_LOOP

         ! Calculate the group centre:
         SELECT CASE (idof)
         CASE (1:3)
           grtemp = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
             edgenum=elements(ap_elnumbers(ielem))%edges(which_local_edges(ielem,idof)))
         CASE (4:6)
           grtemp = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
             edgenum=elements(ap_elnumbers(ielem))%edges(which_local_edges(ielem,idof-3)))
         CASE (7:8)
           grtemp = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx, &
             facenum=elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))
         CASE DEFAULT
           STOP 'IE: invalid local dimension in CBAA_MAKE_FMM_MATRICES.'
         END SELECT
         gr_centre(2) = ap_dims(2) + gr_dim*REAL(CEILING(REAL(grtemp)/REAL(gr_numx))) - gr_dim/2.0
         gr_centre(1) = ap_dims(1) + &
           gr_dim*REAL(grtemp - gr_numx*(CEILING(REAL(grtemp)/REAL(gr_numx))-1)) - gr_dim/2.0
!print *,'gr_centre =',gr_centre

         GAUSS_LOOP: DO qcount = 1,num_qpoints
           ! Calculate the difference between current gausspoint and gr_centre:
           diff_val = 0.0
           diff_val(1:2) = gauss_xyz(qcount,1:2) - gr_centre(1:2)
!print *,'difference =',diff_val

           ! Calculate the basis function factors in the kernels:
           SELECT CASE (idof)
           CASE (1:6)
             basistemp(1:3) = gauss_eval_t2(qcount,idof,1:3)
             divbasistemp   = quad_tri_rules(num_qpoints)%rule(qcount,4)*eval_t1(idof) ! must weigh manually, because this is a constant
           CASE (7:8)
             basistemp(1:3) = gauss_eval_t2(qcount,idof,1:3)
             divbasistemp   = gauss_eval_t1(qcount,idof-6)
           CASE DEFAULT
             STOP 'IE: invalid local dimension in CBAA_MAKE_FMM_MATRICES.'
           END SELECT

           KDIR_LOOP: DO kcount = 1,num_k_dirs
             exptemp = EXP(-cj*k0*SUM(k_dirs(kcount,1:2)*diff_val(1:2)))
             CBAA_elemvec_term1(kcount,row)     = CBAA_elemvec_term1(kcount,row)      &
                                                   + exptemp*divbasistemp
             CBAA_elemvec_term2(kcount,row,1:2) = CBAA_elemvec_term2(kcount,row,1:2)  &
                                                   + exptemp*basistemp(1:2)
           END DO KDIR_LOOP
         END DO GAUSS_LOOP
       END DO FUNCTION_LOOP
     END DO ELEMENT_LOOP

     ! Calculation of [T]:
     ALLOCATE(hankeltable(0:L_max))
     ALLOCATE(exptable(0:L_max))
     DO ordercount = 0,L_max
       ! Calculate the term (-i)^l*(2l+1) and store in look-up table:
       exptable(ordercount) = ((-cj)**ordercount)*(2.0*REAL(ordercount)+1.0)
     END DO

     Grouploop1: DO igroup = 1,max_groups ! Main loops over global group numbers
       IF (gr_matindex(igroup).EQ.0) CYCLE Grouploop1
       
       ! Calculate group centre(igroup):
       gr_centre(2) = ap_dims(2) + gr_dim*REAL(CEILING(REAL(igroup)/REAL(gr_numx))) - gr_dim/2.0
       gr_centre(1) = ap_dims(1) + &
         gr_dim*REAL(igroup - gr_numx*(CEILING(REAL(igroup)/REAL(gr_numx))-1)) - gr_dim/2.0

       g1x = gr_centre(1)   
       g1y = gr_centre(2)   

       Grouploop2: DO jgroup = igroup,max_groups ! calculation properties are symmetric - but NOT [T].
         IF (gr_matindex(jgroup).EQ.0) CYCLE Grouploop2
         IF (NEAR_NEIGHBOURS(igroup,jgroup,gr_dim,gr_numx)) CYCLE Grouploop2
         
         ! Calculate group centre(jgroup):
         gr_centre2(2) = ap_dims(2) + gr_dim*REAL(CEILING(REAL(jgroup)/REAL(gr_numx))) - gr_dim/2.0
         gr_centre2(1) = ap_dims(1) + &
           gr_dim*REAL(jgroup - gr_numx*(CEILING(REAL(jgroup)/REAL(gr_numx))-1)) - gr_dim/2.0

         g2x = gr_centre2(1)   
         g2y = gr_centre2(2)   

!print *, 'gr_centre2 =',gr_centre2


         ! Calculate Xmm':
!        gr_centre = gr_centre - gr_centre2 ! gr_centre(1:2) now contains (x,y) of Xmm'  
         diffx = g1x-g2x
         diffy = g1y-g2y
!         gr_centre = (/ g1x-g2x , g1y-g2y /)

!print *,'difference =',gr_centre
!print *,'ABS(difference) =',SQRT(SUM(gr_centre**2))


         ! Calculate the argument of the Hankel function:
         kXmm = k0*SQRT(diffx**2+diffy**2)

         DO ordercount = 0,L_max
           ! Calculate the spherical Hankel function values and store in look-up table:
           STOP 'Supply routine SPHBES and fix. Temp. hack during SGI port DBD 8 Feb 05'
!           CALL SPHBES(ordercount,kXmm,jn,yn,djn,dyn)
           hankeltable(ordercount) = jn - cj*yn
         END DO


         DO kcount = 1,num_k_dirs
           Leg_x = k0*(k_dirs(kcount,1)*diffx+k_dirs(kcount,2)*diffy)/kXmm
!print *,'Leg_x =',Leg_x,'   kXmm =',kXmm,'     gr-diff =',diffx,',',diffy
           Leg1 = 0.0
           Leg2 = 1.0
           DO ordercount = 0,L_max
             ! Calculate the Legendre polynomial value:
             IF (ordercount.GT.0) THEN
               Leg0 = Leg1
               Leg1 = Leg2
               Leg2 = ((2.0*REAL(ordercount)-1.0)*Leg_x*Leg1 -                      &
                                    (REAL(ordercount)-1.0)*Leg0)/REAL(ordercount)
             END IF
             ! <Leg2> contains the Legendre polynomial of order <ordercount>, 
             ! evaluated at <Leg_x>. The recurrence relation used can be found in
             ! 'Numerical Recipes in FORTRAN' or Arfken-'Math. Methods for Physicists'.
             
             ! Add to the [T] matrix:
             ! add to position (gr_matindex(igroup),gr_matindex(jgroup)):
             CBAA_groupmat(kcount,gr_matindex(igroup),gr_matindex(jgroup)) =              &
                     CBAA_groupmat(kcount,gr_matindex(igroup),gr_matindex(jgroup))        & 
                     + (exptable(ordercount) * hankeltable(ordercount) * Leg2)
             ! add to position (gr_matindex(jgroup),gr_matindex(igroup)):
             Leg_sign = 1.0
             IF ((REAL(ordercount)/2.0-REAL(FLOOR(REAL(ordercount)/2.0))).GT.EPS) Leg_sign = -1.0
             CBAA_groupmat(kcount,gr_matindex(jgroup),gr_matindex(igroup)) =                 &
                     CBAA_groupmat(kcount,gr_matindex(jgroup),gr_matindex(igroup))           & 
                     + (Leg_sign * exptable(ordercount) * hankeltable(ordercount) * Leg2)

!print *,'order =',ordercount
!print *,'Leg_sign =',Leg_sign
!print *,'Leg_x =',Leg_x
!print *,'Leg2 =',Leg2
!print *,'exptable(ordercount) =',exptable(ordercount)
!print *,'hankeltable(ordercount) =',hankeltable(ordercount)
!print *,'product =',(exptable(ordercount) * hankeltable(ordercount) * Leg2)


           END DO
         END DO
       END DO Grouploop2
     END DO Grouploop1

     DEALLOCATE(hankeltable)
     DEALLOCATE(exptable)

     ! Multiply with appropriate scaling factors:
     CBAA_elemvec_term1 = SQRT(2.0)* CBAA_elemvec_term1
     CBAA_elemvec_term2 = cj*SQRT(2.0)*k0* CBAA_elemvec_term2
     CBAA_groupmat = -cj*(k0/((4.0*PI)**2))* CBAA_groupmat

     ! Find max term of groupmat:
     kXmm = ABS(CBAA_groupmat(1,1,1))
     DO kcount = 1,num_k_dirs
       DO igroup = 1,gr_num
         DO jgroup = 1,gr_num
           IF (kXmm.LT.ABS(CBAA_groupmat(kcount,igroup,jgroup))) kXmm = ABS(CBAA_groupmat(kcount,igroup,jgroup))
         END DO
       END DO
     END DO
     print *,'max grgr element =',kXmm

   END SUBROUTINE CBAA_MAKE_FMM_MATRICES
!*******************************************************************************

   SUBROUTINE CBAA_FMM_DEBUGGER
     IMPLICIT NONE
!*******************************************************************************
! Temporary debugging of FMM matrix. Multiply VTVt out and compare with full BE
! matrix.
!*******************************************************************************

     COMPLEX(SPC), ALLOCATABLE, DIMENSION(:,:) :: temp_scalar_mat
     COMPLEX(SPC), ALLOCATABLE, DIMENSION(:,:,:) :: temp_vector_mat
     INTEGER(I4B) :: count,count2,kcount
     REAL(SP) :: errtot,maxerr
     INTEGER(I4B), DIMENSION(2) :: max_err_groups
     INTEGER(I4B) :: CBOUTPUT1 = 25                    ! Unit number for temporary output-file
     INTEGER(I4B) :: fmm_mem,bi_mem,fe_mem             ! Memory

     OPEN(UNIT=CBOUTPUT1,STATUS='REPLACE',FILE='fmm_err.out')

     ALLOCATE(temp_scalar_mat(gr_num,ap_dof))
     ALLOCATE(temp_vector_mat(gr_num,ap_dof,3))
     temp_BE_mat = (0.0,0.0)

     ! Add [Z'] contribution to <temp_BE_mat>:
     DO count = 1,ap_dof
       temp_BE_mat(count,CBAA_BE_colind(CBAA_BE_rowind(count):CBAA_BE_rowind(count+1)-1)) =   &
         temp_BE_mat(count,CBAA_BE_colind(CBAA_BE_rowind(count):CBAA_BE_rowind(count+1)-1))   &
         + CBAA_BE_val(CBAA_BE_rowind(count):CBAA_BE_rowind(count+1)-1)
     END DO

     ! Add the VTVt contribution to <temp_BE_mat>:
     DO kcount = 1,num_k_dirs

print *,'contribution from k-vector ',kcount,'/',num_k_dirs
    
       temp_scalar_mat = (0.0,0.0)
       temp_vector_mat = (0.0,0.0)

       ! Calculate [T][Vt]:
       DO count = 1,gr_num
         DO count2 = 1,gr_num

         temp_scalar_mat(count,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1)) =   &
           temp_scalar_mat(count,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1))   &
           + CBAA_groupmat(kcount,count,count2)*                                                   &
             CONJG(CBAA_elemvec_term1(kcount,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1)))

         temp_vector_mat(count,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1),1:2) =   &
           temp_vector_mat(count,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1),1:2)   &
           - CBAA_groupmat(kcount,count,count2)*                                                       &
             CONJG(CBAA_elemvec_term2(kcount,gredge_colind(gredge_rowind(count2):gredge_rowind(count2+1)-1),1:2))

         END DO
       END DO

       ! Add <kcount>-th VTVt to <temp_BE_mat>:
       DO count = 1,ap_dof

         temp_BE_mat(count,:) = temp_BE_mat(count,:)   &
           + k_dirs(kcount,4)*CBAA_elemvec_term1(kcount,count)*temp_scalar_mat(edge_grnum(count),:) 

         temp_BE_mat(count,:) = temp_BE_mat(count,:)   &
           + k_dirs(kcount,4)*(CBAA_elemvec_term2(kcount,count,1)*temp_vector_mat(edge_grnum(count),:,1))  &
           + k_dirs(kcount,4)*(CBAA_elemvec_term2(kcount,count,2)*temp_vector_mat(edge_grnum(count),:,2))

       END DO

     END DO


print *,'writing out error matrix ...'

     temp_BE_mat = temp_BE_mat - CBAA_BE_mat
     errtot = 0.0
     maxerr = 0.0
     DO count = 1,ap_dof
       write(CBOUTPUT1,*) '>> Row',count
       DO count2 = 1,ap_dof

         if (count.LT.11) then
           write(CBOUTPUT1,*) 'Error =',ABS(temp_BE_mat(count,count2))/ABS(CBAA_BE_mat(count,count2)), &
             '   ABS(ref) =',ABS(CBAA_BE_mat(count,count2)),'   groups = (',edge_grnum(count),',',edge_grnum(count2),')'
         end if
         
         if (maxerr.LT.ABS(temp_BE_mat(count,count2))/ABS(CBAA_BE_mat(count,count2))) then
           maxerr = ABS(temp_BE_mat(count,count2))/ABS(CBAA_BE_mat(count,count2))
           max_err_groups = (/ edge_grnum(count), edge_grnum(count2) /)
         end if

         errtot = errtot + ABS(temp_BE_mat(count,count2))/ABS(CBAA_BE_mat(count,count2))

       END DO
     END DO
     CLOSE(UNIT=CBOUTPUT1)

     print *,'FMM memory (bytes) = ',fmm_mem
     print *,'BI memory (bytes) = ',bi_mem
     print *,'average error =', errtot/REAL(ap_dof**2 - SIZE(CBAA_BE_val))
     print *,'max. error = ', maxerr
     print *,'max_err_groups =',max_err_groups

     STOP

     DEALLOCATE(temp_scalar_mat)
     DEALLOCATE(temp_vector_mat)

   END SUBROUTINE CBAA_FMM_DEBUGGER
!*******************************************************************************

END SUBROUTINE CBAA_SYSMAT

!*******************************************************************************
! External routines:
!*******************************************************************************


SUBROUTINE CBAA_COAX_ZIN_CALC
  USE coax_feed
  USE math_tools, ONLY: PHASE
  USE geometry
  USE matrix
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Calculates the input impedance at the coax feed, according to formulation
! in IEEE Trans. A&P, vol.43, p.1474 - Gong & Volakis.
! Created: May 2000 MMB
! 2001 Jan 21: Made into a interface routine. MMB
!*******************************************************************************
  INTEGER(I4B) :: iedge                         ! counters
  COMPLEX(SPC) :: dof_value,dof_value2
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)     ! The complex constant j
  REAL(SP) :: Z_c                               ! Characteristic impedance of the coax
  COMPLEX(SPC) :: V_plus,V_minus,V_tot,I_plus,I_minus,I_tot ! waves on the coax at the aperture
   
  DO AX_counter = 1,NUM_AX_cards ! cycle through all AX cards

    ! Initialise:
    AXdata(AX_counter)%coax_zin1  = (0.0,0.0) ! Z_in calculated as averages input impeadence
    AXdata(AX_counter)%coax_zin2  = (0.0,0.0) ! Z_in using average coax aperture field
    dof_value  = (0.0,0.0) ! edge field value
    dof_value2 = (0.0,0.0) ! average edge field value

    ! Find the average field in the coax:
    DO iedge = 1,AXdata(AX_counter)%coaxedges_num
      IF (SQRT(SUM((vertices(edges(AXdata(AX_counter)%coaxedges(iedge))%nodes(1))%coord- &
        AXdata(AX_counter)%coaxcentre)**2)).LT.EPS) THEN ! Elemental field points naturally away from centre
        dof_value = x_vec_c(renumbered_e1(AXdata(AX_counter)%coaxedges(iedge)))
      ELSE
        dof_value = - x_vec_c(renumbered_e1(AXdata(AX_counter)%coaxedges(iedge)))
      END IF
      AXdata(AX_counter)%coax_zin1 = AXdata(AX_counter)%coax_zin1 + &
        1/( 2*AXdata(AX_counter)%coax_Iabs*EXP(AXdata(AX_counter)%coax_Iphase*PI/180.0)/ &
        ((AXdata(AX_counter)%coax_b-AXdata(AX_counter)%coax_a)*dof_value) -              &
        2.0*PI*SQRT(AXdata(AX_counter)%coax_eps)/(Z_zero*LOG(AXdata(AX_counter)%coax_b/  &
        AXdata(AX_counter)%coax_a)) ) 

! Temporary test, all values printed should roughly be equal:
print *, dof_value


      dof_value2 = dof_value2 + dof_value
    END DO

    AXdata(AX_counter)%coax_zin1 = AXdata(AX_counter)%coax_zin1/AXdata(AX_counter)%coaxedges_num

    dof_value2 = dof_value2/AXdata(AX_counter)%coaxedges_num ! average coax aperture, radial E-field
    AXdata(AX_counter)%coax_zin2 = 1/( 2*AXdata(AX_counter)%coax_Iabs*EXP(AXdata(AX_counter)%coax_Iphase* &
      PI/180.0)/((AXdata(AX_counter)%coax_b-AXdata(AX_counter)%coax_a)*dof_value2) -                      &
      2.0*PI*SQRT(AXdata(AX_counter)%coax_eps)/(Z_zero*                                                   &
      LOG(AXdata(AX_counter)%coax_b/AXdata(AX_counter)%coax_a)) ) 

    ! Calculate wave values at the aperture:
    Z_c     = Z_zero*LOG(AXdata(AX_counter)%coax_b/AXdata(AX_counter)%coax_a)/  &
             ( 2.0*PI*SQRT(AXdata(AX_counter)%coax_eps) )
    V_tot   = dof_value2*(AXdata(AX_counter)%coax_b-AXdata(AX_counter)%coax_a)
    I_plus  = AXdata(AX_counter)%coax_Iabs*EXP(AXdata(AX_counter)%coax_Iphase*PI/180.0)
    V_plus  = Z_c*I_plus
    V_minus = V_tot - V_plus
    I_minus = - V_minus/Z_c
    I_tot   = I_plus + I_minus

    ! Write the result to the output file:
    WRITE (FILEOUT,'(//,20X,A,I4,/)') 'DATA OF THE COAXIAL CURRENT SOURCE NO.',0
    WRITE (FILEOUT,'(20X,4(A13))') 'real part','imag. part','magn.','phase'
    ! Write out the voltage:
    WRITE (FILEOUT,'(1X,A10,1X,A8,1(1X,E12.5))') 'Z_c','in Ohm', Z_c
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Voltage','in V', REAL(V_tot), &
    AIMAG(V_tot), ABS(V_tot), PHASE(V_tot)
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Current','in A', REAL(I_tot), &
    AIMAG(I_tot), ABS(I_tot), PHASE(I_tot)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'V_+','in V', REAL(V_plus),  AIMAG(V_plus)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'V_-','in V', REAL(V_minus), AIMAG(V_minus)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'I_+','in A', REAL(I_plus),  AIMAG(I_plus)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'I_-','in A', REAL(I_minus), AIMAG(I_minus)
    ! Write out the admittance:
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Admitt.','in A/V', REAL(1/AXdata(AX_counter)%coax_zin1), &
                                                 AIMAG(1/AXdata(AX_counter)%coax_zin1), &
                                                 ABS(1/AXdata(AX_counter)%coax_zin1),   &
                                                 PHASE(1/AXdata(AX_counter)%coax_zin1)
    ! Write out the impedance:
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Impedance','in Ohm', REAL(AXdata(AX_counter)%coax_zin1), &
                                                 AIMAG(AXdata(AX_counter)%coax_zin1), &
                                                 ABS(AXdata(AX_counter)%coax_zin1),   &
                                                 PHASE(AXdata(AX_counter)%coax_zin1)
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Impedance','in Ohm', REAL(AXdata(AX_counter)%coax_zin2), &
                                                 AIMAG(AXdata(AX_counter)%coax_zin2), &
                                                 ABS(AXdata(AX_counter)%coax_zin2),   &
                                                 PHASE(AXdata(AX_counter)%coax_zin2)
  END DO

END SUBROUTINE CBAA_COAX_ZIN_CALC
!*******************************************************************************


SUBROUTINE CBAA_COAX_POSTPRO
  USE coax_feed
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI
  USE geometry
  USE math_tools, ONLY: PHASE
  USE matrix
  USE nrtype
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Calculates the input impedance at the coax feed, according to formulation
! in MMB PH.D. thesis.
!
! 2002-02-10: Created. MMB
! 2002-05-09: Generalized to arbitrary element order. MMB.
!*******************************************************************************
  INTEGER(I4B) :: ielem,iface,iport,port_num,mm,dof_num ! counters
  REAL(SP),DIMENSION(ELEM_TRI_MATRIX_SIZE) :: mode_vbf_integrals

  ! Calculate incident voltage waves, characteristic impedances, 
  ! of every coax port. Initialise the total voltage values:
  DO iport = 1,NUM_AX_cards
    AXdata(iport)%Z_c    = LOG(AXdata(iport)%coax_b/AXdata(iport)%coax_a) *     &
	                       SQRT(AXdata(iport)%coax_mu/AXdata(iport)%coax_eps) * &
						   (Z_zero/(2.0*PI))
    AXdata(iport)%V_plus = AXdata(iport)%coax_Iabs *                   &
                           EXP(j*AXdata(iport)%coax_Iphase*PI/180.0) * &
                           AXdata(iport)%Z_c
    AXdata(iport)%V_tot  = (0.0,0.0)
  END DO
   
  ! Cycle trough all faces, element-wise to create the total voltage
  ! in AXdata()%V_minus, after which it is converted to the reflected amplitudes
  ! by subtracting <V_plus>:
  ELEMENT_LOOP: DO ielem = 1,num_elements
	FACE_LOOP: DO iface = 1,4

	  ! Cycle if this face is not part of a coax port:
	  IF (.NOT.faces(elements(ielem)%faces(iface))%coax_aperture) CYCLE FACE_LOOP

      ! Record the port number of the current port face:
      port_num = faces(elements(ielem)%faces(iface))%coaxnumber

      ! Integrate over the face the dot product of ^n X the face basis functions
	  ! with ^n X the coax mode:
	  CALL CBAA_MAKE_COAX_MODE_INTEGRAL(ielem,iface,AXdata(port_num)%normal, &
	                                    AXdata(port_num)%coaxcentre,         &
										AXdata(port_num)%coax_a,             &
										mode_vbf_integrals)

      ! Add to the running total for this port:
      DO mm = 1,ELEM_TRI_MATRIX_SIZE
	    dof_num = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,mm)
        IF (dof_num.GT.0) THEN
		  AXdata(port_num)%V_tot = AXdata(port_num)%V_tot + &
		                           x_vec_c(dof_num)*mode_vbf_integrals(mm)
		END IF
	  END DO
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

  ! Scale to obtain the total voltages, and then calculate the final,
  ! reflected amplitudes. Then write out to output file:
  DO iport = 1,NUM_AX_cards

    ! Total voltage:
	AXdata(iport)%V_tot = (1.0/(2.0*PI*AXdata(iport)%coax_a)) * &
	                      AXdata(iport)%V_tot

    ! Reflected voltage amplitude:
    AXdata(iport)%V_minus = AXdata(iport)%V_tot - AXdata(iport)%V_plus

    ! Input impedance:
    AXdata(iport)%coax_zin1 = AXdata(iport)%V_tot/(AXdata(iport)%V_plus/AXdata(iport)%Z_c - &
	                          AXdata(iport)%V_minus/AXdata(iport)%Z_c)
	
	! Write the result to the output file:
    WRITE (FILEOUT,'(//,20X,A,I4,/)') 'DATA OF THE COAXIAL CURRENT SOURCE NO.',0
    WRITE (FILEOUT,'(20X,4(A13))') 'real part','imag. part','magn.','phase'
    ! Write out the voltage:
    WRITE (FILEOUT,'(1X,A10,1X,A8,1(1X,E12.5))') 'Z_c','in Ohm', AXdata(iport)%Z_c
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Voltage','in V', REAL(AXdata(iport)%V_tot), &
    AIMAG(AXdata(iport)%V_tot), ABS(AXdata(iport)%V_tot), PHASE(AXdata(iport)%V_tot)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'V_+','in V', REAL(AXdata(iport)%V_plus),  AIMAG(AXdata(iport)%V_plus)
    WRITE (FILEOUT,'(1X,A10,1X,A8,2(1X,E12.5))') 'V_-','in V', REAL(AXdata(iport)%V_minus), AIMAG(AXdata(iport)%V_minus)
    ! Write out the reflection coefficient:
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Ref. coef.','in 1', REAL(AXdata(iport)%V_minus/AXdata(iport)%V_plus), &
                                                 AIMAG(AXdata(iport)%V_minus/AXdata(iport)%V_plus), &
                                                 ABS(AXdata(iport)%V_minus/AXdata(iport)%V_plus),   &
                                                 PHASE(AXdata(iport)%V_minus/AXdata(iport)%V_plus)
    ! Write out the admittance:
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Admitt.','in A/V', REAL(1.0/AXdata(iport)%coax_zin1), &
                                                 AIMAG(1/AXdata(iport)%coax_zin1), &
                                                 ABS(1/AXdata(iport)%coax_zin1),   &
                                                 PHASE(1/AXdata(iport)%coax_zin1)
    ! Write out the impedance:
    WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Impedance','in Ohm', REAL(AXdata(iport)%coax_zin1), &
                                                 AIMAG(AXdata(iport)%coax_zin1), &
                                                 ABS(AXdata(iport)%coax_zin1),   &
                                                 PHASE(AXdata(iport)%coax_zin1)
  END DO

END SUBROUTINE CBAA_COAX_POSTPRO
!*******************************************************************************


SUBROUTINE CBAA_FIELDCALC(xob,yob,zob,E_xyz,H_xyz)
  USE CBAA_data
  USE feminterface, ONLY: FEM_FIELDCALC, &
                          LOCAL_TO_GLOBAL_INDEX_TRI
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! Calculates E and H field in the CBAA_ANALYSIS case, at the location specified.
! 19 Jan 2001 - MMB
!*******************************************************************************
  REAL(SP), INTENT(IN) :: xob,yob,zob
  COMPLEX(SPC), DIMENSION(3), INTENT(OUT) :: E_xyz,H_xyz

  INTEGER(I4B), DIMENSION(3) :: tempedges 
  REAL(SP), DIMENSION(25,3) :: gauss_xyz            ! xyz co-ord.'s of triangular quadrature points
  REAL(SP), DIMENSION(25,8,3) :: gauss_t2           ! Values of ^z cross basis func.s at <gauss_xyz>
  LOGICAL(LGT) :: test_free1
  COMPLEX(SPC), DIMENSION(3) :: F_ap                ! Aperture value of ^z cross (basis func), weighted with facet area and gauss weight
  REAL(SP) :: dxob,dyob,dzob,sqrra
  COMPLEX(SPC) :: Gre0,gterm1,gterm2,gfactor1
  COMPLEX(SPC) :: ddxGre0,ddyGre0,ddzGre0,d2dx2Gre0,d2dy2Gre0,d2dxdyGre0,  &
                  d2dxdzGre0,d2dydzGre0             ! Derivatives of Green function
  INTEGER(I4B) :: inside_el                         ! element number in which the near field point is 
  INTEGER(I4B) :: ielem,row,ifunc                   ! counters
  REAL(SP), DIMENSION(20,3) :: func_values
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)         ! The complex constant j
  INTEGER(I4B) :: jj                                ! counters
  INTEGER(I4B) :: num_qpoints, qcount               ! number of quadrature points for triangular surface integration
  INTEGER(I4B) :: face_hierarchal_order
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE)   :: temp_quad_term1
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE,3) :: temp_quad_term2

  ! Choose the number of quadrature points:
  num_qpoints = 7

  ! Initialize:
  H_xyz = (0.0,0.0) 
  E_xyz = (0.0,0.0)

  z_choice: IF (zob.GE.0.0) THEN

    Element_integration_loop: DO ielem = 1,num_apelements

      ! Check if the facet does contain any dofs, and cycle to next aperture 
      ! element if this one is constrained:
      tempedges = which_local_edges(ielem,1:3)
      test_free1 = edges(elements(ap_elnumbers(ielem))%edges(tempedges(1)))%free.OR. &
                   edges(elements(ap_elnumbers(ielem))%edges(tempedges(2)))%free.OR. &
                   edges(elements(ap_elnumbers(ielem))%edges(tempedges(3)))%free.OR. &
                   faces(elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))%free
      IF (.NOT.test_free1) CYCLE Element_integration_loop

      ! Choose the dimension of the elemental contributions (number af face functions):
      face_hierarchal_order = MAX_ORDER(ap_elnumbers(ielem),which_local_face(ielem))

      ! Calculate the xyz-quadrature points and the dof values at these points:
      SELECT CASE (face_hierarchal_order)
      CASE (1) ! CT/LN
        CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,3, &
		                                  gauss_xyz(1:num_qpoints,1:3), &
		                                  temp_quad_term1(1:num_qpoints,1:3),temp_quad_term2(1:num_qpoints,1:3,1:3))
        gauss_t2(1:num_qpoints,1:3,1:3) = temp_quad_term2(1:num_qpoints,1:3,1:3)
      CASE (2) ! LT/QN
        CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,8, &
		                                  gauss_xyz(1:num_qpoints,1:3), &
		                                  temp_quad_term1(1:num_qpoints,1:8),temp_quad_term2(1:num_qpoints,1:8,1:3))
        gauss_t2(1:num_qpoints,1:8,1:3) = temp_quad_term2(1:num_qpoints,1:8,1:3)
      CASE DEFAULT
        STOP 'IE: Invalid hierarchal order in CBAA_FIELDCALC.'
      END SELECT
       
      ! Now perform the integral over the face area:
      Gauss_rad_integration: DO qcount = 1,num_qpoints
        ! Calculate Green function related values for this quadrature point:
        dxob = xob - gauss_xyz(qcount,1)
        dyob = yob - gauss_xyz(qcount,2)
        dzob = zob - gauss_xyz(qcount,3)
        sqrra = SQRT(dxob**2 + dyob**2 + dzob**2) 
        Gre0 = EXP(-cj*k0*sqrra)/(4*pi*sqrra)
        gfactor1 = (cj*k0/sqrra) + 1/(sqrra**2)
        gterm1 = (cj*k0/(sqrra**3)) + 2.0/(sqrra**4)
        gterm2 = (3*cj*k0/(sqrra**5)) + 8.0/(sqrra**6)
        ! derivatives:
        ddxGre0      = - Gre0*dxob*gfactor1
        ddyGre0      = - Gre0*dyob*gfactor1
        ddzGre0      = - Gre0*dzob*gfactor1
        d2dx2Gre0    = Gre0*( (dxob**2)*(gfactor1**2 + gterm1) - gfactor1 )
        d2dy2Gre0    = Gre0*( (dyob**2)*(gfactor1**2 + gterm1) - gfactor1 )
        d2dxdyGre0   = Gre0*dxob*dyob*( gfactor1**2 + gterm1 )
        d2dxdzGre0   = Gre0*dxob*dzob*( gfactor1**2 + gterm1 )
        d2dydzGre0   = Gre0*dyob*dzob*( gfactor1**2 + gterm1 )
        ! Add the contribution of every basis function at this quadrature point to the grand total:
        Function_loop: DO jj = 1,8
          row = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(ielem),which_local_face(ielem),jj)
          IF (row.GT.0) THEN ! check that this is a dof (free variable)
            F_ap = x_vec_c(row)*gauss_t2(qcount,jj,1:3)
            H_xyz(1) = H_xyz(1) + F_ap(1)*(Gre0 + d2dx2Gre0/(k0**2)) + F_ap(2)*d2dxdyGre0/(k0**2)
            H_xyz(2) = H_xyz(2) + F_ap(2)*(Gre0 + d2dy2Gre0/(k0**2)) + F_ap(1)*d2dxdyGre0/(k0**2) 
            H_xyz(3) = H_xyz(3) + F_ap(1)*d2dxdzGre0/(k0**2) + F_ap(2)*d2dydzGre0/(k0**2)
            E_xyz(1) = E_xyz(1) - F_ap(2)*ddzGre0
            E_xyz(2) = E_xyz(2) + F_ap(1)*ddzGre0
            E_xyz(3) = E_xyz(3) + F_ap(2)*ddxGre0 - F_ap(1)*ddyGre0 
          END IF
        END DO Function_loop
      END DO Gauss_rad_integration
    END DO Element_integration_loop 

    ! Multiply by factors standing outside the integrals:
    H_xyz = -cj*2.0*(k0/Z_zero)*H_xyz
    E_xyz = -2.0*E_xyz
  ELSE ! z_choice
    CALL FEM_FIELDCALC(xob,yob,zob,E_xyz,H_xyz)
  END IF z_choice

END SUBROUTINE CBAA_FIELDCALC
!*******************************************************************************


SUBROUTINE CBAA_FIELDCALC_BI(r_vec,E_xyz,H_xyz,elem_special,num_qpts_special)
  USE CBAA_data
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! Calculates E and H field in the CBAA_ANALYSIS case, at the location specified.
! With the boundary integral.
!
! The special element and number of quadrature points signifies the following:
! -> If num_qpts_special=0: the contribution of <elem_special> is ignored
! -> If num_qpts_special>0: the contribution of <elem_special> is integrated with
!                           <num_qpts_special> quadrature points.
!
! 2002-05-14: Created from CABB_FIELDCALC. MMB.
! 2002-05-16: Added the special option. MMB.
!*******************************************************************************
  REAL(SP), DIMENSION(3), INTENT(IN) :: r_vec
  COMPLEX(SPC), DIMENSION(3), INTENT(OUT) :: E_xyz,H_xyz
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: elem_special,num_qpts_special

  INTEGER(I4B) :: ielem,row,ifunc                   ! counters
  INTEGER(I4B) :: num_qpoints, qcount               ! number of quadrature points for triangular surface integration
  INTEGER(I4B), DIMENSION(3) :: tempedges 
  LOGICAL(LGT) :: test_free1,special_element
  REAL(SP) :: dxob,dyob,dzob,sqrra
  COMPLEX(SPC) :: Gre0,gterm1,gterm2
  COMPLEX(SPC) :: ddxGre0,ddyGre0,ddzGre0,d2dx2Gre0,d2dy2Gre0,d2dz2Gre0,d2dxdyGre0, &
                  d2dxdzGre0,d2dydzGre0             ! Derivatives of Green function
  COMPLEX(SPC), DIMENSION(3,3) :: delaccent_times_IGre0,dyadic_Gre0
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,3)                      :: quad_xyz
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE)   :: quad_term1
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE,3) :: quad_term2
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)         ! The complex constant j

  ! Check if special element is requested:
  special_element = .FALSE.
  IF (PRESENT(elem_special).AND.PRESENT(num_qpts_special)) special_element = .TRUE.

  ! Choose the number of quadrature points:
  num_qpoints = 7 ! (note that this is also set elsewhere in this routine)

  ! Initialize:
  H_xyz = (0.0,0.0) 
  E_xyz = (0.0,0.0)

  Element_integration_loop: DO ielem = 1,num_apelements

    ! Check the special element:
	IF (special_element) THEN
	  IF ((ap_elnumbers(ielem).EQ.elem_special).AND.(num_qpts_special.EQ.0)) CYCLE Element_integration_loop
	  IF ((ap_elnumbers(ielem).EQ.elem_special).AND.(num_qpts_special.GT.0)) num_qpoints = num_qpts_special
    END IF    

    ! Check if the facet does contain any dofs, and cycle to next aperture 
    ! element if this one is constrained:
    tempedges = which_local_edges(ielem,1:3)
    test_free1 = edges(elements(ap_elnumbers(ielem))%edges(tempedges(1)))%free.OR. &
                 edges(elements(ap_elnumbers(ielem))%edges(tempedges(2)))%free.OR. &
                 edges(elements(ap_elnumbers(ielem))%edges(tempedges(3)))%free.OR. &
                 faces(elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))%free
    IF (.NOT.test_free1) CYCLE Element_integration_loop

    ! Calculate the xyz-quadrature points and the dof values at these points:
    CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,ELEM_TRI_MATRIX_SIZE, &
		                              quad_xyz(1:num_qpoints,1:3), &
		                              quad_term1(1:num_qpoints,1:ELEM_TRI_MATRIX_SIZE), &
									  quad_term2(1:num_qpoints,1:ELEM_TRI_MATRIX_SIZE,1:3))
       
    ! Now perform the integral over the face area:
    QUAD_rad_integration: DO qcount = 1,num_qpoints

      ! Calculate Green function related values for this quadrature point:
      dxob       = r_vec(1) - quad_xyz(qcount,1)
      dyob       = r_vec(2) - quad_xyz(qcount,2)
      dzob       = r_vec(3) - quad_xyz(qcount,3)
      sqrra      = SQRT(dxob**2 + dyob**2 + dzob**2) 
      Gre0       = EXP(-cj*k0*sqrra)/(4*PI*sqrra)
      gterm1     = (cj*k0/sqrra) + 1/(sqrra**2)
      gterm2     = (cj*k0/(sqrra**3)) + 2.0/(sqrra**4)
      ddxGre0    = - dxob * gterm1 * Gre0
      ddyGre0    = - dyob * gterm1 * Gre0
      ddzGre0    = - dzob * gterm1 * Gre0
      d2dx2Gre0  = (((dxob**2)*(gterm1**2 + gterm2)) - gterm1) * Gre0
      d2dy2Gre0  = (((dyob**2)*(gterm1**2 + gterm2)) - gterm1) * Gre0
      d2dz2Gre0  = (((dzob**2)*(gterm1**2 + gterm2)) - gterm1) * Gre0
      d2dxdyGre0 = dxob * dyob * (gterm1**2 + gterm2) * Gre0
      d2dxdzGre0 = dxob * dzob * (gterm1**2 + gterm2) * Gre0
      d2dydzGre0 = dyob * dzob * (gterm1**2 + gterm2) * Gre0

      dyadic_Gre0(1,1) = Gre0 + (1.0/(k0**2)) * d2dx2Gre0
      dyadic_Gre0(2,2) = Gre0 + (1.0/(k0**2)) * d2dy2Gre0
      dyadic_Gre0(3,3) = Gre0 + (1.0/(k0**2)) * d2dz2Gre0
      dyadic_Gre0(2,1) = (1.0/(k0**2)) * d2dxdyGre0
      dyadic_Gre0(3,1) = (1.0/(k0**2)) * d2dxdzGre0
      dyadic_Gre0(3,2) = (1.0/(k0**2)) * d2dydzGre0
	  dyadic_Gre0(1,2) = dyadic_Gre0(2,1)
      dyadic_Gre0(1,3) = dyadic_Gre0(3,1)
      dyadic_Gre0(2,3) = dyadic_Gre0(3,2)

      delaccent_times_IGre0(1,1) = (0.0,0.0)
      delaccent_times_IGre0(2,2) = (0.0,0.0)
      delaccent_times_IGre0(3,3) = (0.0,0.0)
      delaccent_times_IGre0(2,1) = - ddzGre0
      delaccent_times_IGre0(3,1) =   ddyGre0
      delaccent_times_IGre0(3,2) = - ddxGre0
      delaccent_times_IGre0(1,2) =   ddzGre0
      delaccent_times_IGre0(1,3) = - ddyGre0
      delaccent_times_IGre0(2,3) =   ddxGre0

      ! Add the contribution of every basis function at this 
	  ! quadrature point to totals:
      Function_loop: DO ifunc = 1,ELEM_TRI_MATRIX_SIZE
        row = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(ielem),which_local_face(ielem),ifunc)
        IF (row.GT.0) THEN ! check that this is a dof (free variable)
          H_xyz = H_xyz + x_vec_c(row)*MATMUL(quad_term2(qcount,ifunc,1:3),dyadic_Gre0)
          E_xyz = E_xyz + x_vec_c(row)*MATMUL(quad_term2(qcount,ifunc,1:3),delaccent_times_IGre0)
        END IF
      END DO Function_loop

    END DO QUAD_rad_integration

    ! Check the special element:
	IF (special_element) THEN
	  IF ((ap_elnumbers(ielem).EQ.elem_special).AND.(num_qpts_special.GT.0)) num_qpoints = 7
    END IF    

  END DO Element_integration_loop 

  ! Multiply by factors standing outside the integrals:
  H_xyz = cj*2.0*(k0/Z_zero)*H_xyz
  E_xyz = 2.0*E_xyz

END SUBROUTINE CBAA_FIELDCALC_BI
!*******************************************************************************


SUBROUTINE MAKE_ERM_CBAA_BI_LINEAR(elem,local_face,num_vbfs,num_qpts,BI_linear,elem_special,num_qpts_special)
  USE basis_function, ONLY: POINT_EVALUATE_FACE_FUNCTIONS
  USE frequency_data
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! This routine creates the elemental contribution to the ERM local problem's
! linear form from the CBAA BI of the previous approximate field.
!
! <num_vbfs> has no effect currently. The full vector is always returned.
!
! 2002-05-16: Created. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face,num_vbfs,num_qpts
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE), INTENT(OUT) :: BI_linear
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: elem_special,num_qpts_special
  
  INTEGER(I4B) :: iquad,ifunc
  REAL(SP), DIMENSION(4,4) :: vmat
  REAL(SP), DIMENSION(4) :: simp,cart
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes
  COMPLEX(SPC), DIMENSION(3) :: E_field,H_field
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE) :: term1
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE,3) :: term2

  ! Calculate the vertices matrix:
  CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vmat)

  ! Calculate for simplex assignment:
  tempfacenodes = LOCAL_FACENODES(local_face)

  simp = 0.0 ! init

  BI_linear = (0.0,0.0) ! init

  QUAD_LOOP: DO iquad = 1,num_qpts
  
    ! Calculate the xyz position of the quadrature point:
    simp(tempfacenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
    cart = MATMUL(vmat,simp)
  
    ! Evaluate the approximate field:
    STOP 'Code in MAKE_ERM_CBAA_BI_LINEAR still needs fixing... DBD 8 Feb 05'
!    Following line edited temporarily during porting to SGI.
!    CALL CBAA_FIELDCALC_BI(cart(1:3),E_field,H_field,elem,0)    
    
	! Scale and carry out double (^z \times):
    H_field = -j*k0*Z_zero*H_field
	H_field(3) = (0.0,0.0)
    H_field = - H_field

    ! Evaluate the basis functions at the current quadrature point:
    CALL POINT_EVALUATE_FACE_FUNCTIONS(elem,local_face,cart(1:3),ELEM_TRI_MATRIX_SIZE,term1,term2)

    ! Add to total:
    FUNC_LOOP: DO ifunc = 1,ELEM_TRI_MATRIX_SIZE
      BI_linear(ifunc) = quad_tri_rules(num_qpts)%rule(iquad,4)*SUM(H_field*term2(ifunc,1:3))
	END DO FUNC_LOOP
  END DO QUAD_LOOP

  ! Scale by area (last step in quadrature):
  BI_linear = FACE_AREA(elem,local_face)*BI_linear

END SUBROUTINE MAKE_ERM_CBAA_BI_LINEAR
!*******************************************************************************


SUBROUTINE CBAA_FMM_ALLOCATE
  USE CBAA_data
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! Allocates the FMM matrix storage and administrative structures. The administrative 
! structures are also filled. Later this subroutine may also iterate to find the 
! optimal gr_dim.
! Created in August 2000 MMB
!*******************************************************************************

  INTEGER(I4B) iedge,ielem,iap_dof                          ! counters
  INTEGER(I4B) :: gl_edgenum,gl_facenum,gl_groupnum         ! Global edge/face/group number
  INTEGER(I4B) :: gr_numy                                   ! number of ap_dof-occupied groups
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: row_dofs       ! temp store of colind of a row
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: temp_colind    ! temp store of CBAA_BE_colind
  INTEGER(I4B) :: igroup,grindex,colcount,tempint,tempint2  ! counters
  INTEGER(I4B) :: alloc_val,val_tot                         ! counters
  INTEGER(I4B) :: this_group                                ! temp
  REAL(SP) :: kdir_phi,kdir_theta,dkdir_phi,dkdir_theta     ! used for kdir specification
  INTEGER(I4B) :: kdir_ppoints,kdir_tpoints,kdir_torder     ! number of points in theta/phi directions, order of G.-Leg. rule
  INTEGER(I4B) :: k_integration_scheme                      ! temp: flag to choose the spherical integration scheme
                                                            ! (phi-theta) 1 = trap-trap, 2 = trap-gauss_Legendre
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: Gauss_Legendre ! Table of Gauss-Legendre points and weights
  REAL(SP) :: L_errlevel,L_k,L_alpha,L_m,L_c                ! Calculation of L_max (see MMB:'Electromagnetics')

  ! 30 rules with 1..30 points respectively, (30,30,1) = points, (30,30,2) = weights.
  ALLOCATE(Gauss_Legendre(40,40,2)) 
  Gauss_Legendre = 0.0
  Gauss_Legendre(10,1:10,1) = (/  -0.97390652851717,  -0.86506336668899,  -0.67940956829902, &
       -0.43339539412925,  -0.14887433898163,   0.14887433898163,   0.43339539412925,   0.67940956829902, &
          0.86506336668899,   0.97390652851717 /) 
  Gauss_Legendre(10,1:10,2) = (/   0.06667134430869,   0.14945134915058,   0.21908636251599, &
        0.26926671930998,   0.29552422471476,   0.29552422471475,   0.26926671930999,   0.21908636251599, &
           0.14945134915058,   0.06667134430869 /)

  Gauss_Legendre(11,1:11,1) = (/  -0.97822865814605,  -0.88706259976810,  -0.73015200557404, &
       -0.51909612920681,  -0.26954315595234,                  0.0,   0.26954315595234,   0.51909612920681, &
          0.73015200557404,   0.88706259976810,   0.97822865814605/)
  Gauss_Legendre(11,1:11,2) = (/   0.05566856711618,   0.12558036946489,   0.18629021092776, &
        0.23319376459196,   0.26280454451027,   0.27292508677789,   0.26280454451024,   0.23319376459199, &
           0.18629021092774,   0.12558036946490,   0.05566856711617/)

  Gauss_Legendre(12,1:12,1) = (/  -0.98156063424672,  -0.90411725637048,  -0.76990267419430, &
       -0.58731795428662,  -0.36783149899818,  -0.12523340851147,   0.12523340851147,   0.36783149899818, &
          0.58731795428662,   0.76990267419430,   0.90411725637048,   0.98156063424672 /)
  Gauss_Legendre(12,1:12,2) = (/   0.04717533638651,   0.10693932599533,   0.16007832854333, &
        0.20316742672309,   0.23349253653832,   0.24914704581345,   0.24914704581334,   0.23349253653841, &
           0.20316742672303,   0.16007832854336,   0.10693932599532,   0.04717533638651 /)

  Gauss_Legendre(13,1:13,1) = (/  -0.98418305471859,  -0.91759839922299,  -0.80157809073328, &
       -0.64234933944036,  -0.44849275103644,  -0.23045831595513,                  0.0,   0.23045831595513, &
          0.44849275103644,   0.64234933944036,   0.80157809073328,   0.91759839922299,   0.98418305471859 /)
  Gauss_Legendre(13,1:13,2) = (/   0.04048400476529,   0.09212149983781,   0.13887351021968, &
        0.17814598076209,   0.20781604753664,   0.22628318026323,   0.23255155323052,   0.22628318026319, &
           0.20781604753670,   0.17814598076205,   0.13887351021971,   0.09212149983779,   0.04048400476529 /)

  Gauss_Legendre(14,1:14,1) = (/  -0.98628380869689,  -0.92843488366343,  -0.82720131506986, &
       -0.68729290481164,  -0.51524863635817,  -0.31911236892789,  -0.10805494870734,   0.10805494870734, &
          0.31911236892789,   0.51524863635817,   0.68729290481164,   0.82720131506986,   0.92843488366343, &
             0.98628380869689 /)
  Gauss_Legendre(14,1:14,2) = (/   0.03511946033169,   0.08015808716003,   0.12151857068738, &
        0.15720316715882,   0.18553839747734,   0.20519846372174,   0.21526385346298,   0.21526385346301, &
           0.20519846372173,   0.18553839747732,   0.15720316715885,   0.12151857068736,   0.08015808716005, &
              0.03511946033169 /)

  Gauss_Legendre(15,1:15,1) = (/  -0.98799251802043,  -0.93727339240083,  -0.84820658341032, &
       -0.72441773136022,  -0.57097217260853,  -0.39415134707757,  -0.20119409399743,                  0.0, &
           0.20119409399743,   0.39415134707757,   0.57097217260853,   0.72441773136022,   0.84820658341032, &
              0.93727339240083,   0.98799251802043 /)
  Gauss_Legendre(15,1:15,2) = (/   0.03075324199612,   0.07036604748802,   0.10715922046739, &
        0.13957067792591,   0.16626920581715,   0.18616100001556,   0.19843148532693,   0.20257824192593, &
           0.19843148532665,   0.18616100001596,   0.16626920581676,   0.13957067792620,   0.10715922046721, &
             0.07036604748811,   0.03075324199609 /)

  Gauss_Legendre(16,1:16,1) = (/  -0.98940093499140,  -0.94457502307367,  -0.86563120238763, &
       -0.75540440835496,  -0.61787624440270,  -0.45801677765721,  -0.28160355077926,  -0.09501250983764, &
           0.09501250983764,   0.28160355077926,   0.45801677765721,   0.61787624440270,   0.75540440835496, &
             0.86563120238763,   0.94457502307367,   0.98940093499140 /)
  Gauss_Legendre(16,1:16,2) = (/   0.02715245941199,   0.06225352393758,   0.09515851168441, &
        0.12462897125337,   0.14959598881866,   0.16915651939313,   0.18260341504623,   0.18945061045471, &
           0.18945061045441,   0.18260341504633,   0.16915651939325,   0.14959598881847,   0.12462897125355, &
              0.09515851168429,   0.06225352393765,   0.02715245941197 /)

  Gauss_Legendre(17,1:17,1) = (/  -0.99057547531407,  -0.95067552176944,  -0.88023915372651, &
       -0.78151400389704,  -0.65767115921655,  -0.51269053708654,  -0.35123176345386,  -0.17848418149585, &
                         0.0,   0.17848418149585,   0.35123176345386,   0.51269053708654,   0.65767115921655, &
                            0.78151400389704,   0.88023915372651,   0.95067552176944,   0.99057547531407 /)
  Gauss_Legendre(17,1:17,2) = (/   0.02414830286864,   0.05545952937329,   0.08503614831855, &
        0.11188384719197,   0.13513636846971,   0.15404576107597,   0.16800410215688,   0.17656270536693, &
           0.17944647035612,   0.17656270536693,   0.16800410215688,   0.15404576107594,   0.13513636846978, &
              0.11188384719187,   0.08503614831864,   0.05545952937324,   0.02414830286866 /)
    
  Gauss_Legendre(18,1:18,1) = (/  -0.99156516841999,  -0.95582394957356,  -0.89260246649552, &
       -0.80370495897374,  -0.69168704305977,  -0.55977083107416,  -0.41175116146280,  -0.25188622569151, &
         -0.08477501304174,   0.08477501304174,   0.25188622569151,   0.41175116146280,   0.55977083107416, &
            0.69168704305977,   0.80370495897374,   0.89260246649552,   0.95582394957356,   0.99156516841999 /)
  Gauss_Legendre(18,1:18,2) = (/   0.02161601352662,   0.04971454889301,   0.07642573025985, &
        0.10094204409975,   0.12255520671787,   0.14064291466542,   0.15468467512988,   0.16427648374387, &
           0.16914238296361,   0.16914238296412,   0.16427648374328,   0.15468467513054,   0.14064291466474, &
              0.12255520671849,   0.10094204409925,   0.07642573026019,   0.04971454889283,   0.02161601352668 /)
    
  Gauss_Legendre(19,1:19,1) = (/  -0.99240684384518,  -0.96020815213095,  -0.90315590361877, &
       -0.82271465653477,  -0.72096617733613,  -0.60054530466145,  -0.46457074137601,  -0.31656409996362, &
         -0.16035864564023,                  0.0,   0.16035864564023,   0.31656409996362,   0.46457074137601, &
            0.60054530466145,   0.72096617733613,   0.82271465653477,   0.90315590361877,   0.96020815213095, &
               0.99240684384518 /)
  Gauss_Legendre(19,1:19,2) = (/   0.01946178822953,   0.04481422676922,   0.06904454272783, &
        0.09149002163678,   0.11156664553150,   0.12875396255515,   0.14260670215809,   0.15276604208098, &
           0.15896884337922,   0.16105444986356,   0.15896884337866,   0.15276604208171,   0.14260670215759, &
              0.12875396255528,   0.11156664553166,   0.09149002163650,   0.06904454272810,   0.04481422676906, &
                 0.01946178822959 /)
    
  Gauss_Legendre(20,1:20,1) = (/  -0.99312859918642,  -0.96397192727432,  -0.91223442825573, &
       -0.83911697181872,  -0.74633190646216,  -0.63605368072565,  -0.51086700195110,  -0.37370608871536, &
         -0.22778585114165,  -0.07652652113350,   0.07652652113350,   0.22778585114165,   0.37370608871536, &
            0.51086700195110,   0.63605368072565,   0.74633190646216,   0.83911697181872,   0.91223442825573, &
               0.96397192727432,   0.99312859918642 /)
  Gauss_Legendre(20,1:20,2) = (/   0.01761400714033,   0.04060142979901,   0.06267204833313, &
        0.08327674157899,   0.10193011981703,   0.11819453195636,   0.13168863846167,   0.14209610929839, &
           0.14917298649815,   0.15275338710321,   0.15275338715642,   0.14917298645102,   0.14209610933570, &
              0.13168863843485,   0.11819453197424,   0.10193011980576,   0.08327674158573,   0.06267204832939, &
                 0.04060142980075,   0.01761400713985 /)
    
  Gauss_Legendre(21,1:21,1) = (/  -0.99375217062046,  -0.96722683856476,  -0.92009933415443, &
       -0.85336336457845,  -0.76843996347942,  -0.66713880419534,  -0.55161883588804,  -0.42434212020722, &
         -0.28802131680243,  -0.14556185416089,                0.0,   0.14556185416089,   0.28802131680243, &
            0.42434212020722,   0.55161883588804,   0.66713880419534,   0.76843996347943,   0.85336336457845, &
               0.92009933415443,   0.96722683856476,   0.99375217062046 /)
  Gauss_Legendre(21,1:21,2) = (/   0.01601722825505,   0.03695378978228,   0.05713442539871, &
        0.07610011368087,   0.09344442337534,   0.10879729927822,   0.12183141590949,   0.13226893881298, &
           0.13988739457831,   0.14452440422605,   0.14608113340636,   0.14452440422304,   0.13988739458215, &
              0.13226893881045,   0.12183141591029,   0.10879729927830,   0.09344442337523,   0.07610011368069, &
                 0.05713442539905,   0.03695378978200,   0.01601722825515 /)

  Gauss_Legendre(22,1:22,1) = (/  -0.99429458549330,  -0.97006049780519,  -0.92695677222523, &
       -0.86581257768978,  -0.78781680599617,  -0.69448726318007,  -0.58764040350850,  -0.46935583798667, &
         -0.34193582089203,  -0.20786042668823,  -0.06973927331972,   0.06973927331972,   0.20786042668823, &
            0.34193582089203,   0.46935583798667,   0.58764040350850,   0.69448726318007,   0.78781680599617, &
               0.86581257768978,   0.92695677222523,   0.97006049780519,   0.99429458549330 /)
  Gauss_Legendre(22,1:22,2) = (/   0.01462799530009,   0.03377490160194,   0.05229333508555, &
        0.06979646854167,   0.08594160607107,   0.10041414459807,   0.11293229592496,   0.12325237696366, &
           0.13117350464054,   0.13654149847750,   0.13925187274978,   0.13925187292673,   0.13654149831470, &
              0.13117350477871,   0.12325237685448,   0.11293229600659,   0.10041414453937,   0.08594160611187, &
                 0.06979646851466,   0.05229333510187,   0.03377490159384,   0.01462799530237 /)

  Gauss_Legendre(23,1:23,1) = (/  -0.99476933496836,  -0.97254247128854,  -0.93297108675496, &
       -0.87675235831393,  -0.80488840160091,  -0.71866136313725,  -0.61960987576236,  -0.50950147784630, &
         -0.39030103803026,  -0.26413568097034,  -0.13325682429847,                0.0,   0.13325682429847, &
            0.26413568097034,   0.39030103803026,   0.50950147784630,   0.61960987576236,   0.71866136313725, &
               0.80488840160091,   0.87675235831393,   0.93297108675496,   0.97254247128854,   0.99476933496836 /)
  Gauss_Legendre(23,1:23,2) = (/   0.01341185946023,   0.03098800589998,   0.04803767169786, &
        0.06423242148180,   0.07928141160174,   0.09291576636091,   0.10489209105171,   0.11499664070985, &
          0.12304908379504,   0.12890572267060,   0.13246203899391,   0.13365457250519,   0.13246203917365, &
             0.12890572234809,   0.12304908420113,   0.11499664028128,   0.10489209145472,   0.09291576601374, &
                0.07928141187809,   0.06423242128041,   0.04803767182755,   0.03098800583300,   0.01341185947950 /)

  Gauss_Legendre(24,1:24,1) = (/  -0.99518722006497,  -0.97472855580952,  -0.93827455216393, &
       -0.88641552690256,  -0.82000198602422,  -0.74012419155605,  -0.64809365194580,  -0.54542147138625, &
         -0.43379350762650,  -0.31504267969615,  -0.19111886747361,  -0.06405689286261,   0.06405689286261, &
            0.19111886747361,   0.31504267969615,   0.43379350762650,   0.54542147138625,   0.64809365194580, &
               0.74012419155605,   0.82000198602422,   0.88641552690256,   0.93827455216393,   0.97472855580952, &
                  0.99518722006497 /)
  Gauss_Legendre(24,1:24,2) = (/   0.01234122984242,   0.02853138859521,   0.04427743877673, &
        0.05929858491871,   0.07334648158193,   0.08619016111680,   0.09761865278321,   0.10744426917873, &
           0.11550566922221,   0.12167047157564,   0.12583745781807,   0.12793819383296,   0.12793819682124, &
              0.12583745498676,   0.12167047412151,   0.11550566704372,   0.10744427095864,   0.09761865139033, &
                 0.08619016216160,   0.07334648083451,   0.05929858542047,   0.04427743847263,   0.02853138874597, &
                    0.01234122980001 /)

  Gauss_Legendre(25,1:25,1) = (/  -0.99555697015776,  -0.97666392057791,  -0.94297457208622, &
       -0.89499199742080,  -0.83344262887921,  -0.75925926304811,  -0.67356636845102,  -0.57766293025168, &
         -0.47300273144262,  -0.36117230580995,  -0.24386688372095,  -0.12286469261071,                0.0, &
            0.12286469261071,   0.24386688372095,   0.36117230580995,   0.47300273144262,   0.57766293025168, &
               0.67356636845102,   0.75925926304811,   0.83344262887921,   0.89499199742080,   0.94297457208622, &
                  0.97666392057791,   0.99555697015776 /)
  Gauss_Legendre(25,1:25,2) = (/   0.01139379860708,   0.02635498688799,   0.04093915556615, &
        0.05490469735730,   0.06803833302528,   0.08014070014094,   0.09102826316914,   0.10053594696464, &
           0.10851962742842,   0.11485825544210,   0.11945576778392,   0.12224243852009,   0.12317605803348, &
              0.12224243919823,   0.11945576657806,   0.11485825694876,   0.10851962583750,   0.10053594848044, &
                 0.09102826182688,   0.08014070125777,   0.06803833215386,   0.05490469798503,   0.04093915516459, &
                    0.02635498709465,   0.01139379854772 /)

  Gauss_Legendre(26,1:26,1) = (/  -0.99588570073304,  -0.97838544702884,  -0.94715906543980, &
       -0.90263786285020,  -0.84544594237955,  -0.77638594894672,  -0.69642726039801,  -0.60669229301824, &
         -0.50844071482503,  -0.40305175512345,  -0.29200483948593,  -0.17685882035690,  -0.05923009342931, &
            0.05923009342931,   0.17685882035690,   0.29200483948593,   0.40305175512345,   0.50844071482503, &
               0.60669229301824,   0.69642726039801,   0.77638594894672,   0.84544594237955,   0.90263786285020, &
                  0.94715906543980,   0.97838544702884,   0.99588570073304 /)
  Gauss_Legendre(26,1:26,2) = (/   0.01055137247535,   0.02441785078003,   0.03796238490282, &
        0.05097582273167,   0.06327404896598,   0.07468414767371,   0.08504589571635,   0.09421379953849, &
           0.10205916147303,   0.10847184045917,   0.11336181642495,   0.11666044369064,   0.11832141506069, &
              0.11832141550720,   0.11666044317933,   0.11336181703952,   0.10847183974352,   0.10205916225062, &
                0.09421379875775,   0.08504589644177,   0.07468414704828,   0.06327404946527,   0.05097582236737, &
                   0.03796238513729,   0.02441785065909,   0.01055137251011 /)

  Gauss_Legendre(27,1:27,1) = (/  -0.99617926242942,  -0.97992347713087,  -0.95090055655163, &
       -0.90948232147704,  -0.85620790770937,  -0.79177163915352,  -0.71701347370025,  -0.63290797197684, &
         -0.54055156456450,  -0.44114825175375,  -0.33599390363827,  -0.22645936543948,  -0.11397258560953, &
                         0.0,   0.11397258560953,   0.22645936543948,   0.33599390363827,   0.44114825175375, &
                            0.54055156456450,   0.63290797197684,   0.71701347370025,   0.79177163915352, &
                               0.85620790770937,   0.90948232147704,   0.95090055655163,   0.97992347713087, &
                                  0.99617926242942 /)
  Gauss_Legendre(27,1:27,2) = (/   0.00979899557777,   0.02268623238055,   0.03529705318540, &
        0.04744941371414,   0.05898353370023,   0.06974882975360,   0.07960485867054,   0.08842317086199, &
           0.09608871173659,   0.10250165679123,   0.10757826364244,   0.11125251323504,   0.11347631925792, &
              0.11422089514680,   0.11347631865066,   0.11125251431291,   0.10757826230045,   0.10250165820385, &
                 0.09608871038702,   0.08842317207538,   0.07960485762691,   0.06974883061413,   0.05898353302611, &
                    0.04744941420515,   0.03529705286711,   0.02268623254625,   0.00979899552983 /)

  Gauss_Legendre(28,1:28,1) = (/  -0.99644249644426,  -0.98130316813229,  -0.95425927784387, &
       -0.91563302798554,  -0.86589252207191,  -0.80564137096605,  -0.73561087802848,  -0.65665109404250, &
         -0.56972047180434,  -0.47587422495721,  -0.37625151608931,  -0.27206162763497,  -0.16456928213340, &
           -0.05507928988403,   0.05507928988403,   0.16456928213340,   0.27206162763497,   0.37625151608931, &
              0.47587422495721,   0.56972047180434,   0.65665109404250,   0.73561087802848,   0.80564137096605, &
                 0.86589252207191,   0.91563302798554,   0.95425927784387,   0.98130316813229,   0.99644249644426 /)
  Gauss_Legendre(28,1:28,2) = (/   0.00912428160093,   0.02113211400117,   0.03290142715172, &
        0.04427293656342,   0.05510734008512,   0.06527293442030,   0.07464619932439,   0.08311343557742, &
           0.09057172385653,   0.09693067919702,   0.10211294754192,   0.10605578283427,   0.10871118026232, &
              0.11004701881097,   0.11004701396062,   0.10871118487183,   0.10605577870369,   0.10211295099357, &
                 0.09693067653196,   0.09057172574755,   0.08311343434907,   0.07464620004935,   0.06527293403902, &
                    0.05510734025539,   0.04427293650701,   0.03290142715821,   0.02113211400769,   0.00912428159754 /)

  Gauss_Legendre(29,1:29,1) = (/   -0.99667943813982,  -0.98254551553387,  -0.95728558486224, &
       -0.92118024006054,  -0.87463780177302,  -0.81818548865140,  -0.75246285141442,  -0.67821453773872, &
         -0.59628179706397,  -0.50759295516018,  -0.41315288816184,  -0.31403163786990,  -0.21135228616587, &
           -0.10627823013268,                0.0,   0.10627823013268,   0.21135228616587,   0.31403163786990, &
              0.41315288816184,   0.50759295516018,   0.59628179706397,   0.67821453773872,   0.75246285141442, &
                 0.81818548865140,   0.87463780177302,   0.92118024006054,   0.95728558486224,   0.98254551553387, &
                    0.99667943813982 /)
  Gauss_Legendre(29,1:29,2) = (/    0.00851690621500,   0.01973206923479,   0.03074053169285, &
        0.04140200147277,   0.05159490345270,   0.06120300113099,   0.07011803704343,   0.07823820666702, &
           0.08547239630468,   0.09173759919369,   0.09696401001951,   0.10109108272695,   0.10407351140562, &
              0.10587594992052,   0.10647958346671,   0.10587596374542,   0.10407348549588,   0.10109111766346, &
                 0.09696396974561,   0.09173764118770,   0.08547235562978,   0.07823824378824,   0.07011800490835, &
                    0.06120302752581,   0.05159488302551,   0.04140201610381,   0.03074052237167,   0.01973207401840, &
                       0.00851690484311 /)

  Gauss_Legendre(30,1:30,1) = (/   -0.99689348034812,  -0.98366813322554,  -0.96002185270296, &
       -0.92620005787715,  -0.88256052873842,  -0.82956576643007,  -0.76777743002873,  -0.69785049578452, &
         -0.62052618256726,  -0.53662414827065,  -0.44703376952745,  -0.35270472552093,  -0.25463692617151, &
           -0.15386991360830,  -0.05147184255532,   0.05147184255532,   0.15386991360830,   0.25463692617151, &
              0.35270472552093,   0.44703376952745,   0.53662414827065,   0.62052618256726,   0.69785049578452, &
                 0.76777743002873,   0.82956576643007,   0.88256052873842,   0.92620005787715,   0.96002185270296, &
                    0.98366813322554,   0.99689348034812 /)
  Gauss_Legendre(30,1:30,2) = (/    0.00796819580904,   0.01846644941363,   0.02878475462996, &
        0.03879911613800,   0.04840277472189,   0.05749303390173,   0.06597436852008,   0.07375582369324, &
           0.08075605351686,   0.08689962906354,   0.09212267051664,   0.09636861006131,   0.09959351499881, &
              0.10176233797274,   0.10285265553524,   0.10285270116633,   0.10176229387377,   0.09959355625463, &
                 0.09636857258953,   0.09212270367168,   0.08689960040898,   0.08075607773264,   0.07375580370501, &
                    0.06597438457613,   0.05749302143563,   0.04840278396789,   0.03879910972288,   0.02878475862396, &
                       0.01846644739513,   0.00796819638311 /)

  ! Specify the FMM parameters:
  ! * gr_dim (group x-y dimensions) and thus D_max (max. lenght within a group)
  ! * min. far interaction distance
  ! * L_max
  ! * Max. numbers of groups in various dimensions and in total
  gr_dim = 0.32*2.0*pi/k0_fmm
  D_max = SQRT(2.0)*gr_dim
  gr_numx = CEILING((ap_dims(3)-ap_dims(1))/gr_dim)
  gr_numy = CEILING((ap_dims(4)-ap_dims(2))/gr_dim)
  max_groups = gr_numx*gr_numy
!  FarInteract_min = (k0_fmm*D_max + 5.0*LOG(k0_fmm*D_max+PI))/k0_fmm
  FarInteract_min = 2.5*D_max

  print *, 'x_a =', FarInteract_min/D_max

!  L_max = CEILING((k0_fmm*D_max) + 5.0*LOG(k0_fmm*D_max+PI)) !eq.24 - Rokhlin:Pedestrian...
  ! Calculate L_max: (MMB: 'Electromagnetics' ...)
  L_errlevel = 1.0E-4
  L_m = 7.826
  L_c = -1.92*LOG10(L_errlevel) - 0.138
  L_alpha = (1.85*ATAN(2.8*((FarInteract_min/D_max)-1.74)) + 2.455) * 0.8**(-LOG10(L_errlevel))
  L_k = (37.2*EXP(-2.42*SQRT((FarInteract_min/D_max)-1)) - 1.4) * 1.3**(-LOG10(L_errlevel))
  L_max = CEILING(L_k*EXP(-L_alpha*D_max) + L_m*D_max + L_c) - 1

  print *, 'L_max =', L_max
    
  ALLOCATE(edge_grnum(ap_dof))
  ALLOCATE(gr_matindex(max_groups))
  gr_matindex = 0 ! matrix operation

  ! Fill the data structures gr_matindex and edge_grnum:
  ! gr_matindex => is first set to 0/1 if a global group does contain one 
  !                or more dof's. Then the 1's are changed to indices into
  !                the [T]-matrix (as in VTVt).
  ! edge_grnum  => first filled with the global group num's that the
  !                ap_dof's belong to, then overwritten with the group 
  !                matrix co-ordinates.
  ! These data structures are set up by first cycling through all aperture edges
  ! and faces, in order to cover all possible dof's.
  Face_loop: DO ielem = 1,num_apelements
    gl_facenum = elements(ap_elnumbers(ielem))%faces(which_local_face(ielem))
    IF (faces(gl_facenum)%free) THEN ! face contains ap_dof's
      gl_groupnum = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx,facenum=gl_facenum)
      ! Record the (global) group number of the associated face-based dofs and flag the group
      ! as occupied:
      SELECT CASE(faces(gl_facenum)%order)
      CASE(1)
        CONTINUE
      CASE(2)  
        gr_matindex(gl_groupnum) = 1 ! flag this group as occupied
        edge_grnum(renumbered_f1(gl_facenum)) = gl_groupnum
        edge_grnum(renumbered_f2(gl_facenum)) = gl_groupnum
      CASE DEFAULT
        STOP 'IE: Hierarchal order not supported in CBAA_FMM_ALLOCATE.'
      END SELECT
    END IF

    Edge_loop: DO iedge = 1,3
      gl_edgenum = elements(ap_elnumbers(ielem))%edges(which_local_edges(ielem,iedge))
      IF (edges(gl_edgenum)%free) THEN ! edge is an ap_dof
        gl_groupnum = GLOBAL_GROUP_NUMBER(ap_dims(1),ap_dims(2),gr_dim,gr_numx,edgenum=gl_edgenum)
        ! Flag this group as occupied:
        gr_matindex(gl_groupnum) = 1
        ! Record the (global) group number of the associated dof(s):
        SELECT CASE(edges(gl_edgenum)%order)
        CASE(1)
          edge_grnum(renumbered_e1(gl_edgenum)) = gl_groupnum
        CASE(2)
          edge_grnum(renumbered_e1(gl_edgenum)) = gl_groupnum
          edge_grnum(renumbered_e2(gl_edgenum)) = gl_groupnum
        CASE DEFAULT
          STOP 'IE: Hierarchal order not supported in CBAA_FMM_ALLOCATE.'
        END SELECT
      END IF

    END DO Edge_loop
  END DO Face_loop

  ! Fill gr_matindex with the index values:
  grindex = 0
  DO igroup = 1,max_groups
    IF (gr_matindex(igroup).EQ.1) THEN
      grindex = grindex + 1
      gr_matindex(igroup) = grindex
    END IF
  END DO
  gr_num = grindex ! dimension of [T]

  print *,'occupied groups =',gr_num

  ! Convert the global group values in edge_grnum to gr_matindex values:
  DO iedge = 1,ap_dof
    edge_grnum(iedge) = gr_matindex(edge_grnum(iedge))
  END DO

  ! Allocate and fill the sparse storage admin for the [Vt] matrix:
  ALLOCATE(gredge_colind(ap_dof))
  ALLOCATE(gredge_rowind(gr_num+1))
  gredge_rowind(1) = 1
  grindex = 0
  DO igroup = 1,gr_num
    DO iedge = 1,ap_dof
      IF (edge_grnum(iedge).EQ.igroup) THEN
        grindex = grindex + 1
        gredge_colind(grindex) = iedge
      END IF
    END DO
    gredge_rowind(igroup+1) = grindex + 1
  END DO

  ! Allocation of [Z'] storage and admin stuff:
  alloc_val = ap_dof*INT(REAL(ap_dof)/REAL(gr_num))
  ALLOCATE(CBAA_BE_colind(alloc_val))    ! initial estimate
  ALLOCATE(temp_colind(alloc_val))       ! used for resizing colind 
  ALLOCATE(CBAA_BE_rowind(ap_dof+1))     ! Final allocation
  CBAA_BE_rowind(1) = 1
  ALLOCATE(row_dofs(ap_dof))             ! temp variable
  val_tot = 0
  Z_dof_loop: DO iap_dof = 1,ap_dof ! count over all aperture dofs

    ! Find global group number of this aperture dof:
    GROUP_LOOP1: DO igroup = 1,max_groups
      IF (gr_matindex(igroup).EQ.edge_grnum(iap_dof)) THEN
        this_group = igroup
        EXIT GROUP_LOOP1
      END IF
    END DO GROUP_LOOP1 

    colcount = 0
    DO igroup = 1,max_groups ! count over global group numbers
      IF ((gr_matindex(igroup).GT.0).AND.NEAR_NEIGHBOURS(this_group,igroup,gr_dim,gr_numx)) THEN
        tempint = gredge_rowind(gr_matindex(igroup)+1) - &
                  gredge_rowind(gr_matindex(igroup))
        row_dofs(colcount+1:colcount+tempint) = &
          gredge_colind(gredge_rowind(gr_matindex(igroup)):gredge_rowind(gr_matindex(igroup)+1)-1)
        colcount = colcount + tempint
      END IF
    END DO
    CBAA_BE_rowind(iap_dof+1) = &
      CBAA_BE_rowind(iap_dof) + colcount
    ! Write to the colind vector:
    IF (val_tot+colcount.GT.alloc_val) THEN ! not enough space exists
      temp_colind = CBAA_BE_colind
      print *, 'Increasing FMM allocation'
      alloc_val = alloc_val + ap_dof*INT(REAL(ap_dof)/REAL(gr_num))
      DEALLOCATE(CBAA_BE_colind)
      ALLOCATE(CBAA_BE_colind(alloc_val))
      CBAA_BE_colind(1:val_tot) = temp_colind(1:val_tot)
      DEALLOCATE(temp_colind)
      ALLOCATE(temp_colind(alloc_val))
    END IF
    CBAA_BE_colind(val_tot+1:val_tot+colcount) = row_dofs(1:colcount)
    val_tot = val_tot + colcount
  END DO Z_dof_loop

  print *, 'val_tot =',val_tot
  temp_colind = CBAA_BE_colind
  DEALLOCATE(CBAA_BE_colind)
  ALLOCATE(CBAA_BE_colind(val_tot))
  CBAA_BE_colind = temp_colind(1:val_tot)
  DEALLOCATE(temp_colind)
  DEALLOCATE(row_dofs)

  ! Allocate and calculate the spherical integration vectors:
  k_integration_scheme = 3

  SELECT CASE(k_integration_scheme)

  CASE(1) ! Trapezoidal in phi and theta:
    kdir_ppoints = 2*(L_max+1)
    kdir_tpoints = L_max+1
    dkdir_phi = 2.0*PI/REAL(kdir_ppoints)
    dkdir_theta = 0.5*PI/REAL(kdir_tpoints)
    num_k_dirs = kdir_ppoints*kdir_tpoints ! a global variable
    ALLOCATE(k_dirs(num_k_dirs,4)) ! / x , y , z , weight /

    DO tempint = 1,kdir_ppoints
      kdir_phi = REAL(tempint)*dkdir_theta-dkdir_theta/2.0
       
      DO tempint2 = 1,kdir_tpoints
        kdir_theta = REAL(tempint2)*dkdir_phi-dkdir_phi/2.0

        ! ( x,y,z-Direction of ^k,Weight of ^k-direction ):
        k_dirs((tempint-1)*kdir_tpoints+tempint2,1:4) =                                        &
          (/ SIN(kdir_theta)*COS(kdir_phi) , SIN(kdir_theta)*SIN(kdir_phi) , COS(kdir_theta) , &
             dkdir_phi*(COS(kdir_theta-dkdir_theta/2.0)-COS(kdir_theta+dkdir_theta/2.0)) /)
        IF (ABS(k_dirs((tempint-1)*kdir_tpoints+tempint2,3)).GT.EPS) &
          k_dirs((tempint-1)*kdir_tpoints+tempint2,4) = 2.0*k_dirs((tempint-1)*kdir_tpoints+tempint2,4)
      END DO
    END DO

  CASE(2) ! Trapezoidal in phi and Gauss-Legendre in theta(0<theta<pi/2) (See Delves and Mohamed):
    kdir_ppoints = 2*L_max
    kdir_tpoints = L_max+1 !L_max+1 !CEILING(REAL(L_max+1)/2.0)
    num_k_dirs = kdir_ppoints*kdir_tpoints ! a global variable
    ALLOCATE(k_dirs(num_k_dirs,4)) ! / x , y , z , weight /

    DO tempint = 1,kdir_tpoints
      kdir_theta = 0.25*PI*(1.0 + Gauss_Legendre(kdir_tpoints,tempint,1))
      
      DO tempint2 = 1,kdir_ppoints
        kdir_phi = REAL(tempint2)*2.0*PI/REAL(kdir_ppoints)

        ! ( x,y,z-Direction of ^k,Weight of ^k-direction and weight):
        k_dirs((tempint-1)*kdir_ppoints+tempint2,1:4) =                                        &
          (/ SIN(kdir_theta)*COS(kdir_phi) , SIN(kdir_theta)*SIN(kdir_phi) , COS(kdir_theta) , &
             (PI**2)*SIN(kdir_theta)*Gauss_Legendre(kdir_tpoints,tempint,2)/REAL(kdir_ppoints) /)
      END DO
    END DO

  CASE(3) ! Trapezoidal in phi and Gauss-Legendre in theta(0<theta<pi) (See Delves and Mohamed):
    kdir_ppoints = CEILING(1.6*REAL(L_max) + 1.6)
    kdir_torder = L_max + 6
    kdir_tpoints = CEILING(0.5*REAL(kdir_torder)) !L_max+1 !CEILING(REAL(L_max+1)/2.0)
    num_k_dirs = kdir_ppoints*kdir_tpoints ! a global variable
    ALLOCATE(k_dirs(num_k_dirs,4)) ! / x , y , z , weight /

    print *, 'kdir_tpoints =',kdir_tpoints 
    print *, 'kdir_ppoints =',kdir_ppoints 
    print *, 'num_kdirs =',num_k_dirs

    DO tempint = 1,kdir_tpoints
      kdir_theta = 0.5*PI*(1.0 + Gauss_Legendre(kdir_torder,tempint,1))
     
      DO tempint2 = 1,kdir_ppoints
        kdir_phi = REAL(tempint2)*2.0*PI/REAL(kdir_ppoints)

        ! ( x,y,z-Direction of ^k,Weight of ^k-direction and weight):
        k_dirs((tempint-1)*kdir_ppoints+tempint2,1:4) =                                        &
          (/ SIN(kdir_theta)*COS(kdir_phi) , SIN(kdir_theta)*SIN(kdir_phi) , COS(kdir_theta) , &
             (PI**2)*SIN(kdir_theta)*Gauss_Legendre(kdir_torder,tempint,2)/REAL(kdir_ppoints) /)

        IF (ABS(k_dirs((tempint-1)*kdir_ppoints+tempint2,3)).GT.EPS) &
          k_dirs((tempint-1)*kdir_ppoints+tempint2,4) = 2.0*k_dirs((tempint-1)*kdir_ppoints+tempint2,4)

      END DO
    END DO

  END SELECT

  DEALLOCATE(Gauss_Legendre) 

  ! Checking the integration rule just created:
  print *, '4*PI =', SUM(k_dirs(1:num_k_dirs,4))  

  ! Allocate storage for [Z'], [V], [T] values:
  ALLOCATE(CBAA_BE_val(val_tot))
  ALLOCATE(CBAA_groupmat(num_k_dirs,gr_num,gr_num))
  ALLOCATE(CBAA_elemvec_term1(num_k_dirs,ap_dof))    ! scalar values
  ALLOCATE(CBAA_elemvec_term2(num_k_dirs,ap_dof,3))  ! vector values

  print *,'fmm_mem in bytes =', 8*(val_tot + 4*num_k_dirs*ap_dof + num_k_dirs*gr_num*gr_num)
  print *,'bi_mem in bytes  =', 8*ap_dof*ap_dof

END SUBROUTINE CBAA_FMM_ALLOCATE
!*******************************************************************************


SUBROUTINE MATVECPROD_FMM(input,output)
  USE CBAA_data
  IMPLICIT NONE
!*******************************************************************************
! Calculates the matrix-vector product of the FMM BI matrix factors with a given
! vector of the same dimension.
! MMB 11 Apr 2001. Copied/modified code from the old CBAA_BCG_SOLVE.
!*******************************************************************************
  COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: input    ! x in Ax=y
  COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: output  ! y in Ax=y

  INTEGER(I4B) :: vec_dim,mulcount,kcount     ! dimension of the vectors (x & y); counters
  COMPLEX(SPC),ALLOCATABLE,DIMENSION(:) :: temp_fmm_scalar
  COMPLEX(SPC),ALLOCATABLE,DIMENSION(:,:) :: temp_fmm_vector
  COMPLEX(SPC),ALLOCATABLE,DIMENSION(:) :: input_work
  LOGICAL(LGT) :: only_near_interact

  vec_dim = SIZE(input) ! establish the vector length (equal to <ap_dof>, but this is more general)
  ALLOCATE(input_work(vec_dim))  
  ALLOCATE(temp_fmm_scalar(gr_num))
  ALLOCATE(temp_fmm_vector(gr_num,2)) ! 2=>only x&y vector components in the CBAA aperture.
  only_near_interact = .FALSE.

  ! Add the near-interaction contribution from [Z']:
  ! (<output> is not initialized, because all values are assigned here.)
  DO mulcount = 1,vec_dim
    output(mulcount) = SUM(CBAA_BE_val(CBAA_BE_rowind(mulcount):CBAA_BE_rowind(mulcount+1)-1)* &
      input(CBAA_BE_colind(CBAA_BE_rowind(mulcount):CBAA_BE_rowind(mulcount+1)-1)))
  END DO

  ! calculate [V][T][Vt] contribution:
  ! Conjugation beforehand, to improve speed, the input will be conjugated again, after the product,
  ! to restore it to its original, correct value:
  input_work = CONJG(input)

  DO kcount = 1,num_k_dirs ! spherical quadrature
    IF (only_near_interact) CYCLE

    ! calculate temp_fmm_... = [Vt][pk]
    ! For speed improvement: CONJG the input vector and CONJG the product afterwards
    ! in order to decrease the number of conjugations:
    DO mulcount = 1,gr_num
      temp_fmm_scalar(mulcount) =                                                     &
        SUM(CBAA_elemvec_term1(kcount,gredge_colind(gredge_rowind(mulcount):          &
        gredge_rowind(mulcount+1)-1)) *                                               &
        input_work(gredge_colind(gredge_rowind(mulcount):                             &
        gredge_rowind(mulcount+1)-1)))
  
      temp_fmm_vector(mulcount,1) =                                                   &
        SUM(CBAA_elemvec_term2(kcount,gredge_colind(gredge_rowind(mulcount):          &
        gredge_rowind(mulcount+1)-1),1) *                                             &
        input_work(gredge_colind(gredge_rowind(mulcount):                             &
        gredge_rowind(mulcount+1)-1)))

      temp_fmm_vector(mulcount,2) =                                                   &
        SUM(CBAA_elemvec_term2(kcount,gredge_colind(gredge_rowind(mulcount):          &
        gredge_rowind(mulcount+1)-1),2) *                                             &
        input_work(gredge_colind(gredge_rowind(mulcount):                             &
        gredge_rowind(mulcount+1)-1)))
    END DO
    ! now the product is conjugated:
    temp_fmm_scalar = CONJG(temp_fmm_scalar)
    temp_fmm_vector = - CONJG(temp_fmm_vector) ! minus sign to counteract the erronous
                                               ! minus sign that is intruduced by 
                                               ! conjugation of the scaling factor in Jin,
                                               ! eq. 9.136 (sf = sqrt(-2*k_0^2)).
               
    ! calculate temp_fmm_... = [T][Vt][pk] 
    temp_fmm_scalar      = MATMUL(CBAA_groupmat(kcount,:,:),temp_fmm_scalar)
    temp_fmm_vector(:,1) = MATMUL(CBAA_groupmat(kcount,:,:),temp_fmm_vector(:,1))
    temp_fmm_vector(:,2) = MATMUL(CBAA_groupmat(kcount,:,:),temp_fmm_vector(:,2))
  
    ! Calculate [V][T][Vt][pk] for this k-direction and add to the total:
    DO mulcount = 1,vec_dim
      output(mulcount) = output(mulcount) + k_dirs(kcount,4)*                                   &
        ( CBAA_elemvec_term1(kcount,mulcount)*temp_fmm_scalar(edge_grnum(mulcount)) +           &
        SUM(CBAA_elemvec_term2(kcount,mulcount,1:2)*temp_fmm_vector(edge_grnum(mulcount),1:2)) )
    END DO
  END DO

  DEALLOCATE(temp_fmm_scalar)
  DEALLOCATE(temp_fmm_vector)

  DEALLOCATE(input_work)  

print *,'exiting into MATVECPRD_FMM...'


END SUBROUTINE MATVECPROD_FMM
!*******************************************************************************

SUBROUTINE CBAA_MAKE_COAX_A_AND_B(flag)
  USE coax_feed
  USE feminterface, ONLY: CONVERTCOR
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Calculates the b-vector and A-matrix contributions for the coax feed 
! technique presented in IEEE Trans. A&P, vol.43, p.1474 - Gong & Volakis.
! <flag> indicates which contribution must be added (A_mat. or b-vec.).
! Created: May 2000 MMB
! 2001 Jan 21: Made into a interface routine. MMB
!*******************************************************************************
  INTEGER(I4B),INTENT(IN) :: flag            ! 0 => add the A-matrix contribution
                                             ! 1 => add the b_vector contribution

  COMPLEX(SPC) :: Amat_tempterm,bvec_tempterm
  INTEGER(I4B) :: iedge,jedge,column,row,iport   ! counters
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)      ! The complex constant j
  REAL(SP), DIMENSION(3) :: unit_vec             ! unit vec. to describe inner conductor direction

  DO AX_counter = 1,NUM_AX_cards ! cycle through all AX cards
  
    DO iedge = 1,AXdata(AX_counter)%coaxedges_num
      row = renumbered_e1(AXdata(AX_counter)%coaxedges(iedge))

      ! Add b-vector contribution:
      bvec_tempterm = cj*2.0*k0*Z_zero*(AXdata(AX_counter)%coax_b-AXdata(AX_counter)%coax_a)*    &
                        AXdata(AX_counter)%coax_Iabs* &
                        EXP(cj*AXdata(AX_counter)%coax_Iphase*PI/180.0)/AXdata(AX_counter)%coaxedges_num
      IF (SQRT(SUM((vertices(edges(AXdata(AX_counter)%coaxedges(iedge))%nodes(2))%coord- &
        AXdata(AX_counter)%coaxcentre)**2)).LT.EPS) THEN 
        ! Edge vector (iedge) points toward coax centre
        bvec_tempterm = - bvec_tempterm
      END IF 

      IF (flag.EQ.1) THEN ! this contribution was requested
        b_vec_c(row) = b_vec_c(row) + bvec_tempterm
      END IF

      DO jedge = 1,AXdata(AX_counter)%coaxedges_num
        column = renumbered_e1(AXdata(AX_counter)%coaxedges(jedge))

        ! Add A-matrix contribution:
        Amat_tempterm = cj*2.0*PI*k0*SQRT(AXdata(AX_counter)%coax_eps)* &
          ((AXdata(AX_counter)%coax_b-AXdata(AX_counter)%coax_a)**2)/((AXdata(AX_counter)%coaxedges_num)* &
          LOG(AXdata(AX_counter)%coax_b/AXdata(AX_counter)%coax_a))
        IF (SQRT(SUM((vertices(edges(AXdata(AX_counter)%coaxedges(iedge))%nodes(2))%coord- &
          AXdata(AX_counter)%coaxcentre)**2)).LT.EPS) THEN 
          ! Edge vector (iedge) points toward coax centre
          Amat_tempterm = - Amat_tempterm
        END IF  
        IF (SQRT(SUM((vertices(edges(AXdata(AX_counter)%coaxedges(jedge))%nodes(2))%coord- &
          AXdata(AX_counter)%coaxcentre)**2)).LT.EPS) THEN 
          ! Edge vector (jedge) points toward coax centre
          Amat_tempterm = - Amat_tempterm
        END IF  

        IF (row.EQ.column) THEN ! add to the diagonal
          IF (flag.EQ.0) THEN ! this contribution was requested
            IF (SPARSE) THEN
              Asparse_c(CONVERTCOR(row,column)) = Asparse_c(CONVERTCOR(row,column)) + Amat_tempterm
            ELSE
              A_mat_c(row,column) = A_mat_c(row,column) + Amat_tempterm



print *, Amat_tempterm


           END IF
          END IF
        END IF

      END DO
    END DO

    ! Report this action in the output file:
    IF (flag.EQ.1) THEN ! only report in the b-vectro contribution case
      WRITE (FILEOUT,'(//,20X,A,/)') 'COAXIAL EXCITATION AT THE FEM MESH BOUNDARY'
      WRITE (FILEOUT,'(1X,A,I4)')    'Number of source:                  N = ', 0
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Frequency in Hz:                FREQ = ', frequency
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Wavelength in m:              LAMBDA = ', lambda0
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Matched load current in A:      |I0| = ', AXdata(AX_counter)%coax_Iabs
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Phase in deg.:               ARG(I0) = ', AXdata(AX_counter)%coax_Iphase
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Inner conductor rad. in m:         a = ', AXdata(AX_counter)%coax_a
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Outer conductor rad. in m:         b = ', AXdata(AX_counter)%coax_b
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Relative eps_r of coax:          EPS = ', AXdata(AX_counter)%coax_eps
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Location of aperture center in m:  X = ', AXdata(AX_counter)%coaxcentre(1)
      WRITE (FILEOUT,'(1X,A,E12.5)') '                                   Y = ', AXdata(AX_counter)%coaxcentre(2)
      WRITE (FILEOUT,'(1X,A,E12.5)') '                                   Z = ', AXdata(AX_counter)%coaxcentre(3)
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Inner conductor len. in m:       LEN = ', AXdata(AX_counter)%coaxlen
      unit_vec = 0.0
      unit_vec(AXdata(AX_counter)%coaxdir) = 1.0
      WRITE (FILEOUT,'(1X,A,E12.5)') 'Inner conductor +direction:       UX = ', unit_vec(1)
      WRITE (FILEOUT,'(1X,A,E12.5)') '                                  UY = ', unit_vec(2)
      WRITE (FILEOUT,'(1X,A,E12.5)') '                                  UZ = ', unit_vec(3)
    END IF

  END DO

END SUBROUTINE CBAA_MAKE_COAX_A_AND_B
!*******************************************************************************


SUBROUTINE CBAA_MAKE_COAX_AMATRIX(k0)
  USE B_matrix, ONLY: B_MAKE_HIERARCHAL
  USE coax_feed
  USE feminterface, ONLY: CONVERTCOR,LOCAL_TO_GLOBAL_INDEX_TRI
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! Assembles the elemental, coax port contributions into the A-matrix. The general,
! dominant mode coax formulation is used. MMB Ph.D. 2002
! 2002-02-07: Created. MMB
! Changed DBD 07 March 2002
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k0  ! wavenumber at which calculations must be performed

! Next line changed DBD 07 March 2002
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,ELEM_TRI_MATRIX_SIZE) :: Bs 
  LOGICAL(LGT) :: port_face_found
  INTEGER(I4B) :: ielem,jface,port_num
  INTEGER(I4B) :: mm,nn,row,column,indpos
  INTEGER(I4B), DIMENSION(3) :: localfaceedges
  REAL(SP), DIMENSION(3) :: normal    ! Outward directed unit normal on port.
  REAL(SP) :: kc                      ! coax porpogation constant
  COMPLEX(SPC) :: gamma               ! constant multiplier in elemental contribution

  ! Compute the elemental [B] matrix if this element has a face on a port. 
  ! Note that an element should NOT have more than one face on a port,
  ! nor have faces on more than one port.

  ! Cycle through all element faces and add the contributions of those 
  ! that lie on a port:
  ELEMENT_LOOP: DO ielem = 1,num_elements

    ! Initialize for this element:
    port_face_found = .FALSE.

    FACE_LOOP: DO jface = 1,4
      
	  ! Cycle if not part of a coax aperture:
	  IF (.NOT.faces(elements(ielem)%faces(jface))%coax_aperture) CYCLE FACE_LOOP

      ! Error check and record the port face's data:
      IF (port_face_found) CALL ERROR_FEMFEKO(1,4504,int1=ielem)
      port_face_found = .TRUE.
      port_num = faces(elements(ielem)%faces(jface))%coaxnumber
      kc = k0*SQRT(AXdata(port_num)%coax_mu)*SQRT(AXdata(port_num)%coax_eps)
      gamma = j*kc/AXdata(port_num)%coax_mu

      ! Calculate the elemental port face matrix:
      CALL B_MAKE_HIERARCHAL(ielem,jface,AXdata(port_num)%normal,Bs)

! TEMP!!!!
!print *, Bs(1:3,1:3)
!pause
!print *,'gamma =',gamma
!print *,'port_num =',port_num


      ! Now add <Bs> to the system matrix:
      DO mm = 1,8 ! row counter

        ! Assign the global row:
        row = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,jface,mm)

        DO nn = 1,8 ! column counter

          ! Assign the global column:
          column = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,jface,nn)

          ! Add this elemental matrix element to the system matrix:
          IF ((row.GT.0).AND.(column.GT.0)) THEN ! this is two dofs:
            IF (SPARSE) THEN
              indpos = CONVERTCOR(row,column)
              Asparse_c(indpos) = Asparse_c(indpos) + gamma*Bs(mm,nn)
            ELSE
              A_mat_c(row,column)  = A_mat_c(row,column) + gamma*Bs(mm,nn)


! TEMP!!!!
!print *,Bs(mm,nn)


            END IF
          END IF

        END DO
      END DO
    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE CBAA_MAKE_COAX_AMATRIX
!*******************************************************************************


SUBROUTINE CBAA_MAKE_COAX_BVECTOR(k0)
  USE coax_feed
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI
  USE geometry
  USE matrix
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Assembles the elemental, coax port contributions into the b-vector. The general,
! dominant mode coax formulation is used. MMB Ph.D. 2002.
!
! 2002-02-07: Created. MMB
! 2002-05-09: Generalized to arbitrary order elements. MMB.
!*******************************************************************************
  REAL(SP), INTENT(IN) :: k0  ! wavenumber at which calculations must be performed

  INTEGER(I4B) :: ielem,iface,port_num,iquad,mm,dof_num ! counters
  REAL(SP),DIMENSION(ELEM_TRI_MATRIX_SIZE) :: mode_vbf_integrals
  COMPLEX(SPC),DIMENSION(ELEM_TRI_MATRIX_SIZE) :: face_bvec_contrib
  INTEGER(I4B),DIMENSION(3) :: facenodes
  COMPLEX(SPC) :: b_val
  REAL(SP) :: kc ! coax porpogation constant

  ! Cycle trough all faces, element-wise:
  ELEMENT_LOOP: DO ielem = 1,num_elements
	FACE_LOOP: DO iface = 1,4
      
	  ! Cycle if this face is not part of a coax port:
	  IF (.NOT.faces(elements(ielem)%faces(iface))%coax_aperture) CYCLE FACE_LOOP

      ! Record the port number of the current port face:
      port_num = faces(elements(ielem)%faces(iface))%coaxnumber

      ! Integrate over the face the dot product of ^n X the face basis functions
	  ! with ^n X the coax mode:
	  CALL CBAA_MAKE_COAX_MODE_INTEGRAL(ielem,iface,AXdata(port_num)%normal, &
	                                    AXdata(port_num)%coaxcentre,         &
										AXdata(port_num)%coax_a,             &
										mode_vbf_integrals)

      ! Weight the coax mode so that it corresponds to the incident electric 
	  ! field:
      face_bvec_contrib = AXdata(port_num)%coax_Iabs *                               &
                          EXP(j*AXdata(port_num)%coax_Iphase*PI/180.0) *             &
					      (Z_zero/(2.0*PI*AXdata(port_num)%coax_a)) *                &
					      SQRT(AXdata(port_num)%coax_mu/AXdata(port_num)%coax_eps) * &
					      mode_vbf_integrals

      ! Scale by factors outside the integral:
      kc = k0*SQRT(AXdata(port_num)%coax_mu)*SQRT(AXdata(port_num)%coax_eps)
      face_bvec_contrib = 2*j*(kc/AXdata(port_num)%coax_mu) * face_bvec_contrib

      ! Add the contribution of the current face to the excitation (b) vector:
      DO mm = 1,ELEM_TRI_MATRIX_SIZE
		dof_num = LOCAL_TO_GLOBAL_INDEX_TRI(ielem,iface,mm)
        IF (dof_num.GT.0) THEN
		  b_vec_c(dof_num) = b_vec_c(dof_num) + face_bvec_contrib(mm)
        END IF
	  END DO

    END DO FACE_LOOP
  END DO ELEMENT_LOOP

END SUBROUTINE CBAA_MAKE_COAX_BVECTOR
!*******************************************************************************


SUBROUTINE CBAA_MAKE_COAX_MODE_INTEGRAL(elem,local_face,normal,coaxcenter,radius_a, &
                                        int_values)
  USE basis_function, ONLY: EVALUATE_ELEMENTAL_FUNCTIONS
  USE geometry
  USE math_tools, ONLY: CROSS_PRODUCT,VECTOR_LENGTH
  USE nrtype
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! Integrate the dot product of ^n X (coax mode function) and ^n X (face VBF's),
! over local face <local_face>, of element <elem>. Returns the integral values
! in <int_values>. Max. dim. of <int_values> is 8, which is for the LT/QN case.
!
! 2002-02-09: Created. MMB
! 2002-03-07: Started to provide variable dimensions. DBD.
! 2002-05-09: Generalized to arbitrary order elements by adding variable 
!             dimensions for local number of VBF's. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face
  REAL(SP), DIMENSION(3), INTENT(IN) :: normal,coaxcenter
  REAL(SP), INTENT(IN) :: radius_a
  REAL(SP), DIMENSION(ELEM_TRI_MATRIX_SIZE), INTENT(OUT) :: int_values
  
  INTEGER(I4B) :: num_qpts               ! number of quadrature points
  INTEGER(I4B) :: inode,iquad,ifunc,mm   ! counters
  REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: vbf_values
  REAL(SP),DIMENSION(ELEM_TET_MATRIX_SIZE) :: vbf_integrals,vbf_integrals_temp
  REAL(SP),DIMENSION(4) :: r_vec,s_vec
  REAL(SP),DIMENSION(4,4) :: vertmat,vertmattemp
  INTEGER(I4B),DIMENSION(3) :: facenodes
  REAL(SP),DIMENSION(3) :: coax_mode_term,node_coord
  LOGICAL(LGT) :: node_on_center
  INTEGER(I4B) :: facecenternode
  INTEGER(I4B),DIMENSION(2) :: faceouternodes
  REAL(SP),DIMENSION(2,3) :: newnodes
  REAL(SP) :: templength
  REAL(SP),DIMENSION(2) :: temp_areas
  
  num_qpts      = 7  ! set the number of quadrature points
  int_values    = 0.0 ! vector operation, initialise

  ! Establish vertices matrix for repeated use:
  CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vertmat)
  facenodes = LOCAL_FACENODES(local_face)
  s_vec = 0.0 ! vector assignment

  ! Establish whether one on the three nodes coincides with the coax
  ! center:
  node_on_center = .FALSE. ! init
  DO inode = 1,3 
    node_coord = vertices(elements(elem)%nodes(facenodes(inode)))%coord
    IF (VECTOR_LENGTH(node_coord-coaxcenter).LT.EPS) THEN
	  node_on_center = .TRUE.
      facecenternode = inode
	END IF
  END DO  
  
  vbf_integrals = 0.0 ! vector operation, initialise the integrals
  
  IF (.NOT.node_on_center) THEN ! normal quadrature, no singularity

    QUAD_LOOP: DO iquad = 1,num_qpts
        
      ! Calculate the xyz coords of current quad point:
      s_vec(facenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
      r_vec = MATMUL(vertmat,s_vec)

      ! Evaluate basis functions:
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values)

      ! Cross with the normal:        
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	    vbf_values(ifunc,1:3) = CROSS_PRODUCT(normal,vbf_values(ifunc,1:3))
      END DO

      ! Evaluate ^n X (incident coax mode):
      coax_mode_term = COAX_MODE(coaxcenter,radius_a,r_vec(1:3))
      coax_mode_term = CROSS_PRODUCT(normal,coax_mode_term)

      ! Add to the total 
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
        vbf_integrals(ifunc) = vbf_integrals(ifunc) +                      &
		                       quad_tri_rules(num_qpts)%rule(iquad,4) *    &
		                       SUM(coax_mode_term * vbf_values(ifunc,1:3))
      END DO
    END DO QUAD_LOOP
    vbf_integrals = FACE_AREA(elem,local_face) * vbf_integrals ! final step in quadrature: scale by area

  ELSE ! <node_on_center> is true, the singularity will now be delt with,
       ! by splitting adding two node at the a-radius distance fron the
	   ! coax center, on the two radial edges. These two nodes plus the two
	   ! outer, triangular, facial node, form a quadrilateral thar will be 
	   ! split into two triangles, on which normal quadrature can be used.
    
	! Create faceouternodes datastructure:
    IF (node_on_center) THEN
      SELECT CASE (facecenternode)
      CASE (1)
        faceouternodes = (/ 2, 3 /)
	  CASE (2)
        faceouternodes = (/ 1, 3 /)
      CASE (3)
        faceouternodes = (/ 1, 2 /)
      END SELECT
    END IF
    facecenternode = facenodes(facecenternode)
    faceouternodes = facenodes(faceouternodes)

    ! Calculate the two new nodes' coordinates and the new face areas:
	DO inode = 1,2
      node_coord = vertices(elements(elem)%nodes(faceouternodes(inode)))%coord
	  templength = VECTOR_LENGTH(coaxcenter-node_coord)
	  newnodes(inode,1:3) =   (radius_a/templength)*node_coord               &
	                        + ((templength-radius_a)/templength)*coaxcenter
      SELECT CASE (inode)
	  CASE (1)
	    temp_areas(1) = ((templength-radius_a)/templength) * &
		                FACE_AREA(elem,local_face)
	  CASE (2)
	    temp_areas(2) = ((templength-radius_a)/templength) * &
		                (FACE_AREA(elem,local_face) - temp_areas(1))
      END SELECT	  
	END DO

    ! Replace center node with new node 1:
    vertmat(1:3,facecenternode) = newnodes(1,1:3)
    
    vbf_integrals_temp = 0.0 ! vector operation, initialise the integrals

    QUAD_LOOP_1: DO iquad = 1,num_qpts
      ! Calculate the xyz coords of current quad point:
      s_vec(facenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
      r_vec = MATMUL(vertmat,s_vec)

      ! Evaluate basis functions:
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values)

      ! Cross with the normal:        
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	    vbf_values(ifunc,1:3) = CROSS_PRODUCT(normal,vbf_values(ifunc,1:3))
      END DO

      ! Evaluate ^n X (incident coax mode):
      coax_mode_term = COAX_MODE(coaxcenter,radius_a,r_vec(1:3))
      coax_mode_term = CROSS_PRODUCT(normal,coax_mode_term)

      ! Add to the total 
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
        vbf_integrals_temp(ifunc) = vbf_integrals_temp(ifunc) +                      &
		                            quad_tri_rules(num_qpts)%rule(iquad,4) *    &
		                            SUM(coax_mode_term * vbf_values(ifunc,1:3))
      END DO
    END DO QUAD_LOOP_1
    ! final step in quadrature: scale by area. Then add to total for both faces:
	vbf_integrals = vbf_integrals + temp_areas(1) * vbf_integrals_temp

    ! Replace outer face node 1 with new node 2:
    vertmat(1:3,faceouternodes(1)) = newnodes(2,1:3)
    
    vbf_integrals_temp = 0.0 ! vector operation, initialise the integrals

    QUAD_LOOP_2: DO iquad = 1,num_qpts
      ! Calculate the xyz coords of current quad point:
      s_vec(facenodes) = quad_tri_rules(num_qpts)%rule(iquad,1:3)
      r_vec = MATMUL(vertmat,s_vec)

      ! Evaluate basis functions:
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(elem,r_vec(1:3),0,vbf_values)

      ! Cross with the normal:        
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
	    vbf_values(ifunc,1:3) = CROSS_PRODUCT(normal,vbf_values(ifunc,1:3))
      END DO

      ! Evaluate ^n X (incident coax mode):
      coax_mode_term = COAX_MODE(coaxcenter,radius_a,r_vec(1:3))
      coax_mode_term = CROSS_PRODUCT(normal,coax_mode_term)

      ! Add to the total 
      DO ifunc = 1,ELEM_TET_MATRIX_SIZE
        vbf_integrals_temp(ifunc) = vbf_integrals_temp(ifunc) +                      &
		                            quad_tri_rules(num_qpts)%rule(iquad,4) *    &
		                            SUM(coax_mode_term * vbf_values(ifunc,1:3))
      END DO
    END DO QUAD_LOOP_2
    ! final step in quadrature: scale by area. Then add to total for both faces:
    vbf_integrals = vbf_integrals + temp_areas(2) * vbf_integrals_temp

  END IF

  ! Convert the tetrahedral values to the facial values:
  DO mm = 1,ELEM_TRI_MATRIX_SIZE
    int_values(mm) = vbf_integrals(LOCAL_TRI_TO_LOCAL_TET_INDEX(local_face,mm))
  END DO

END SUBROUTINE CBAA_MAKE_COAX_MODE_INTEGRAL
!*******************************************************************************


FUNCTION COAX_MODE(center,radius_a,r_vec)
  USE math_tools, ONLY: VECTOR_LENGTH
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Returns the vector-valued, coaxial, dominant, normalised, TEM mode function. 
! The coax aperture is described by its center <center> and inner radius 
! <radius_a>. IMPORTANT: this routine assumes that both points does lie in the
! coax port aperture.
! 2002-02-08: Created. MMB
!*******************************************************************************
  REAL(SP), DIMENSION(3), INTENT(IN) :: center,r_vec
  REAL(SP), INTENT(IN) :: radius_a
  REAL(SP), DIMENSION(3) :: COAX_MODE

  COAX_MODE = r_vec - center
  COAX_MODE = (radius_a/(VECTOR_LENGTH(COAX_MODE)**2)) * COAX_MODE

END FUNCTION COAX_MODE
!*******************************************************************************


SUBROUTINE CBAA_MAKE_CURRENTPROBE_BVECTOR(prx0,pry0,prz0,prdir,prlen,prrad,prabs,prphase)
  USE basis_function, ONLY: EVALUATE_ELEMENTAL_FUNCTIONS
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TET
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Adds the contribution of a current probe in the FEM mesh to the b-vector. 
! The probe's information is passed via parameters and is described below.
! MMB 18 Jan 2001
!*******************************************************************************
  REAL(SP), INTENT(IN) :: prabs          ! current magnitude
  REAL(SP), INTENT(IN) :: prphase        ! current phase in radians
  REAL(SP), INTENT(IN) :: prlen          ! length of probe
  REAL(SP), INTENT(IN) :: prx0,pry0,prz0 ! xyz of starting point (current defined as flowing away from this point)
  REAL(SP), INTENT(IN) :: prrad          ! radius of probe
  INTEGER(I4B), INTENT(IN) :: prdir      ! x/y/z : 1/2/3 direction of probe

  REAL(SP) :: prrad_temp                   ! radius of probe used in this routine
  INTEGER(I4B) :: prnum                    ! circumference-wise discretization
  INTEGER(I4B) :: prdisc                   ! length-wise discretization
  INTEGER(I4B) :: iint,jint,probe_counter, &
                  ifunc,row                ! counters
  REAL(SP) :: prx,pry,prz                  ! position integration variables
  INTEGER(I4B) :: inside_el                ! element number in which the near field point is located
  REAL(SP), DIMENSION(3) :: r_vec          ! integration point xyz co-ordinates for EVALUATE_ELEMENTAL_FUNCTIONS
! Next line changed 07 Mar 2002 DBD.
  REAL(SP), DIMENSION(ELEM_TET_MATRIX_SIZE,3) :: func_values ! result of EVALUATE_ELEMENTAL_FUNCTIONS
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)! The complex constant j
  REAL(SP),DIMENSION(3) :: unit_vec        ! represents the feed direction

  ! Calculate the number of integration points, with the following guidelines:
  ! - Use a discretization of lambda0/1000
  ! - At least 20 points in the lenght-wise integration
  ! - All current is assumed to flow on the probe surface
  prdisc = CEILING(1000.0*prlen/lambda0)
  prrad_temp = prrad
  IF (prdisc.LT.20) THEN
    prdisc = 20
  END IF
  prnum = CEILING(2.0*PI*prrad_temp*REAL(prdisc)/prlen)
  IF (prnum.LE.1) THEN ! the probe is infinitely thin and located at the tube center
    prnum = 1
    prrad_temp = 0.0
  END IF

  ! Calculate integral: (eq. 9.147 - jin)
  Probe_lateral_integration: DO jint = 1,prnum ! Integrating around the circumference
    SELECT CASE(prdir) ! increment the integration position
    CASE(1) 
      pry = pry0 + prrad_temp*COS(jint*2*PI/prnum)
      prz = prz0 + prrad_temp*SIN(jint*2*PI/prnum)
    CASE(2) 
      prz = prz0 + prrad_temp*COS(jint*2*PI/prnum)
      prx = prx0 + prrad_temp*SIN(jint*2*PI/prnum)
    CASE(3) 
      prx = prx0 + prrad_temp*COS(jint*2*PI/prnum)
      pry = pry0 + prrad_temp*SIN(jint*2*PI/prnum)
    END SELECT

    Probe_lengthwise_integration: DO iint = 0,prdisc-1 ! Integrating along the length
      SELECT CASE(prdir) ! increment the integration position
      CASE(1) 
        prx = prx0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
      CASE(2) 
        pry = pry0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
      CASE(3) 
        prz = prz0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
      END SELECT

      ! Find the element that the point belongs to:
      inside_el = XYZ_TO_ELNUM(prx,pry,prz)

      ! Check that the point is inside the mesh:
      IF (inside_el.EQ.0) CYCLE Probe_lengthwise_integration ! Need to create a warning here
             
      ! Evaluate all the basis functions at this point:
      r_vec = (/ prx , pry , prz /)
      CALL EVALUATE_ELEMENTAL_FUNCTIONS(inside_el,r_vec,0,func_values)

      ! Add the contribution of this integration point to the b-vector:
      basis_function_loop: DO ifunc = 1,ELEM_TET_MATRIX_SIZE
        row = LOCAL_TO_GLOBAL_INDEX_TET(inside_el,ifunc)
        IF (row.GT.0) THEN ! this is a degree of freedom of the system
          b_vec_c(row) = b_vec_c(row) &
                         - cj*k0*Z_zero*(prlen/REAL(prdisc))*func_values(ifunc,prdir)* &
                           (prabs/REAL(prnum))*EXP(cj*prphase)
        END IF
      END DO basis_function_loop

    END DO Probe_lengthwise_integration
  END DO Probe_lateral_integration

  ! Write the exitation information to the output file:
  WRITE (FILEOUT,'(//,20X,A,/)') 'EXCITATION BY CURRENT PROBE SOURCE IN THE FEM MESH'
  WRITE (FILEOUT,'(1X,A,I4)')    'Number of current source:         N = ', 0
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Frequency in Hz:               FREQ = ', frequency
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Wavelength in m:             LAMBDA = ', lambda0
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Short circuit current in A:    |I0| = ', prabs
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Phase in deg.:              ARG(I0) = ', prphase
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Probe length in m:            prlen = ', prlen
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Probe radius in m:            prrad = ', prrad
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Location of the excit. in m:      X = ', prx0
  WRITE (FILEOUT,'(1X,A,E12.5)') '                                  Y = ', pry0
  WRITE (FILEOUT,'(1X,A,E12.5)') '                                  Z = ', prz0
  unit_vec = 0.0
  unit_vec(prdir) = 1.0
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Positive feed direction:          X = ', unit_vec(1)
  WRITE (FILEOUT,'(1X,A,E12.5)') '                                  Y = ', unit_vec(2)
  WRITE (FILEOUT,'(1X,A,E12.5)') '                                  Z = ', unit_vec(3)

END SUBROUTINE CBAA_MAKE_CURRENTPROBE_BVECTOR
!*******************************************************************************

SUBROUTINE CBAA_MAKE_PLANEWAVE_BVECTOR(t_pw,p_pw,e_pwE,mag_pwE,phase_pw)
  USE CBAA_data
  USE feminterface, ONLY: LOCAL_TO_GLOBAL_INDEX_TRI
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Adds to the b-vector an incident plane wave exitation in the CBAA case.
!*******************************************************************************
  REAL(SP), INTENT(IN) :: t_pw,p_pw,e_pwE,mag_pwE,phase_pw ! properties of the incident plane wave     
                                                           ! theta,phi,eta and mag/phase(radians) at the origin

  INTEGER(I4B), DIMENSION(3) :: tempedges      ! temp storage
  LOGICAL(LGT) :: test_free1                   ! does a face have associated dofs?
  REAL(SP) :: mag_pwH,e_pwH                    ! H-field quantities
  REAL(SP), DIMENSION(25,3) :: gauss_xyz       ! xyz co-ord.'s of triangular quadrature points
  REAL(SP), DIMENSION(25,8,3) :: gauss_t2      ! Values of ^z cross basis func.s at <gauss_xyz>
  COMPLEX(SPC), DIMENSION(3) :: H_inc          ! Incident H-value
  COMPLEX(SPC), DIMENSION(8) :: bvec_temp      ! max. 8 functions in a face
  INTEGER(I4B) :: bs_dimension                 ! number of basis functions in a face
  INTEGER(I4B) :: row,ielem,jj                 ! Counters
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)    ! The complex constant j
  INTEGER(I4B) :: num_qpoints, qcount          ! number of quadrature points for triangular surface integration
  INTEGER(I4B) :: face_order                   ! hierarchal order of current face
  REAL(SP), DIMENSION(25,ELEM_TRI_MATRIX_SIZE)   :: temp_quad_term1
  REAL(SP), DIMENSION(25,ELEM_TRI_MATRIX_SIZE,3) :: temp_quad_term2

  ! Choose the number of quadrature points:
  num_qpoints = 7

  ! Calculate the equivalent H-field desribing the plane wave:
  e_pwH = e_pwE + PI/2.0
  mag_pwH = mag_pwE/Z_zero

  Element_bvector_loop: DO ielem = 1,num_apelements
    ! Check if the facet does contain any dofs, and cycle to next aperture 
    ! element if this one is constrained:
    tempedges = which_local_edges(ielem,1:3)
    test_free1 = &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(1)))%free.OR. &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(2)))%free.OR. &
      edges(elements(ap_elnumbers(ielem))%edges(tempedges(3)))%free.OR. &
      faces(elements(ap_elnumbers(ielem))%faces(which_local_face(ielem)))%free
    IF (.NOT.test_free1) CYCLE Element_bvector_loop

    ! Find the current aperture face's hierarchal order:
    face_order = MAX_ORDER(ap_elnumbers(ielem),which_local_face(ielem))
    
    ! Choose the dimension of the vector {b^s} (see Jin, eq.9.137):
    SELECT CASE (face_order)
    CASE (1) !CT/LN
      bs_dimension = 3 ! 3 E1 in the face
    CASE (2) ! LT/QN
      bs_dimension = 8 ! 3 E1, 3 E2, 1 F1 and 1 F2
    END SELECT

    ! Calculate the xyz-quadrature points and the basis function values 
    ! at these points:
    SELECT CASE (face_order)
    CASE (1) ! CT/LN
      CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,3, &
	                                    gauss_xyz(1:num_qpoints,1:3),temp_quad_term1(1:num_qpoints,1:3), &
									    temp_quad_term2(1:num_qpoints,1:3,1:3))
      gauss_t2(1:num_qpoints,1:3,1:3) = temp_quad_term2(1:num_qpoints,1:3,1:3)
    CASE (2) ! LT/QN
      CALL EVALUATE_FACE_FUNCTIONS_QUAD(ap_elnumbers(ielem),which_local_face(ielem),num_qpoints,8, &
	                                    gauss_xyz(1:num_qpoints,1:3),temp_quad_term1(1:num_qpoints,1:8), &
									    temp_quad_term2(1:num_qpoints,1:8,1:3))
      gauss_t2(1:num_qpoints,1:8,1:3) = temp_quad_term2(1:num_qpoints,1:8,1:3)
	END SELECT
       
    ! Now perform the integral over the face area (see Jin, eq. 9.137):
    bvec_temp = (0.0,0.0) ! initialize
    DO qcount = 1,num_qpoints

      ! Calculate incident field value at this quadrature point:
      H_inc(1) = -COS(e_pwH)*COS(t_pw)*COS(p_pw)-SIN(e_pwH)*SIN(p_pw)
      H_inc(2) = -COS(e_pwH)*COS(t_pw)*SIN(p_pw)+SIN(e_pwH)*COS(p_pw)
      H_inc(3) = COS(e_pwH)*SIN(t_pw)
      H_inc = EXP(cj*k0*(  SIN(t_pw)*COS(p_pw)*gauss_xyz(qcount,1)    &
                         + SIN(t_pw)*SIN(p_pw)*gauss_xyz(qcount,2)    &
                         + COS(t_pw)*gauss_xyz(qcount,3)  )) * H_inc
      H_inc = mag_pwH * EXP(cj*phase_pw) * H_inc

      ! Add to the total for this element:
      bvec_temp(1:bs_dimension) = bvec_temp(1:bs_dimension)                 &
                            + H_inc(1)*gauss_t2(qcount,1:bs_dimension,1)    &
                            + H_inc(2)*gauss_t2(qcount,1:bs_dimension,2)    &
                            + H_inc(3)*gauss_t2(qcount,1:bs_dimension,3)
    END DO

    ! Now add to the b-vector:
    DO jj = 1,bs_dimension
      row = LOCAL_TO_GLOBAL_INDEX_TRI(ap_elnumbers(ielem),which_local_face(ielem),jj)
      IF (row.GT.0) THEN ! check that this is a dof (free variable)
        b_vec_c(row) = b_vec_c(row) - 2.0*cj*k0*Z_zero*bvec_temp(jj) 
      END IF
    END DO

  END DO Element_bvector_loop

  ! Write the exitation information to the output file:
  WRITE (FILEOUT,'(//,20X,A,/)') 'EXCITATION BY PLANE LINEAR POLARISED ELECTROMAGNETIC WAVE'
  WRITE (FILEOUT,'(1X,A,I4)')    'Number of excitation:             N = ', 0
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Frequency in Hz:               FREQ = ', frequency
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Wavelength in m:             LAMBDA = ', lambda0
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Direction of incidence:       THETA = ', t_pw*180.0/PI
  WRITE (FILEOUT,'(1X,A,E12.5)') '                                PHI = ', p_pw*180.0/PI
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Direction of polarization:      ETA = ', e_pwE*180.0/PI
  WRITE (FILEOUT,'(1X,A,E12.5)') 'Direction of propagation:    BETA0X = ', -k0*SIN(t_pw)*COS(p_pw)
  WRITE (FILEOUT,'(1X,A,E12.5)') '                             BETA0Y = ', -k0*SIN(t_pw)*SIN(p_pw)
  WRITE (FILEOUT,'(1X,A,E12.5)') '                             BETA0Z = ', -k0*COS(t_pw)
  WRITE (FILEOUT,'(1X,A,E12.5,A,E12.5)') 'Field strength in V/m:        |E0X| = ', 0.0,'   ARG(E0X) = ',0.0
  WRITE (FILEOUT,'(1X,A,E12.5,A,E12.5)') '(Phase in deg.)               |E0Y| = ', 0.0,'   ARG(E0X) = ',0.0
  WRITE (FILEOUT,'(1X,A,E12.5,A,E12.5)') '                              |E0Z| = ', 0.0,'   ARG(E0X) = ',0.0

END SUBROUTINE CBAA_MAKE_PLANEWAVE_BVECTOR
!*******************************************************************************

SUBROUTINE CBAA_MATRIX_ALLOCATE
  USE CBAA_data
  USE feminterface, ONLY: MATRIX_SPARSE_ALLOCATE
  USE frequency_data
  USE geometry
  USE matrix
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! MMB 14 Dec 2000:
! Allocate the appropriate storage for the CBAA system matrix according to
! the Hierarchal Order, FMM, Sparse flags.
!*******************************************************************************

  IF (.NOT.SPARSE) THEN
    ALLOCATE(A_mat_c(dof,dof)) ! Full storage for the system matrix
  ELSE
    CALL MATRIX_SPARSE_ALLOCATE ! Allocate sparse FE storage
    
	IF (ap_dof.GT.0) THEN

	  IF (CBAA_FMM_storage) THEN
        k0_fmm = (2.0*PI/c_0)*FRdata%freq0 ! use the first frequency value as FMM allocation parameter.
        CALL CBAA_FMM_ALLOCATE ! (FMM is currently only for CT/LN elements.)
        IF (CBAA_FMM_debug) THEN 
          ! FMM debug mode: both FMM and standard BI matrices are stored.
          ALLOCATE(CBAA_BE_mat(ap_dof,ap_dof))
          ALLOCATE(temp_BE_mat(ap_dof,ap_dof))
        END IF
      ELSE
        ALLOCATE(CBAA_BE_mat(ap_dof,ap_dof)) ! Full BI matrix storage
      END IF 

    END IF
  END IF

  ALLOCATE(x_vec_c(dof)) ! Storage for the solution vector        
  ALLOCATE(b_vec_c(dof)) ! Storage for the exitation vector

END SUBROUTINE CBAA_MATRIX_ALLOCATE
!*******************************************************************************

SUBROUTINE CBAA_PREPROCESSING
  USE CBAA_data
  USE geometry
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! - Set CBAA FMM flags.
! - Find the max aperture dimensions.
! 22 Jan 2001 MMB
!*******************************************************************************
  INTEGER(I4B) :: iface,ielem,iedge               ! counters
  INTEGER(I4B), DIMENSION(3) :: tempedges         ! temp storage
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes     ! temp storage
  REAL(SP), DIMENSION(3) :: tempvec,tempvec2      ! temp storage
  REAL(SP), DIMENSION(3) :: pl_norm               ! basis vectors for plane's geometry
  REAL(SP), DIMENSION(4) :: pl_eq                 ! equation for plane: ax+by+cz+d=0
  REAL(SP) :: minx, miny, maxx, maxy              ! Used in determining aperture dimensions                 
  LOGICAL(LGT) :: first_face                      ! flag in aperture dims search

  ! Storage and matrix solution:
  CBAA_FMM_storage = .FALSE.
  CBAA_FMM_debug = .FALSE.
  first_face = .TRUE.

  Element_loop: DO ielem = 1,num_elements
    Face_loop: DO iface = 1,4
      Aperture_Test: IF (faces(elements(ielem)%faces(iface))%CBAA_aperture) THEN
          
        ! find the global numbers of the face's 3 nodes
        tempfacenodes = faces(elements(ielem)%faces(iface))%nodes

        ! Set up aperture information necessary for FMM:
        IF (CBAA_FMM_storage) THEN
          minx = MIN(vertices(tempfacenodes(1))%coord(1),vertices(tempfacenodes(2))%coord(1), &
                     vertices(tempfacenodes(3))%coord(1))
          miny = MIN(vertices(tempfacenodes(1))%coord(2),vertices(tempfacenodes(2))%coord(2), &
                     vertices(tempfacenodes(3))%coord(2))
          maxx = MAX(vertices(tempfacenodes(1))%coord(1),vertices(tempfacenodes(2))%coord(1), &
                     vertices(tempfacenodes(3))%coord(1))
          maxy = MAX(vertices(tempfacenodes(1))%coord(2),vertices(tempfacenodes(2))%coord(2), &
                     vertices(tempfacenodes(3))%coord(2))
          IF (first_face) THEN
            ap_dims = (/ minx, miny, maxx, maxy /)
          ELSE
            IF (minx.LT.ap_dims(1)) ap_dims(1) = minx
            IF (miny.LT.ap_dims(2)) ap_dims(2) = miny
            IF (maxx.GT.ap_dims(3)) ap_dims(3) = maxx
            IF (maxy.GT.ap_dims(4)) ap_dims(4) = maxy
          END IF
        END IF
        first_face = .FALSE.

      END IF Aperture_Test
    END DO Face_loop
  END DO Element_loop

  IF (CBAA_FMM_storage) PRINT *, 'Max. aperture dimensions:',ap_dims
  print *, 'num_apelements = ',num_apelements

END SUBROUTINE CBAA_PREPROCESSING
!*******************************************************************************


SUBROUTINE CBAA_COUNT_APERTURE_DOFS
  USE CBAA_data
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Counts the number of aperture degrees of freedom.
! 
! 22 Jan 2001 MMB
! 2002-05-15: Added capabilty for up to QTQN elements. MMB.
!*******************************************************************************
  INTEGER(I4B) :: iface,iedge                     ! counters

  ! Count the degrees of freedom in the aperture:
  ap_dof = 0
  DO iedge = 1,num_edges
    IF (edges(iedge)%CBAA_aperture.AND.edges(iedge)%free) THEN
      SELECT CASE(edges(iedge)%order)
      CASE(1)
        IF (edges(iedge)%mixed) THEN
		  ap_dof = ap_dof + 1 ! CT/LN
		ELSE
		  ap_dof = ap_dof + 2 ! LT/LN
		END IF
      CASE(2)
        IF (edges(iedge)%mixed) THEN
		  ap_dof = ap_dof + 2 ! LT/QN
		ELSE
		  ap_dof = ap_dof + 3 ! QT/QN
		END IF
      CASE DEFAULT
        STOP 'IE: invalid order in CBAA_COUNT_APERTURE_DOFS.'
      END SELECT      
    END IF
  END DO
  DO iface = 1,num_faces
    IF (faces(iface)%CBAA_aperture.AND.faces(iface)%free) THEN
      SELECT CASE(faces(iface)%order)
      CASE(1) ! CT/LN and LT/LN
        ! No action required
      CASE(2) ! LT/QN
        IF (faces(iface)%mixed) THEN
		  ap_dof = ap_dof + 2 ! LT/QN
		ELSE
		  ap_dof = ap_dof + 3 ! QT/QN
		END IF
      CASE DEFAULT
        STOP 'IE: invalid order in CBAA_COUNT_APERTURE_DOFS.'
      END SELECT      
    END IF
  END DO

  print *, 'ap_dof = ',ap_dof

END SUBROUTINE CBAA_COUNT_APERTURE_DOFS
!*******************************************************************************


SUBROUTINE CBAA_PROBE_VOLTAGE_CALC(prx0,pry0,prz0,prdir,prlen,prabs,prphase,prvolt)
  USE frequency_data
  USE math_tools, ONLY: PHASE
  USE nrtype
  USE problem_info
  USE unit_numbers
  IMPLICIT NONE
!*******************************************************************************
! Integrates the E-field along the current probe's center to calculate the input 
! impedance.
! 21 Jan 2001 MMB
!*******************************************************************************
  REAL(SP), INTENT(IN) :: prx0,pry0,prz0  ! xyz of starting point
  INTEGER(I4B), INTENT(IN) :: prdir       ! x/y/z : 1/2/3 direction of probe
  REAL(SP), INTENT(IN) :: prlen           ! length of probe
  REAL(SP), INTENT(IN) :: prabs           ! Abs value of the current
  REAL(SP), INTENT(IN) :: prphase         ! phase of the current in radians
  COMPLEX(SPC), INTENT(OUT) :: prvolt     ! probe voltage

  INTEGER(I4B) :: iint                        ! counters
  REAL(SP) :: prx,pry,prz                     ! temp integration variables
  INTEGER(I4B) :: prdisc                      ! discretization for integration porposes
  COMPLEX(SPC), DIMENSION(3) :: E_xyz,H_xyz   ! arguments for CBAA_FIELDCALC
  COMPLEX(SPC) :: complex1                    ! General temporary variable

  ! Check that this is a CBAA analysis:
  IF (.NOT.(CBAA_ANALYSIS)) STOP 'IE: Probe voltage can only be calculated for CBAA'

  ! Calculate the number of integration points, with the following guidelines:
  ! - Use a discretization of lambda0/1000
  ! - At least 20 points in the lenght-wise integration
  prdisc = CEILING(1000.0*prlen/lambda0)
  IF (prdisc.LT.20) THEN
    prdisc = 20
  END IF

  prvolt = (0.0,0.0)  ! initialize
  prx = prx0          !     "
  pry = pry0          !     "
  prz = prz0          !     "

  ! Calculate integral: (eq. 9.147 - jin)
  Probe_lengthwise_integration: DO iint = 0,prdisc-1 ! Integrating along the length
    SELECT CASE(prdir) ! increment the integration position
    CASE(1) 
      prx = prx0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
    CASE(2) 
      pry = pry0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
    CASE(3) 
      prz = prz0 + prlen/(2.0*prdisc) + iint*prlen/prdisc
    END SELECT
    CALL CBAA_FIELDCALC(prx,pry,prz,E_xyz,H_xyz)
    prvolt = prvolt - E_xyz(prdir)*prlen/REAL(prdisc)
  END DO Probe_lengthwise_integration

  ! Write the exitation results to the output file:
  WRITE (FILEOUT,'(//,20X,A,I4,/)') 'DATA OF THE PROBE CURRENT SOURCE NO.',0
  WRITE (FILEOUT,'(20X,4(A13))') 'real part','imag. part','magn.','phase'
  ! Write out the voltage:
  WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Voltage','in V', REAL(prvolt), &
  AIMAG(prvolt), ABS(prvolt), PHASE(prvolt)
  ! Write out the admittance:
  complex1 = prabs*EXP(j*prphase)/prvolt
  WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Admitt.','in A/V', REAL(complex1), AIMAG(complex1), &
                                                ABS(complex1), PHASE(complex1)
  ! Write out the impedance:
  complex1 = 1.0/complex1     
  WRITE (FILEOUT,'(1X,A10,1X,A8,4(1X,E12.5))') 'Impedance','in Ohm', REAL(complex1), AIMAG(complex1), &
                                                ABS(complex1), PHASE(complex1)

END SUBROUTINE CBAA_PROBE_VOLTAGE_CALC
!*******************************************************************************

 
FUNCTION GLOBAL_GROUP_NUMBER(x_origin,y_origin,dim,numx,edgenum,facenum)
!*******************************************************************************
! Returns the global group number to which the global edge <edgenum> belongs. 
! Also recieves the number of groups in the x-direction and the group dimension
! as well as the origin of the group axes.
!*******************************************************************************
  USE CBAA_data
  USE geometry
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype         
  IMPLICIT NONE

  INTEGER(I4B), INTENT(IN) :: numx
  REAL(SP), INTENT(IN) :: x_origin,y_origin,dim   
  INTEGER(I4B), OPTIONAL, INTENT(IN) :: edgenum,facenum
  INTEGER(I4B) :: GLOBAL_GROUP_NUMBER
  REAL(SP), DIMENSION(2) :: xy_coord
  REAL(SP) :: gr_xco,gr_yco

  IF (PRESENT(edgenum)) THEN
    ! Place edge center co-ordinate in xy_coord:
    xy_coord = 0.5*( vertices(edges(edgenum)%nodes(1))%coord(1:2) + &
                     vertices(edges(edgenum)%nodes(2))%coord(1:2) )   
  END IF

  IF (PRESENT(facenum)) THEN
    ! Place face center co-ordinate in xy_coord:
    xy_coord = 0.33333333*( vertices(faces(facenum)%nodes(1))%coord(1:2) + &
                            vertices(faces(facenum)%nodes(2))%coord(1:2) + &
                            vertices(faces(facenum)%nodes(3))%coord(1:2) )   
  END IF

  ! Calculate the global group number:
  gr_xco = ABS((xy_coord(1)-x_origin)/dim)
  gr_yco = ABS((xy_coord(2)-y_origin)/dim)
  GLOBAL_GROUP_NUMBER = numx*FLOOR(gr_yco) + FLOOR(gr_xco) + 1
  
END FUNCTION GLOBAL_GROUP_NUMBER
!*******************************************************************************

FUNCTION NEAR_NEIGHBOURS(gr1,gr2,dim,numx)
!*******************************************************************************
! TRUE if the global group numbers gr1 and gr2 are near neighbours in the FMM.
!*******************************************************************************
  USE CBAA_data
  USE frequency_data
  USE matrix
  USE unit_numbers
  USE problem_info
  USE nrtype         
  IMPLICIT NONE

  INTEGER(I4B), INTENT(IN) :: gr1,gr2,numx
  REAL(SP), INTENT(IN) :: dim
  LOGICAL(LGT) :: NEAR_NEIGHBOURS
  REAL(SP) :: gr1gr2_dist
  INTEGER(I4B), DIMENSION(2) :: gr1co,gr2co

  gr1co(2) = CEILING(REAL(gr1)/REAL(numx))
  gr1co(1) = gr1 - numx*(gr1co(2)-1)
  gr2co(2) = CEILING(REAL(gr2)/REAL(numx))
  gr2co(1) = gr2 - numx*(gr2co(2)-1)
  gr1gr2_dist = dim*SQRT(REAL(SUM((gr2co-gr1co)**2)))

  IF ( FarInteract_min .GT. gr1gr2_dist ) THEN
    NEAR_NEIGHBOURS = .TRUE.
  ELSE
    NEAR_NEIGHBOURS = .FALSE.
  END IF

END FUNCTION NEAR_NEIGHBOURS
!*******************************************************************************


SUBROUTINE CBAA_BOUNDARY_SEARCH
  USE CBAA_data
  USE geometry
  USE nrtype
  IMPLICIT NONE
  !*******************************************************************************
  ! This routine searches through all the boundary faces for the edges/faces in 
  ! the CBAA aperture and flags them thus.
  !*******************************************************************************
  INTEGER(I4B) :: iface,ielem
  INTEGER(I4B), DIMENSION(3) :: tempedges, gfortran_tmp1, gfortran_tmp2
  REAL(SP), DIMENSION(3) :: z_tmp           ! z coordinates of nodes per face.
  
  ! Allocate storage for the aperture information:
  ALLOCATE(ap_elnumbers(num_elements))
  ALLOCATE(which_local_edges(num_elements,3))
  ALLOCATE(which_local_face(num_elements))
  
  num_apelements = 0
  
  ELEMENT_LOOP: DO ielem = 1,num_elements
     FACE_LOOP: DO iface = 1,4 ! Face-wise search to find edges in the quad.
      FACE_CONNECT_TEST: IF (elements(ielem)%connect2elem(iface).EQ.0) THEN     
         
         ! Check whether in z=0 plane:
         gfortran_tmp1 = (/1,2,3/)
         gfortran_tmp2 = GLOBAL_FACENODES(ielem,iface)
         z_tmp(gfortran_tmp1) = vertices(gfortran_tmp2)%coord(3)
         IF ((ABS(z_tmp(1)).LE.EPS).AND. &
              (ABS(z_tmp(2)).LE.EPS).AND. &
              (ABS(z_tmp(3)).LE.EPS)) THEN ! the face is in the CBAA aperture
            
            tempedges = LOCAL_FACEEDGES(iface)
            
            ! Flag the three edges and face as part of the CBAA aperture:
            edges(elements(ielem)%edges(tempedges))%CBAA_aperture = .TRUE.
            faces(elements(ielem)%faces(iface))%CBAA_aperture = .TRUE.
            
            ! Set up necessary aperture information:
            num_apelements = num_apelements + 1
            ap_elnumbers(num_apelements) = ielem
            tempedges = LOCAL_FACEEDGES(iface)
            which_local_edges(num_apelements,1:3) = tempedges
            which_local_face(num_apelements) = iface
            
         END IF
      END IF FACE_CONNECT_TEST
   END DO FACE_LOOP
END DO ELEMENT_LOOP

END SUBROUTINE CBAA_BOUNDARY_SEARCH
!*******************************************************************************


SUBROUTINE COAX_BOUNDARY_SEARCH(whitney_approx)
  USE coax_feed
  USE geometry
  USE math_tools, ONLY: FIND_IN_LIST, VECTOR_LENGTH
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine searches through all the boundary faces for the edges/faces in 
! the coax aperture(s) and flags them thus. It also records the global edge
! numbers of the radial coax edges.
! Note: Assumes 6 radial edges in every coax aperture.
!*******************************************************************************
  LOGICAL(LGT), INTENT(IN) :: whitney_approx      ! coax formulation to use

  INTEGER(I4B) :: iface,iedge,ielem,listpos,this_edge,iport
  REAL(SP), DIMENSION(3) :: pl_norm               ! basis vectors for plane's geometry
  REAL(SP), DIMENSION(4) :: pl_eq                 ! equation for plane: ax+by+cz+d=0
  REAL(SP), DIMENSION(3) :: tempvec,tempvec2      ! temp storage
  INTEGER(I4B), DIMENSION(3) :: tempedges,tempnodes
  REAL(SP), DIMENSION(3) :: node1_coord,node2_coord ! temp storage
  LOGICAL(LGT) :: in_aperture                     ! is face in coax aperture?

  ! Flag coax aperture edges and faces:
  DO iport = 1,NUM_AX_cards ! cycle through all AX cards

    ! Initialise:
    IF (whitney_approx) THEN
	  AXdata(iport)%coaxedges_num = 0
      AXdata(iport)%coaxedges     = 0 ! vector assignment
    END IF
       
    ELEMENT_LOOP: DO ielem = 1,num_elements
      FACE_LOOP: DO iface = 1,4
        
		! Cycle if face not on boundary:
		IF (elements(ielem)%connect2elem(iface).NE.0) CYCLE FACE_LOOP 

        ! Cycle if face not in coax aperture:
		CALL FACE_IN_CIRCLE(AXdata(iport)%normal,AXdata(iport)%coaxcentre,AXdata(iport)%coax_b,elements(ielem)%faces(iface),in_aperture)
		IF (.NOT.in_aperture) CYCLE FACE_LOOP

        tempedges = LOCAL_FACEEDGES(iface)
		faces(elements(ielem)%faces(iface))%coax_aperture     = .TRUE.
        edges(elements(ielem)%edges(tempedges))%coax_aperture = .TRUE.
		faces(elements(ielem)%faces(iface))%coaxnumber        = iport
        edges(elements(ielem)%edges(tempedges))%coaxnumber    = iport
                
        ! If the face is in the center conductor, then it must be 
		! a PEC boundary:
		CALL FACE_IN_CIRCLE(AXdata(iport)%normal,AXdata(iport)%coaxcentre,AXdata(iport)%coax_a,elements(ielem)%faces(iface),in_aperture)
		IF (in_aperture) THEN
		  faces(elements(ielem)%faces(iface))%PEC     = .TRUE.
          edges(elements(ielem)%edges(tempedges))%PEC = .TRUE.
		END IF

        IF (whitney_approx) THEN

          DO iedge = 1,3
            this_edge   = elements(ielem)%edges(tempedges(iedge))
            node1_coord = vertices(edges(this_edge)%nodes(1))%coord(1:3)
            node2_coord = vertices(edges(this_edge)%nodes(2))%coord(1:3)
                
            IF ((VECTOR_LENGTH(node1_coord-AXdata(iport)%coaxcentre).LT.EPS).OR. &
                (VECTOR_LENGTH(node2_coord-AXdata(iport)%coaxcentre).LT.EPS)) THEN ! edge does connect to coax center

              CALL FIND_IN_LIST(AXdata(iport)%coaxedges,6,this_edge,listpos)
              IF (listpos.EQ.0) THEN ! edge number not yet recorded
                AXdata(iport)%coaxedges_num = AXdata(iport)%coaxedges_num + 1
                AXdata(iport)%coaxedges(AXdata(iport)%coaxedges_num) = &
                  this_edge ! record glabal edge numbers of aperture, radial edges
                print *, 'another radial coax aperture edge...'
              END IF
            END IF
          END DO
        END IF
      END DO FACE_LOOP
    END DO ELEMENT_LOOP

    ! Find centre conductor edges and constrain them as Etan = 0:
    DO iedge = 1,num_edges
      tempvec  = vertices(edges(iedge)%nodes(1))%coord(1:3)
      tempvec2 = vertices(edges(iedge)%nodes(2))%coord(1:3)
          
      IF ((VECTOR_LENGTH(tempvec-AXdata(iport)%coaxcentre) + &
          VECTOR_LENGTH(tempvec-AXdata(iport)%coaxend)     - &
          ABS(AXdata(iport)%coaxlen)).LT.EPS) THEN ! first node is on centre conductor
            
        IF ((VECTOR_LENGTH(tempvec2-AXdata(iport)%coaxcentre) + &
          VECTOR_LENGTH(tempvec2-AXdata(iport)%coaxend)      - &
          ABS(AXdata(iport)%coaxlen)).LT.EPS) THEN ! second node is on centre conductor

          edges(iedge)%PEC = .TRUE.
          print *, 'another constrained edge...'

        END IF
      END IF
    END DO

  END DO ! AX card loop

END SUBROUTINE COAX_BOUNDARY_SEARCH
!*******************************************************************************


SUBROUTINE CBAA_MAKE_BI_MATRIX_ELEMENTAL(elem_outer,local_face_outer,num_vbfs_outer,num_qpts_outer, &
                                         elem_inner,local_face_inner,num_vbfs_inner,num_qpts_inner, &
										 BI_mat)
  USE frequency_data
  USE nrtype
  USE problem_info
  IMPLICIT NONE
!*******************************************************************************
! This routine returns the contribution of elemental, facial pair designated by
! <elem_outer,local_face_outer> and <elem_inner,local_face_inner>, to the CBAA
! BI A-matrix contribution as <BI_matrix>. The outer and inner numbers of quadrature
! points must also be specified. Triangular surface quadrature is always used for 
! the outer integral. In case of the self-term, the inner number of 
! quadrature points will designate the number of radial, trapezoidal integration 
! points for cylindrical co-ordinate integration, otherwise it designates the number
! of points for triangular surface quadrature (same as outer intergral).
!
! 2002-05-13: Created from old CBAA_MAKE_BE_AMATRIX. MMB.
!*******************************************************************************
  SAVE
  INTEGER(I4B), INTENT(IN)  :: elem_outer,local_face_outer, &
                               num_vbfs_outer,num_qpts_outer ! Outer integral spec.'s
  INTEGER(I4B), INTENT(IN)  :: elem_inner,local_face_inner, &
                               num_vbfs_inner,num_qpts_inner ! Inner integral spec.'s
  COMPLEX(SPC), DIMENSION(num_vbfs_outer,num_vbfs_inner),   &
                INTENT(OUT) :: BI_mat                        ! Elemental facial pair BI contribution matrix

  INTEGER(I4B) :: old_elem_outer = 0
  INTEGER(I4B) :: old_local_face_outer,old_num_vbfs_outer, & 
				  old_num_qpts_outer
  LOGICAL(LGT) :: new_eval_outer
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,3)                      :: quad_xyz_outer
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE)   :: quad_term1_outer
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE,3) :: quad_term2_outer
  INTEGER(I4B) :: iquad,mm,nn ! counters
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE)   :: term1_inner
  COMPLEX(SPC), DIMENSION(ELEM_TRI_MATRIX_SIZE,3) :: term2_inner

  ! Check if new outer elemental face data should be created:
  new_eval_outer = (old_elem_outer      .NE.elem_outer)      .OR. &
                   (old_local_face_outer.NE.local_face_outer).OR. &
			       (old_num_vbfs_outer  .NE.num_vbfs_outer)  .OR. &
			       (old_num_qpts_outer  .NE.num_qpts_outer)

  ! If this is a new outer element, create new tables, containing 
  ! the following values for this element, so that they are calculated 
  ! only once for every set of contiguous equal outer element number calls:
  ! - xyz's of quadrature points in this element, <quad_xyz_outer>
  ! - del dot ^z cross (face VBFs) at the integration points, <quad_term1_outer>.
  ! - ^z cross (face VBFs) at the integration points, <quad_term2_outer>.
  !   (Both weighted with face area and quadrature point weights.)
  NEW_ELEMENT_OUTER: IF (new_eval_outer) THEN
    
	! Save for future use:
	old_elem_outer       = elem_outer
    old_local_face_outer = local_face_outer
    old_num_vbfs_outer   = num_vbfs_outer
    old_num_qpts_outer   = num_qpts_outer

    ! Evaluate the quadrature points and the functions at the points:
    ! (Calculated here, so that it is not recalculated for every 
	! outer(constant):inner combination.)
    CALL EVALUATE_FACE_FUNCTIONS_QUAD(elem_outer,local_face_outer,num_qpts_outer,num_vbfs_outer, &
		                              quad_xyz_outer  (1:num_qpts_outer,1:3),                    &
								      quad_term1_outer(1:num_qpts_outer,1:num_vbfs_outer),       &
		                              quad_term2_outer(1:num_qpts_outer,1:num_vbfs_outer,1:3))
  END IF NEW_ELEMENT_OUTER

  ! Initialise the elemental BI matrix:
  BI_mat = (0.0,0.0)
		 
  ! Perform Quadrature over the outer element surface:
  QUAD_LOOP: DO iquad = 1,num_qpts_outer

    ! Carry out the inner integral:
    IF ((elem_outer.NE.elem_inner).OR.(local_face_outer.NE.local_face_inner)) THEN
      CALL INNER_INTEGRAL_NORMAL(quad_xyz_outer(iquad,1:3),elem_inner,local_face_inner,num_vbfs_inner, &
	                             num_qpts_inner,term1_inner(1:num_vbfs_inner),term2_inner(1:num_vbfs_inner,1:3))
    ELSE
      CALL INNER_INTEGRAL_SELF(quad_xyz_outer(iquad,1:3),elem_outer,elem_inner,local_face_inner,num_vbfs_inner, &
	                           60,term1_inner(1:num_vbfs_inner),term2_inner(1:num_vbfs_inner,1:3))
    END IF

    ! Add the current outer quadrature point's contribution to the elemental matrix:
    DO mm = 1,num_vbfs_outer
      DO nn = 1,num_vbfs_inner
        
		! Add the term1 contribution:
		BI_mat(mm,nn) =   BI_mat(mm,nn) &
		                + 2.0*quad_term1_outer(iquad,mm)*term1_inner(nn)

		! Add the term2 contribution:
		BI_mat(mm,nn) =   BI_mat(mm,nn) &
		                - 2.0*(k0**2)*SUM(quad_term2_outer(iquad,mm,1:3)*term2_inner(nn,1:3))

      END DO
    END DO
  END DO QUAD_LOOP

END SUBROUTINE CBAA_MAKE_BI_MATRIX_ELEMENTAL
!*******************************************************************************


SUBROUTINE INNER_INTEGRAL_NORMAL(r_vec,elem,local_face,num_vbfs,num_qpts,inner_term1,inner_term2)
  USE frequency_data
  USE geometry
  USE nrtype
  USE problem_info
  USE quad_tables
  IMPLICIT NONE
!*******************************************************************************
! Calculates the inner integral values (<inner_t1>,<inner_t2>) for a specific 
! outer variable (vector <r>), when the two elements are not the same. The 
! integral is evaluated by gaussian quadratrue. <elem> is the index into <ap_elnumbers>.
! NOTE: t1 and t2 refer to the two terms in the [P^{st}] sub-matrix definition.
!  
! Updated: The inner BI integral is evaluated fir glbal element <elem>'s local face <local_face>,
! with the outer integral's integration point <r_vec> outside the inner integration domain. The 
! number of VBF's on the face as well as the triangular quadrature rule to be used, are
! specified with <num_vbfs> and <num_qpts>.
!
! 2002-05-12: Extensive revamp to generalize for easy addition of new VBFs. MMB.
!*******************************************************************************
  SAVE
  REAL(SP),DIMENSION(3), INTENT(IN) :: r_vec                      ! cartesian co-ord. of outer integral gausspoint
  INTEGER(I4B), INTENT(IN) :: elem,local_face                     ! inner element global number and local face
  INTEGER(I4B), INTENT(IN) :: num_vbfs,num_qpts                   ! Num. quadrature points to use + num VBFs to evaluate.
  COMPLEX(SPC), DIMENSION(num_vbfs), INTENT(OUT) :: inner_term1   ! inner integral of term 1 in Jin, eq.9.136, for E1
  COMPLEX(SPC), DIMENSION(num_vbfs,3), INTENT(OUT) :: inner_term2 ! inner integral of term 2 in Jin, eq.9.136, for E1

  INTEGER(I4B) :: old_elem = 0
  INTEGER(I4B) :: old_local_face,old_num_vbfs,old_num_qpts
  LOGICAL(LGT) :: new_eval
  REAL(SP) :: abs_R
  COMPLEX(SPC) :: Green_val
  INTEGER(I4B) :: gcount
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,3)                      :: quad_xyz
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE)   :: quad_term1
  REAL(SP), DIMENSION(MAX_QUAD_TRI_RULES,ELEM_TRI_MATRIX_SIZE,3) :: quad_term2

  ! Check if new elemental face data should be created:
  new_eval = (old_elem.NE.elem)            .OR. &
             (old_local_face.NE.local_face).OR. &
			 (old_num_vbfs.NE.num_vbfs)    .OR. &
			 (old_num_qpts.NE.num_qpts)

  ! First, create new tables, containing the following values for this element,
  ! so that they are calculated only once for every outer element number:
  ! - xyz's of quadrature points in this element, <quad_xyz>
  ! - del dot ^z cross (face VBFs) at the integration points, <quad_term1>.
  ! - ^z cross (face VBFs) at the integration points, <quad_term2>.
  !   (Both weighted with face area and quadrature point weights.)
  NEW_ELEMENT: IF (new_eval) THEN
    
	! Save for future use:
	old_elem       = elem
    old_local_face = local_face
    old_num_vbfs   = num_vbfs
    old_num_qpts   = num_qpts

    ! Evaluate the quadrature points and the functions at the points:
	CALL EVALUATE_FACE_FUNCTIONS_QUAD(elem,local_face,num_qpts,num_vbfs,quad_xyz(1:num_qpts,1:3), &
		                              quad_term1(1:num_qpts,1:num_vbfs),quad_term2(1:num_qpts,1:num_vbfs,1:3))
    ! Scaling:
    quad_term1(1:num_qpts,1:num_vbfs)     = (1.0/(4.0*PI)) * quad_term1(1:num_qpts,1:num_vbfs)
    quad_term2(1:num_qpts,1:num_vbfs,1:3) = (1.0/(4.0*PI)) * quad_term2(1:num_qpts,1:num_vbfs,1:3)
  END IF NEW_ELEMENT

  ! Integrate over the triangle:
  inner_term1 = (0.0,0.0) ! initialise
  inner_term2 = (0.0,0.0) ! initialise
  DO gcount = 1,num_qpts
    abs_R = SQRT(SUM((r_vec - quad_xyz(gcount,1:3))**2)) ! Calculate ABS(^r - ^r')
    Green_val = EXP(-j*k0*abs_R)/abs_R                   ! 4*PI*(value of Green function) at current gausspoint     
    inner_term1 = inner_term1 + Green_val * quad_term1(gcount,1:num_vbfs)
    inner_term2 = inner_term2 + Green_val * quad_term2(gcount,1:num_vbfs,1:3)
  END DO
 
END SUBROUTINE INNER_INTEGRAL_NORMAL
!*******************************************************************************


SUBROUTINE INNER_INTEGRAL_SELF(r_vec,outer_elem,elem,local_face,num_vbfs,num_qpts,inner_term1,inner_term2)
  USE basis_function, ONLY: POINT_EVALUATE_FACE_FUNCTIONS
  USE frequency_data
  USE geometry
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Calculates the inner integral values (<inner_t1>,<inner_t2>) for a specific 
! outer variable (vector <r>). The integral is first transformed to the cilidrical
! co-ordinate system before it is evaluated: analytical in the r-direction and
! trapezoidal in the theta-direction. Takes the self-term singularity into account.
! The radial integral is carried out to polynomial order 1 for the (del dot ^z times VBF)
! terms and to polynomial order 2 for the (^z times VBF) terms.
!
! 2000-??-??: Created. MMB.
! 2000-??-??: Addded LT/QN capability. MMB.
! 2002-05-12: Generalized for arbitrary basis functions. MMB.
!*******************************************************************************
  REAL(SP),DIMENSION(3), INTENT(IN) :: r_vec                       ! cartesian co-ord. of outer integral gausspoint
  INTEGER(I4B), INTENT(IN) :: outer_elem,elem,local_face           ! outer&inner element global numbers and inner local face
  INTEGER(I4B), INTENT(IN) :: num_vbfs                             ! number of VBFs to integrate (starting at the lowest order)
  INTEGER(I4B), INTENT(IN) :: num_qpts                             ! number of angular trapezoidal integration points
  COMPLEX(SPC), DIMENSION(num_vbfs), INTENT(OUT) ::  inner_term1   ! inner integral of term 1 in Jin, eq.9.136
  COMPLEX(SPC), DIMENSION(num_vbfs,3), INTENT(OUT) ::  inner_term2 ! inner integral of term 2 in Jin, eq.9.136

  INTEGER(I4B) :: count1,count2,info,jj                 ! counters/miscellaneous
  INTEGER(I4B) :: theta_count,t_points                  ! counter and radial discretazation
  REAL(SP) :: theta_min,dtheta,theta,r1,r2,t1,t2,rtemp  ! Theta related...
  LOGICAL, DIMENSION(4) :: vertical                     ! whether lines are vertical
  REAL(SP), DIMENSION(4) :: vert_pos                    ! x = ? for vertical lines
  REAL(SP), DIMENSION(4,2) :: equations                 ! equation-coeff.'s for lines
  LOGICAL(LGT), DIMENSION(3) :: parallel                ! are lines parallel to int_line
  INTEGER(I4B) :: num_crosses                           ! Has the values: 2 or 3
  INTEGER(I4B) :: first_cross                           ! Check if it is the first crossing examined
  REAL(SP), DIMENSION(3,2) :: crosses                   ! Crossing points: int_line <-> edges
  INTEGER(I4B), DIMENSION(2) :: int_x                   ! index in crosses
  REAL(SP), DIMENSION(4) :: cart_ct                     ! Temp. store coordinates
  INTEGER(I4B), DIMENSION(3) :: tempedges               ! temp storage
  INTEGER(I4B), DIMENSION(2) :: tempnodes               ! temp storage
  COMPLEX(SPC), PARAMETER :: cj = (0.0,1.0)             ! The complex constant j
  INTEGER(I4B) :: num_r_points,rcount                   ! Variables for trapezoidal r-integration:
  REAL(SP) :: delta_r,radius                            ! Variables for trapezoidal r-integration:
  REAL(SP), DIMENSION(3,3) :: rmat                      ! radius matrix for A,B,C calculation
  REAL(SP), DIMENSION(3) :: fvec                        ! function vector for A,B,C calculation
  REAL(SP), DIMENSION(4,4) :: vmat2,simat2              ! simplex related matrices
  REAL(SP),    DIMENSION (:),ALLOCATABLE :: work_local  ! for <rmat> inversion
  INTEGER(I4B),DIMENSION (:),ALLOCATABLE :: ipiv_local  ! ditto.
  REAL(SP), DIMENSION(num_vbfs)   :: point1_term1,point2_term1,point3_term1
  REAL(SP), DIMENSION(num_vbfs)   :: m_term1,c_term1
  REAL(SP), DIMENSION(num_vbfs,3) :: point1_term2,point2_term2,point3_term2
  REAL(SP), DIMENSION(num_vbfs,3) :: A_term2,B_term2,C_term2

  ! Calculate the elemental matrices for co-ordinate evaluation:
  CALL SIMPLEX_COEFFICIENTS(elem,simat2(1:4,4),simat2(1:4,1),simat2(1:4,2),simat2(1:4,3),vmat2)

  ! for <rmat> inversion:
  ALLOCATE(ipiv_local(3))
  ALLOCATE(work_local(3))

  ! Gather edge information and theta_min and dtheta:
  vertical  = .FALSE.
  vert_pos  = 0.0
  equations = 0.0
  theta_min = 0.0
  dtheta    = 0.0
  tempedges = LOCAL_FACEEDGES(local_face)
  DO jj = 1,3

    tempnodes = LOCAL_EDGENODES(tempedges(jj))
    IF (ABS(vmat2(1,tempnodes(1))-vmat2(1,tempnodes(2))).LT.EPS) THEN
      vertical(jj) = .true.
      vert_pos(jj) = vmat2(1,tempnodes(1)) - r_vec(1)
    ELSE
      equations(jj,1) = (vmat2(2,tempnodes(1))-vmat2(2,tempnodes(2)))       &
                         /(vmat2(1,tempnodes(1))-vmat2(1,tempnodes(2)))
      equations(jj,2) = ((vmat2(2,tempnodes(1))-r_vec(2))*               &
                         (vmat2(1,tempnodes(2))-r_vec(1)) -              &
                         (vmat2(2,tempnodes(2))-r_vec(2))*               &
                         (vmat2(1,tempnodes(1))-r_vec(1)))               &
                         /(vmat2(1,tempnodes(2)) - vmat2(1,tempnodes(1)))
    END IF
  
    ! Find theta_min and dtheta:
    IF (outer_elem.EQ.elem) THEN
      theta_min = 0.0
      dtheta    = 2.0*PI
    ELSE
      r1 = SQRT((vmat2(1,tempnodes(1))-r_vec(1))**2 +               &
                (vmat2(2,tempnodes(1))-r_vec(2))**2)
      r2 = SQRT((vmat2(1,tempnodes(2))-r_vec(1))**2 +               &
                (vmat2(2,tempnodes(2))-r_vec(2))**2)
      IF (MIN(r1,r2).GT.EPS) THEN  ! We MAY divide by r1 and r2
  
        t1 = ACOS((vmat2(1,tempnodes(1))-r_vec(1))/r1)
        IF ((vmat2(2,tempnodes(1))-r_vec(2)).LT.0.0) THEN
          t1 = 2*pi - t1
        END IF
        t2 = ACOS((vmat2(1,tempnodes(2))-r_vec(1))/r2)
        IF ((vmat2(2,tempnodes(2))-r_vec(2)).LT.0.0) THEN
          t2 = 2*pi - t2
        END IF
        ! 0 < t1,t2 < 2*pi
  
        theta = MIN(ABS(t1-t2),2*pi-ABS(t1-t2))
        IF (theta.GT.dtheta) THEN
          dtheta = theta
          IF (ABS(t1-t2).GE.pi) THEN  ! taking jump in angle at 0 in account
            theta_min = MAX(t1,t2)
          ELSE
            theta_min = MIN(t1,t2)
          END IF
        END IF
   
      END IF  
    END IF 
  END DO

  ! Assign the number of directions to integrate analytically:
  ! (I.e. the number of angular trapezoidal quadrature points.)
  t_points = num_qpts

  ! Initialize the result variables:
  inner_term1 = (0.0,0.0)
  inner_term2 = (0.0,0.0)

  ! Start the summation over theta:
  theta_loop: DO theta_count = 1,t_points
    theta = theta_min + dtheta/(2.0_SP*(t_points)) + &
              (REAL(theta_count,SP)-1.0_SP)*(dtheta/t_points)

    ! Set up line information
    IF (ABS(COS(theta)).LT.EPS) THEN
      vertical(4) = .true.
      vert_pos(4) = 0.0
    ELSE
      equations(4,1) = TAN(theta)
      ! equations(4,2) is already = 0
    END IF
    parallel = .false.
    DO jj = 1,3
      IF (((ABS(equations(4,1)-equations(jj,1)).LT.EPS).AND.  &
         (.NOT.vertical(jj)).AND.(.NOT.vertical(4))).OR.      &
         (vertical(jj).AND.vertical(4))) THEN
        parallel(jj) = .true.
      END IF
    END DO
            
    ! Find the crossing points with the edge-lines:
    num_crosses = 0
    crosses     = 0.0
    DO jj = 1,3
      IF (.NOT.parallel(jj)) THEN
        num_crosses = num_crosses + 1
        IF (vertical(jj)) THEN ! Vertical edge
          crosses(num_crosses,1) = vert_pos(jj)
          crosses(num_crosses,2) = equations(4,1)*vert_pos(jj)
        ELSE IF (vertical(4)) THEN ! Vertical int_line
          crosses(num_crosses,1) = 0.0_SP
          crosses(num_crosses,2) = equations(jj,2)
        ELSE ! General case, both lines are non-vertical
          crosses(num_crosses,1) = (equations(4,2)-equations(jj,2))     &
                                   /(equations(jj,1)-equations(4,1))
          crosses(num_crosses,2) = (equations(4,1)*equations(jj,2)-     &
                                   equations(jj,1)*equations(4,2))      &
                                   /(equations(4,1)-equations(jj,1))
        END IF
      END IF
    END DO
              
    ! Now find the integration boudries: r1 and r2, with r1 < r2
    ! int_x points to the positions of r1,r2 in crosses
    r1          = 0.0
    r2          = 0.0
    int_x       = 1
    first_cross = 0
    DO jj = 1,num_crosses
      SELECT CASE(num_crosses)
      CASE(2)
        rtemp = SQRT(SUM(crosses(jj,1:2)**2))
        IF (rtemp.GT.r2) THEN
          r1 = r2
          int_x(1) = int_x(2)
          r2 = rtemp
          int_x(2) = jj
        ELSE 
          r1 = rtemp
          int_x(1) = jj
        END IF
      CASE(3)  ! First check whether crossing is on the element
        cart_ct(1) = crosses(jj,1) + r_vec(1)
        cart_ct(2) = crosses(jj,2) + r_vec(2)
        cart_ct(3) = 0.0
        cart_ct(4) = 1.0  
        IF ((SUM(ABS(MATMUL(simat2,cart_ct)))-1).LT.EPS) THEN
		  !print *, 'op rand van driehoek'
          first_cross = first_cross + 1
          rtemp = SQRT(SUM(crosses(jj,1:2)**2))
 
          IF (ABS(first_cross-1).LT.EPS) THEN ! Initialise the two points
            r1 = rtemp
            int_x(1) = jj
            r2 = rtemp
            int_x(2) = jj
          END IF 
            
          IF (rtemp.GT.r2) THEN
            r2 = rtemp
            int_x(2) = jj
          ELSE IF (rtemp.LT.r1) THEN
            r1 = rtemp
            int_x(1) = jj
          END IF
 
        END IF
      END SELECT
    END DO
 
    IF (outer_elem.EQ.elem) THEN ! Elements lie on each other
      IF (ABS(COS(theta)-crosses(int_x(1),1)/r1).LT.EPS) THEN ! r1 is not in the theta direction
        crosses(int_x(2),1:2) = crosses(int_x(1),1:2)
        r2 = r1
      END IF
      crosses(int_x(1),1:2) = 0.0
      r1 = 0.0
    END IF
       
    ! Case when r1 and r2 are VERY close: (will shortly divide by r2-r1)
    IF (ABS(r2-r1).LT.EPS) r2 = r2 + EPS

    ! Evaluate the vector basis functions at the appropriate locations (see MMB PhD 2002):
    cart_ct(3) = 0.0 ! z=0 plane
    cart_ct(1:2) = crosses(int_x(1),1:2) + r_vec(1:2) ! r1
    CALL POINT_EVALUATE_FACE_FUNCTIONS(elem,local_face,cart_ct(1:3),num_vbfs, &
		                               point1_term1,point1_term2)
    cart_ct(1:2) = 0.5*crosses(int_x(1),1:2) + 0.5*crosses(int_x(2),1:2) + r_vec(1:2) ! center of r1,r2
    CALL POINT_EVALUATE_FACE_FUNCTIONS(elem,local_face,cart_ct(1:3),num_vbfs, & 
                                       point2_term1,point2_term2)
    cart_ct(1:2) = crosses(int_x(2),1:2) + r_vec(1:2) ! r2
    CALL POINT_EVALUATE_FACE_FUNCTIONS(elem,local_face,cart_ct(1:3),num_vbfs, &
                                       point3_term1,point3_term2)

	! Calculate the coefficients m,c for (del dot ^z times VBF) (see MMB PhD 2002): 
    m_term1 = (point3_term1(1:num_vbfs) - point1_term1(1:num_vbfs))/(r2-r1)
	c_term1 = (point1_term1(1:num_vbfs)*r2 - point3_term1(1:num_vbfs)*r1)/(r2-r1)

	! Calculate the coefficients A,B,C for (^z times VBF) (see MMB PhD 2002):
    ! (First create the r-matrix and invert it.)
    rmat(1,1:3) = (/ r1**2           , r1          , 1.0 /)
    rmat(2,1:3) = (/ 0.25*(r1+r2)**2 , 0.5*(r1+r2) , 1.0 /)
    rmat(3,1:3) = (/ r2**2           , r2          , 1.0 /)
    CALL SGETRF ((3),(3),rmat,(3),ipiv_local,info)  
    IF(info.NE.0) CALL ERROR_FEMFEKO(1,4300,int1=info) ! Check for errors, shouldn't occur
    CALL SGETRI ((3),rmat,(3),ipiv_local,work_local,(3),info)
    IF(info.NE.0) CALL ERROR_FEMFEKO(1,4301,int1=info) ! Check for errors, shouldn't occur
    DO count1 = 1,num_vbfs
      DO count2 = 1,3
        fvec = (/ point1_term2(count1,count2) , point2_term2(count1,count2) , &
                  point3_term2(count1,count2) /)
        fvec = MATMUL(rmat,fvec)
        A_term2(count1,count2) = fvec(1)
        B_term2(count1,count2) = fvec(2)
        C_term2(count1,count2) = fvec(3)
      END DO
    END DO

    ! Integrate analytically in the r-direction:
    inner_term1 = inner_term1                                                    &
                  - (m_term1/(cj*k0))*( r2*EXP(-cj*k0*r2) - r1*EXP(-cj*k0*r1) +  &
                    (1/(cj*k0))*(EXP(-cj*k0*r2)-EXP(-cj*k0*r1)) )                &
                  - (c_term1/(cj*k0))*(EXP(-cj*k0*r2)-EXP(-cj*k0*r1))
    inner_term2 = inner_term2                                                           &
                  + (A_term2/(-cj*k0))*((r2**2)*EXP(-cj*k0*r2)-(r1**2)*EXP(-cj*k0*r1) - &
                    (2.0/(-cj*k0))*(r2*EXP(-cj*k0*r2) - r1*EXP(-cj*k0*r1)) +            &
                    (2.0/((-cj*k0)**2))*(EXP(-cj*k0*r2)-EXP(-cj*k0*r1)))                &
                  - (B_term2/(cj*k0))*( r2*EXP(-cj*k0*r2) - r1*EXP(-cj*k0*r1) +         &
                    (1/(cj*k0))*(EXP(-cj*k0*r2)-EXP(-cj*k0*r1)) )                       &
                  - (C_term2/(cj*k0))*(EXP(-cj*k0*r2)-EXP(-cj*k0*r1))
  END DO theta_loop
 
  ! Add the 1/(4*PI), Gre0 factor and add the delta_theta weighting:
  inner_term1 = dtheta/(4.0*PI*REAL(t_points)) * inner_term1
  inner_term2 = dtheta/(4.0*PI*REAL(t_points)) * inner_term2

  DEALLOCATE(ipiv_local) !
  DEALLOCATE(work_local) ! for <rmat> inversion

END SUBROUTINE INNER_INTEGRAL_SELF
!*******************************************************************************


SUBROUTINE EVALUATE_FACE_FUNCTIONS_QUAD(elem,local_face,num_qpts,num_vbfs,quad_xyz,quad_term1,quad_term2)
  USE basis_function, ONLY: POINT_EVALUATE_FACE_FUNCTIONS
  USE geometry
  USE nrtype
  USE problem_info
  USE CBAA_data
  USE quad_tables
  IMPLICIT NONE   
!*******************************************************************************
! Returns the following for element number <elem> on its local face <local_face>:
! * xyz's of quadrature points of the requested quadrature rule <num_qpts>on this
!   element's aperture face <quad_xyz>.
! * Del dot ^z cross (VBFs) of the first <num_vbfs> VBFs at the quadrature points,
!   <quad_term1> (weighted with face area and gauss weight).
! * ^z cross (VBFs) of the first <num_vbfs> VBFs at the quadrature points, 
!   <quad_term2> (weighted with face area and gauss weight).
!
! Note: the dimension specification is in F90 standard - see Metcalf 'Fortran 90
! explained', p.278 #2.
!
! MMB - 19 Nov 2000
! 2002-05-11: Generalised for trivial extension to other/higher order VBFs. MMB.
!*******************************************************************************
  INTEGER(I4B), INTENT(IN) :: elem,local_face
  INTEGER(I4B), INTENT(IN) :: num_qpts,num_vbfs
  REAL(SP), DIMENSION(num_qpts,3), INTENT(OUT) :: quad_xyz
  REAL(SP), DIMENSION(num_qpts,num_vbfs), INTENT(OUT) :: quad_term1
  REAL(SP), DIMENSION(num_qpts,num_vbfs,3), INTENT(OUT) :: quad_term2

  INTEGER(I4B) :: gcount
  REAL(SP), DIMENSION(4) :: simp,cart
  REAL(SP), DIMENSION(4,4) :: vmat
  INTEGER(I4B), DIMENSION(3) :: tempfacenodes
  REAL(SP) :: tri_face_area

  ! Calculate the vertices matrix:
  CALL SIMPLEX_COEFFICIENTS(elem,coord_mat=vmat)

  ! Calculate for simplex assignment:
  tempfacenodes = LOCAL_FACENODES(local_face)

  simp = 0.0 ! initialise

  ! Evaluate the basis functions at the Gaussian quadrature points:
  GAUSS_LOOP: DO gcount = 1,num_qpts

    ! Calculate the xyz position of the Gauss points by
    ! assigning the relavant simplex coordinates their values,
    ! the 4th one remains =0, because the point lies on the 
    ! element side (z=0 plane).
    simp(tempfacenodes) = quad_tri_rules(num_qpts)%rule(gcount,1:3)
    cart = MATMUL(vmat,simp)
    quad_xyz(gcount,1:3) = cart(1:3)

    ! Evaluate the basis functions at the current quadrature point:
    CALL POINT_EVALUATE_FACE_FUNCTIONS(elem,local_face,quad_xyz(gcount,1:3),num_vbfs,quad_term1(gcount,1:num_vbfs),&
	     quad_term2(gcount,1:num_vbfs,1:3))

    ! Multiply with the current quadrature weight:
	quad_term1(gcount,1:num_vbfs)     = quad_tri_rules(num_qpts)%rule(gcount,4) &
                                        * quad_term1(gcount,1:num_vbfs)
	quad_term2(gcount,1:num_vbfs,1:3) = quad_tri_rules(num_qpts)%rule(gcount,4) &
	                                    * quad_term2(gcount,1:num_vbfs,1:3)
  END DO GAUSS_LOOP
              
  ! Weight all the values with the face area:
  tri_face_area = FACE_AREA(elem,local_face)
  quad_term1    = tri_face_area * quad_term1
  quad_term2    = tri_face_area * quad_term2

END SUBROUTINE EVALUATE_FACE_FUNCTIONS_QUAD
!*******************************************************************************

END MODULE CBAA_SYS


