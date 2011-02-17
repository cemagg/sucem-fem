! Last changed: 
! 2005.03.23: Remove EA card support, error on sparse eigen-analysis,
! remove some XX_CARD_OLD cruft. NM
!  2003-03-26: New AB and TD cards added, BC card changed slightly.DBD.
! 2002-12-14: New iterative option added to FM card. DBD.
! 2002-05-09: Added <rtemp1> to SUBROUTINE FC_CARD. MMB.
! 2002-04-15: Added upgrade_type variable in SUBROUTINE AD_CARD. MMB.
! 2002-03-21: Added frac_upgrade variable in SUBROUTINE AD_CARD. MMB.
! 2002-03-17: Added AD card. MMB.
! 28 Feb 2002: Output info of element type updated. DBD.
! 21 Feb 2002: FC card changed. DBD. 
! 1 June 2001: 
! 27 Apr 2001: Pre-conditioner capability control from FM card added. DBD. 
! 16 Apr 2001: QMR possibility added  in FM card. MMB
! 9 Apr 2001: AX card added MMB

! 29 March 2001. gw_data added at line 19 and S_params at 518
! DEBUG_GW_PRE control flag added, line 179. 
! 20 March 2001. A9 card extended. Routine READ_LINE extended. 
! 16 March 2001. GW card extended. 
!  6 March 2001. Error corrected on line 770 DBD.


SUBROUTINE MESHIN_FEK_OLD(infilename)
   USE boundary_conditions
   USE coax_feed
   USE far_field_data
   USE frequency_data
   USE geometry
   USE gw_data
   USE inc_field
   USE matrix
   USE near_field_data
   USE nrtype
   USE output_error, ONLY: ERROR_FEMFEKO
   USE probe_feed
   USE problem_info
   USE quad_tables
   USE TD_source
   USE unit_numbers
   USE parseinput
   USE material_properties
   IMPLICIT NONE
!*******************************************************************************
! SUBROUTINE DESCRIPTION
!*******************************************************************************
! This subroutine handles mesh input from old fek files using the original 
! FEMFEKO card definitions (DBD/MMB/FJCM), prior to extensive changing during 2003 during
! full integration in FEKO. 
! The file requirements are detailed below. 
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! D B Davidson
!
!*******************************************************************************
! Last revised:
!*******************************************************************************
! 19 March 1997 by DBD. Tested on approx 1800 element data file; read correctly.
! Interface to FAM material and BC's not yet implemented. 
! Minor revision 8 April (selective unsorted node print added).
!
! Major revisions 28-29 April: 
!
! 1. Mesh input file now read twice to allow proper 
! dynamic memory allocation. All large arrays are now properly
! dynamically allocated without requiring any user knowledge
! of the mesh size.  
!
! 2. Input program control file introduced. Mesh file is now defined
! in this file. Debugging options are now controlled by this file.
! 
! Minor changes 16 May by DBD to debug control names.
! Changes 19 May by DBD to properly support reading and checking of 
! anisotropic materials (with diagonal constitutive tensors).
! Changes 20 May: element_edges re-dimensioned
! Minor addition 3 June: additional debugging variable added and prettied up.
!
! Minor addition 16 June: additional output control added.
! Ditto 23 June.
! 
! Minor addition 26 November 1998 by DBD. 
! Element TYPE also now read in.
!
! Change 24 January 1999: nodes numbered from 1-4, not 0-3. 
! Unused variables tidied up and NAMELISTs moved. 26 Jan 99 DBD.
! 
! Changed March 1999  : input from the feko file format, changed by RHG
! Corrected 21 Apr    : RHG
! Prettied up 27 Apr  : DBD
! Changed May 16 1999 : DBD. Problem with non-advancing READ encountered
!                       during port to SGI O2 under Irix 6.3
!                       List-directed IO used to simplify code.
! Changed May 18 1999 : DBD. Constitutive parameter definition 
!                       entirely via MA card.
! Changed May 27 1999:  Major code restructing to improve readability (use made
!                       of internal subprograms).
!                       Support added for FC and EI cards. Obsolete code
!                       w.r.t. physical size of meshed object (max_dim) removed.
! Changed June 9 1999:  Multiple FE cards checked for.
! Changed June 19&21 1999: New option added to average fields on element
!                       boundaries, using new FA card.
! Changed 31 Aug 1999:  DBD. On request from Ulrich Jakobus, 
!                       MP card changed to MA.
! Changed 17 Nov 1999:  DBD. Error in declarations of eps_r_matrix, mu_r_matrix
!                       corrected. CMPLX statements used in internal 
!                       subroutine MA_CARD (MIPS 7.2.1 compiler complained 
!                       about previous conversion method). 
! Changed late Jan 2000:DBD. Error handling streamlined: BO; BC cards added
!         and early Feb 
!
! Changed 15 Feb 2000:  DBD. A0 card added.
! Changed 17 Feb 2000:  DBD. FF and FR cards added.
! Changed 24 Feb 2000:  DBD. EI card extended.
! Changed 14 Mar 2000:  DBD. Data structures changed.
! Corrected 3 April 2000: DBD. port_counter initialized properly.
! Changed 26 Apr 2000:  DBD/RHG: 
! Changed  2 May 2000:  RHG.  Error checking for EA (ARPACK) card. 
! Changed  8 May 2000:  MMB.  Multiple FF,FE cards possible.   
! Changed 12 May 2000:  MMB/DBD: Minor corrections.
! Changed 17 May 2000:  DBD: for Compaq port.
! Changed 22 June 2000: DBD/MMB Merge.
! Changed 27 June 2000: DBD. FM card extended.
! Changed 12 July 2000: DBD. LT/QN support for GW analysis included. 
! Changed 20 July 2000. DBD. FM card extended further.                     
! Changed 20 Jan 2001.  MMB. Added analysis number system. 
! Changed 2 Feb 2001.   DBD. Added additional LT/QN element type in FC card
!                       and new FX card.
! Changed 9 Apr 2001.   MMB. Added the development-only, AX, coax feed card for CBAA.
! Changed 24 Mar 2003.  DBD. Added some support for A0 excitation in FE part of code. 
!                            Added preliminary Time Domain capabilities. 
! Changed  9 Dec 2003.  DBD. Routine name changed.
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data uses FEKO file format, extended to support Finite Element
! Analysis (FEA). This is done via cards, similar to NEC.
!
! The filename is xxxx.fek (this must be specified in INFILE)
! The BO card specifies at PEC ground plane, using for cavity-backed
! antenna analysis. (Existing FEKO card, with additional funcationality
! in FEKO, not used here). 
! The EI card provides additional control parameters for EIgenanalysis 
! (new card).
! The FA card specifies Field Averaging.
! The FC card specifies is the Finite element analysis Control card
! to be performed (new card).
! The FM card provides additional options for the Finite element Matrix 
! (new card). 
! The MA card describes the MAterial parameters for FEA (new card).
! The TE cards describe TEtrahedral elements (new card).
!
! Other cards and data are not presently used by femfeko.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Echoes input data to file 'outfilename' - does not do any sorting.
! or processing of input geometrical data. This is done in subsequent routines. 
!
!*******************************************************************************
   INTEGER(I4B) i_node,i, itemp, itemp1, itemp2, itemp3, itemp4, itemp5,&
                i_vert, ipos, i_mat,  i_elem, i_row, i_col, i_BC
                !  misc. temp. storage
   INTEGER(I4B) lowest_material, highest_material ! lowest/highest MA labels 
                                                  ! (last is same as MAX_MATERIALS)
   INTEGER(I4B) highest_material_referenced       ! highest TE label
   INTEGER(I4B) node_counter, element_counter, BC_counter, &
                port_counter ! counters
   REAL(SP), DIMENSION(3) :: xyztemp1, xyztemp2, xyztemp3, xyztemp4
   REAL(SP) :: rtemp1,rtemp2,rtemp3,rtemp4,rtemp5,rtemp6,rtemp7
   REAL(SP), DIMENSION(4,3) :: xyz_pos
   REAL(SP), DIMENSION(4) :: parms
   COMPLEX(SPC), DIMENSION(3,3) :: eps_r_matrix, mu_r_matrix ! Full matrices
  
   INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: temp_vertices, temp_nodes

   CHARACTER(FILENAMELENGTH), INTENT(IN) :: infilename

   CHARACTER(8) DATE
   CHARACTER(10) TIME
   CHARACTER(15) string, string1, filestring
   CHARACTER(5) stemp1, stemp2, stemp3, stemp4, stemp5
! Correction DBD 9 Dec 2003
   LOGICAL(LGT) new_node, EA_card_defined, PL_card_defined, QU_card_defined, TD_card_defined !,KU_card_defined
   INTEGER(I4B) NUM_VERTICES

   CALL GENERAL_DEFAULTS ! Set up a number of defaults.

   ! Write the output file header:
   ! PROGRAM INFORMATION
   WRITE(FILEOUT,'(20X,A)') 'PROGRAM INFORMATION'
   WRITE(FILEOUT,'(/,1X,A,A,/,1X,A,A,/,1X,A,/,1X,A,A,/,1X,A,A,/,1X,A,A)') &
                      'Program name:     ','FEMFEKO',          &
                      'Description:      ','3 Dimensional ElectroMagnetic', &
                      '                  Finite Element Method code', &
                      'Version number:   ',VERSION_NUMBER,     &
                      'Revision number:  ',REVISION_NUMBER,    &
                      'Last changed:     ',LAST_CHANGED 
   WRITE(FILEOUT,'(1X,A,A,3(/,1X,A))')                              &
                      'Authors:          ','D.B.Davidson, M.M.Botha & R.H.Geschke',              &
                      '                  (Dept E&E Engr., Univ. of Stellenbosch, South Africa)', &
                      '                  F.J.C.Meyer',                                           &
                      '                  (EMSS Pty Ltd, Stellenbosch, South Africa)'

   CALL DATE_AND_TIME(DATE,TIME)
   WRITE(FILEOUT,'(1X,A,A)') &
                      'Date of run(ymd): ',DATE 
   WRITE(FILEOUT,'(1X,A,A,A,A)') &
                      'Time run started: ',TIME(1:2),'h',TIME(3:4)
   ! ASSOCIATED FILES
   WRITE(FILEOUT,'(//,20X,A)') 'ASSOCIATED FILES'
   WRITE(FILEOUT,'(/,1X,A,A)') 'Control data read from file: ', infilename
   WRITE(FILEOUT,'(1X,A,A)')   'Mesh read from file:         ', meshfilename
   WRITE(FILEOUT,'(1X,A,A)')   'Output data written to file: ', outfilename

   CALL MESH_SIZE ! Does a dummy read to find number of elements, materials etc.

   CALL DATA_ALLOCATE ! Dynamically allocates memory for geometry variables.

   EA_card_defined = .FALSE. 
!*******************************************************************************
!**** Data input section ****
!*******************************************************************************
   ! Now re-open the file and store the data 
        
   OPEN (UNIT=MESHIN,STATUS='OLD',file=meshfilename)

   !(Re-)initialisation of counters****************************
   node_counter = 1
   element_counter = 0

   lowest_material  = 0            ! Keeps track of material numbers
   highest_material = 0 
   highest_material_referenced = 0 

   BC_counter  = 0
   port_counter= 0

   FE_counter = 0          ! These cards are effectively counted again; these
   FF_counter = 0          ! are used as indices for storage.
   A0_counter = 0
   A8_counter = 0
   A9_counter = 0
   AB_counter = 0
   AX_counter = 0
   DP_counter = 0
   MA_counter = 0

   !******************************************************   

   !Read header of FEK file.
   HEADER_LOOP2: DO
     READ(MESHIN,'(A9)') string
     IF (string.EQ.'BEGIN_FEM') THEN
       EXIT HEADER_LOOP2
     END IF
   END DO HEADER_LOOP2

   CARD_LOOP2: DO 
     ! Read initially without advancing to next record,
     ! (or equivalently backspacing to start of entire record); 
     ! this enables the type of card to be determined.
       READ(MESHIN,'(A5)',ADVANCE='NO') string 
     CARD_TEST: IF (string.EQ.'A0') THEN
       CALL A0_CARD ! Incident plane wave
     ELSE IF (string.EQ.'A8') THEN
       CALL A8_CARD ! CBAA current probe excitation
     ELSE IF (string.EQ.'A9') THEN
       CALL A9_CARD ! CBAA current probe excitation
! Added DBD 25 March 2003
     ELSE IF (string.EQ.'AB') THEN        
       CALL AB_CARD ! Absording boundary condition card
! End added DBD 25 March 2003
     ELSE IF (string.EQ.'AD') THEN  
       CALL AD_CARD ! ADaptive FEM control card
     ELSE IF (string.EQ.'AX') THEN
       CALL AX_CARD ! CBAA coax excitation (experimental)
      ELSE IF (string.EQ.'BC') THEN
       CALL BC_CARD ! Boundary condition
     ELSE IF (string.EQ.'BO') THEN  
       CALL BO_CARD ! Ground plane
! Added DBD 28 March 2003
     ELSE IF (string.EQ.'DP') THEN        
       CALL DP_CARD ! Define point card. 
! End added DBD 28 March 2003
     ELSE IF (string.EQ.'EA') THEN  
        CALL EA_CARD ! Eigenanalysis ARPACK control card.
     ELSE IF (string.EQ.'EI') THEN  
        PRINT*, "EI_card now handled in input.dat, please remove from .fek file."
        STOP
     ELSE IF (string.EQ.'FA') THEN  
       CALL FA_CARD ! Field averaging command
     ELSE IF (string.EQ.'FC') THEN  
        PRINT*, "FC card now handled in input.dat, please remove from .fek file."
        STOP
     ELSE IF (string.EQ.'FE') THEN  
       CALL FE_CARD ! Near field card. 
     ELSE IF (string.EQ.'FF') THEN  
       CALL FF_CARD ! FEKO Far field card. 
     ELSE IF (string.EQ.'FM') THEN  
!!$        CALL FM_CARD ! Finite element Matrix card
        PRINT*, "FM card now handled in input.dat, please remove from .fek file."
        STOP
     ELSE IF (string.EQ.'FR') THEN  
        PRINT*, "FR card now handled in input.dat, please remove from .fek file."
        STOP
     ELSE IF (string.EQ.'FX') THEN  
        PRINT*, "FX card now handled in input.dat, please remove from .fek file."
        STOP
     ELSE IF (string.EQ.'GW') THEN  
       CALL GW_CARD ! Guided wave options card
     ELSE IF (string.EQ.'MA') THEN  
        PRINT*, "MA card now handled in aux input file, please remove from .fek file."
        STOP
! Added DBD 25 March 2003 and 05 June 2003
     ELSE IF (string.EQ.'PL') THEN        
       CALL PL_CARD ! Perfectly matched layer card
     ELSE IF (string.EQ.'QU') THEN        
       CALL QU_CARD ! Cuboid card
     ELSE IF (string.EQ.'TD') THEN        
       CALL TD_CARD ! Time Domain analysis card
! End added DBD 25 March 200305 June 2003 
     ELSE IF (string.EQ.'TE') THEN        
       CALL TE_CARD ! Tetrahedral element card
     ELSE IF (string.EQ.'END_F') THEN ! End of the file card.
       EXIT CARD_LOOP2
     ELSE ! Comment cards or presently unimplemented card.
       READ (MESHIN,'(A15)') string 
     END IF CARD_TEST
   END DO CARD_LOOP2


   NUM_analysis = analysis_counter  ! only now the number of analysis are known
   num_nodes = node_counter-1       ! now the exact number of nodes is known 
    
   !************************************************************************
   !WRITE(FILEOUT,*) 'This shows the arrays vertices and nodes***************'
   !WRITE(FILEOUT,*) 'nodes**************************************************' 
   !  DO i = 1,num_elements
   !   WRITE(FILEOUT,*) elements(i)%nodes
   !  END DO
   !  
   !WRITE(FILEOUT,*) 'vertices***********************************************' 
   !  DO i = 1,num_nodes
   !     WRITE(FILEOUT,*) vertices(i)%coord
   !  END DO
   !WRITE(FILEOUT,*) 'the number of elements: ', element_counter 
   !***********************************************************************
   
   CLOSE (UNIT=MESHIN)
     

   CALL ERROR_CHECKING


   RETURN

CONTAINS ! Internal subprograms.

!*******************************************************************************

SUBROUTINE DATA_ALLOCATE

!*******************************************************************************
!**** Data allocation section **** 
!*******************************************************************************
! At this stage only the number of elements are known.
! For the number of nodes, edges and faces, worst case estimations are made.
!*******************************************************************************
  MAX_NODES =    num_elements*4! Worst case estimation of maximum number of nodes
  MAX_ELEMENTS = num_elements  ! Exact number of elements
  MAX_EDGES = MAX_ELEMENTS*6   ! Maximum number of edges & faces.
  MAX_FACES = MAX_ELEMENTS*4   ! Note: worst case - 
                               ! 6 edges and 4 faces per tet.               
  MAX_BOUNDARIES = BC_counter
  NUM_A0_cards = A0_counter
  NUM_A8_cards = A8_counter
  NUM_A9_cards = A9_counter
  NUM_AX_cards = AX_counter
  NUM_FE_cards = FE_counter
  NUM_FF_cards = FF_counter
  NUM_MA_cards = MA_counter
  MAX_PORTS = port_counter

  ! See module DATASTRC for definitions. 

  ! Vertices, elements, egdes, faces:  
  PRINT*, MAX_ELEMENTS, MAX_FACES, MAX_EDGES, MAX_NODES
  IF (max_elements > 0) THEN
     CALL ALLOCATE_GEOMETRY(alloc_vertices=MAX_NODES, alloc_edges=MAX_EDGES, &
          alloc_faces=MAX_FACES, alloc_elements=MAX_ELEMENTS)
  END IF
  
  ! Material parameters.
  ALLOCATE (A0data(NUM_A0_cards+1))   ! A0 card parameters
  ALLOCATE (A8data(NUM_A8_cards+1))   ! A8 card parameters
  ALLOCATE (A9data(NUM_A9_cards+1))   ! A8 card parameters
  ALLOCATE (AXdata(NUM_AX_cards+1))   ! AX card parameters
  ALLOCATE (FEdata(NUM_FE_cards+1))   ! FE card parameters
  ALLOCATE (FFdata(NUM_FF_cards+1))   ! FF card parameters
  ALLOCATE (MAdata(NUM_MA_cards+1))   ! FF card parameters

  ! Boundary conditions and ports.
  ALLOCATE (BCs(MAX_BOUNDARIES))
  ALLOCATE (ports(MAX_PORTS))

! Start DBD extension 29 March
  ! S parameter matrix.
  ALLOCATE (S_params(MAX_PORTS,MAX_PORTS))
! End DBD extension 29 March

! DBD addition 28 March 2003
  ALLOCATE(DPoints(num_DPoints))
  ALLOCATE(ABCs(num_ABCs))
! End DBD addition 28 March 2003


END SUBROUTINE DATA_ALLOCATE

!*******************************************************************************


SUBROUTINE GENERAL_DEFAULTS
!*******************************************************************************
! Input program control data & defaults 
!*******************************************************************************

  ! Set up defaults for namelists
  !   NAMELIST/PROBLEM/ 
  ! although they don't belong to this namelist? Perhaps in an earlier life?
  FIELD_AVERAGING = .false.
  eps_field_averaging = EPS 


  ! Set up defaults for other variables read in from .fek file.
  NUM_VERTICES = 4  !  This is for the 3D case

  ! Defaults for (triangular) quadrature
  ! gauss_points is defined in module quad_tables
  ! NM: These are initialised in input_problem_info now
!!$  SELECT CASE (HIERARCHAL_ORDER)
!!$  CASE(1)
!!$    gauss_points = 3
!!$  CASE(2)
!!$    gauss_points = 4
!!$  END SELECT

  ! TE10 default
  num_wg_modes_m = 1
  num_wg_modes_n = 0
  
  QU_card_defined = .false.
  TD_card_defined = .false.
  
  CALL geometry_init

END SUBROUTINE GENERAL_DEFAULTS

!*******************************************************************************

SUBROUTINE ERROR_CHECKING
  USE scattering_analysis_data
!*******************************************************************************
! Do some error checking on the data, especially material parameters.
!*******************************************************************************
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: port_check
  INTEGER(I4B) ABC_num
!!$  IF (highest_material_referenced.GT.highest_material) THEN
!!$    CALL ERROR_FEMFEKO(1,4027,int1=highest_material_referenced)
!!$  END IF     

  IF (SPARSE .AND. REAL_EIGEN_ANALYSIS) THEN
    IF (.NOT.(EA_card_defined)) THEN
      CALL ERROR_FEMFEKO(1,4028)
    END IF
  END IF

  IF (SPARSE .AND. CMPLX_EIGEN_ANALYSIS) THEN
      CALL ERROR_FEMFEKO(1,4029) ! This will soon be possible
  END IF 
 
  IF ((CMPLX_EIGEN_ANALYSIS.OR.REAL_EIGEN_ANALYSIS).AND.(NUM_FF_cards.GT.0)) THEN
    CALL ERROR_FEMFEKO(1,4044)
  END IF

  IF ((CMPLX_EIGEN_ANALYSIS.OR.REAL_EIGEN_ANALYSIS).AND.(NUM_A8_cards.GT.0)) THEN
    CALL ERROR_FEMFEKO(1,4045)
  END IF

  IF ((CMPLX_EIGEN_ANALYSIS.OR.REAL_EIGEN_ANALYSIS).AND.(NUM_A0_cards.GT.0)) THEN
    CALL ERROR_FEMFEKO(1,4046)
  END IF

  IF (CBAA_ANALYSIS) THEN
    ! Check that nodes(z) <= 0 in CBAA case:
    ! Check if a ground plane is defined
    IF(.NOT.PEC_ground) CALL ERROR_FEMFEKO(1,4031)
  END IF

! Changed DBD 16 March 01

  IF (GW_ANALYSIS) THEN
  ! Test present geometrical assumptions: (some are now tested in WG_ORIENTATION)
    IF (num_ports.LT.1) THEN 
      CALL ERROR_FEMFEKO(1,4036)
    END IF
    IF (num_A9_cards.NE.num_ports) THEN
      CALL ERROR_FEMFEKO(1,4050)
    END IF

  END IF

! End DBD changes 16 March 01

! Comment changed DBD 1 June 01
  ! Check that ports are correctly labeled.
  IF (GW_ANALYSIS) THEN
    ALLOCATE (port_check(NUM_ports))
    port_check = 0
! Next line removed DBD 1 June 01. Code handles multiple ports. 
!    IF (NUM_ports.NE.2) CALL ERROR_FEMFEKO(1,4050)
    DO i_BC = 1,NUM_BCs
      IF (BCs(i_BC)%type.EQ.1) THEN
! Next line changed DBD 1 June 01: port_num.LT.1  (was LT.0). 
        IF ((BCs(i_BC)%port_num.GT.NUM_ports).OR.(BCs(i_BC)%port_num.LT.1)) THEN
          CALL ERROR_FEMFEKO(1,4051,int1=i_BC,int2=BCs(i_BC)%port_num)
        END IF
        port_check(BCs(i_BC)%port_num) = 1
      END IF
    END DO
    IF (SUM(port_check).NE.NUM_ports) CALL ERROR_FEMFEKO(1,4052)
    DEALLOCATE (port_check)
  END IF

  ! Check that illegal cards are not present in the GW case:
  IF (GW_ANALYSIS) THEN
    IF (NUM_A0_cards.GT.0) CALL ERROR_FEMFEKO(1,4055)
    IF (NUM_A8_cards.GT.0) CALL ERROR_FEMFEKO(1,4056)
    IF (NUM_FF_cards.GT.0) CALL ERROR_FEMFEKO(1,4057)
  END IF

! Start DBD additions 19 Feb and 25 Apr 2001
  IF(ELEMENT_TYPE.EQ.2.AND.SCALE_BY_EDGE_LENGTH) THEN
    CALL ERROR_FEMFEKO(1,4049)
  END IF 
! DBD change 24 Feb 2002.
  IF(ELEMENT_TYPE.GT.1.AND.HIERARCHAL_ORDER.GT.1.AND..NOT.CUBATURE) THEN
    CALL ERROR_FEMFEKO(1,4050)
  END IF
! End 24 Feb 02 change.
 
  IF(USE_PRE_CONDITIONER.AND.(SOLVER_TYPE.EQ.2.OR.SOLVER_TYPE.EQ.3)) THEN
    CALL ERROR_FEMFEKO(1,4082)
  END IF
! End DBD additions 19 Feb and 25 Apr 2001


! Added DBD 28 March and 8 May 2003
  IF(TD_ANALYSIS.AND..NOT.QU_card_defined) THEN
    CALL ERROR_FEMFEKO(1,4089)
  END IF
!  IF(TD_ANALYSIS.AND.SCAT_FIELD.AND..NOT.KU_card_defined) THEN
!    CALL ERROR_FEMFEKO(1,4097)
!  END IF
  IF(TD_ANALYSIS) THEN
    IF (num_ABCs.NE.6) THEN
      CALL ERROR_FEMFEKO(1,4089)
    END IF 
    DO ABC_num = 1,num_ABCs
      IF (ABCs(ABC_num)%label.LT.1.OR.ABCs(ABC_num)%label.GT.6) THEN
        CALL ERROR_FEMFEKO(1,4089)
	  END IF
      IF (ABCs(ABC_num)%label.NE.ABC_num) THEN
	    CALL ERROR_FEMFEKO(1,4089)
	  END IF
    END DO
	IF (PML_PRESENT.AND..NOT.SCAT_FIELD) THEN
	  CALL ERROR_FEMFEKO(1,4099)
	END IF
  END IF 
! End added DBD 28 March and 8 May 2003


END SUBROUTINE ERROR_CHECKING

!*******************************************************************************


SUBROUTINE MESH_SIZE
!*******************************************************************************
! Determine size of mesh, number of materials, boundaries etc.
!*******************************************************************************

   element_counter = 0  
   MA_counter = 0
   A0_counter = 0
   A8_counter = 0
   A9_counter = 0
! Added DBD 26 March 2003
   AB_counter = 0
! End added DBD 26 March 2003
! Added MMB 9 Apr
   AX_counter = 0
! end MMB
   BC_counter = 0
! Added DBD 28 March 2003
   DP_counter = 0
! End added DBD 28 March 2003
   FE_counter = 0
   FF_counter = 0
   port_counter = 0
   ! The following is a "dummy" read to establish data sizes. 
   OPEN (UNIT=MESHIN,STATUS='OLD',file=meshfilename)
   ! Read and discard until the first TE card is encountered or  
   ! the "end of FEM" block marker.
   CARD_LOOP1: DO 
       READ (MESHIN,'(A5)',ADVANCE='NO') string 
     IF (string.EQ.'END_F') THEN ! End of FEM block
       EXIT CARD_LOOP1
     ELSE IF(string.EQ.'TE') THEN 
       CALL READ_LINE(string,itemp,itemp1,itemp2,itemp3,itemp4, & 
     xyztemp1(1),xyztemp1(2),xyztemp1(3))
       READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp2
       READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp3
       READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp4
       element_counter = element_counter+1
     ELSE IF(string.EQ.'MA') THEN
       MA_counter = MA_counter+1
       CALL READ_LINE(string,itemp) ! Read material label 
       IF (itemp.GT.MAX_MATERIALS) THEN
         MAX_MATERIALS = itemp
       END IF
     ELSE IF(string.EQ.'BC') THEN
       BC_counter = BC_counter+1
! Error corrected DBD: itemp1 -> itemp
       CALL READ_LINE(string,itemp) ! Read boundary condition card
       IF (itemp.EQ.1) THEN
         port_counter = port_counter + 1
       END IF
     ELSE IF(string.EQ.'A0') THEN
       A0_counter = A0_counter+1
       CALL READ_LINE(string) ! Read incident plane wave card
     ELSE IF(string.EQ.'A8') THEN
       A8_counter = A8_counter+1
       CALL READ_LINE(string) ! Read CBAA current probe card
! Added MMB 9 Apr
     ELSE IF(string.EQ.'AX') THEN
       AX_counter = AX_counter+1
       CALL READ_LINE(string) ! Read coax excitation card
! end MMB
     ELSE IF(string.EQ.'FE') THEN
       FE_counter = FE_counter+1
       CALL READ_LINE(string) ! Read near field card
     ELSE IF(string.EQ.'FF') THEN
       FF_counter = FF_counter+1
       CALL READ_LINE(string) ! Read far field card
     ELSE IF(string.EQ.'A9') THEN
       A9_counter = A9_counter+1
       CALL READ_LINE(string) ! Read port card
     ELSE IF(string.EQ.'AB') THEN
       AB_counter = AB_counter+1
       CALL READ_LINE(string) ! Read absorbing boundary condition card
     ELSE IF(string.EQ.'DP') THEN
       DP_counter = DP_counter+1
       CALL READ_LINE(string) ! Read define point card.
     ELSE IF(string.EQ.'PL') THEN
       PML_PRESENT = .TRUE.
       CALL READ_LINE(string) ! Read define point card.
     ELSE
       READ(MESHIN,'(A5)') string ! Otherwise read card.
     END IF
   END DO CARD_LOOP1
   num_elements = element_counter 
   num_BCs   = BC_counter
! Added DBD 26 March 2003
   num_ABCs  = AB_counter
! End added DBD 26 March 2003
   num_ports = port_counter
   num_DPoints = DP_counter
   
   CLOSE (UNIT=MESHIN)
END SUBROUTINE MESH_SIZE

!*******************************************************************************


SUBROUTINE A0_CARD
!***************************************************
! The A0 card: FEKO incident plane wave card. Optional. 
!***************************************************
! Format I1 (I2) I3 I4 (I5) R1 R2 R3 R4 R5 R6 R7
!***************************************************
! I1 = 0  : New excitation.
!    = 1  : Additional excitation - not yet implemented in FEMFEKO.
! I2      : Not used, but should have zero "place-holders".
! I3      : Number of theta angles. 
! I4      : Number of phi angles.
! I5      : Not used, but should have zero "place-holders".
! R1      : Magnitude of Einc (V/m).
! R2      : Phase of Einc  (degrees)
! R3      : Theta angle of incidence (degrees).
! R4      : Phi angle of incidence (degrees).
! R5      : Eta angle of polarization (degrees).
! R6      : Delta theta (degrees).
! R7      : Delta phi (degrees).
!
!***************************************************

  A0_counter = A0_counter+1
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5, & 
    rtemp1,rtemp2,rtemp3,rtemp4,rtemp5,rtemp6,rtemp7)
  SELECT CASE(itemp1)
  CASE(0:1)
    A0data(A0_counter)%anfl = itemp1
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4000,int1=itemp1,int2=A0_counter)
  END SELECT
  A0data(A0_counter)%NThei = itemp3
  A0data(A0_counter)%NPhii = itemp4
  A0data(A0_counter)%EiR1  = rtemp1
  A0data(A0_counter)%EiR2  = rtemp2
  A0data(A0_counter)%EiR3  = rtemp3
  A0data(A0_counter)%EiR4  = rtemp4
  A0data(A0_counter)%EiR5  = rtemp5
  A0data(A0_counter)%DThei = rtemp6
  A0data(A0_counter)%DPhii = rtemp7

  ! Force this to be a new excitation if it is the first one read:
  IF (analysis_counter.EQ.0) THEN
    IF (A0data(A0_counter)%anfl.NE.0) THEN
      A0data(A0_counter)%anfl = 0
      CALL ERROR_FEMFEKO(0,4001,int1=A0_counter)
    END IF
  END IF

  ! Enforce dtheta & dphi = 0 if there is only one observation point, this simplifies later processing
  IF (A0data(A0_counter)%NThei*A0data(A0_counter)%NPhii.EQ.1) THEN
    A0data(A0_counter)%DThei = 0.0 ! Must add a warning that this is enforced
    A0data(A0_counter)%DPhii = 0.0 ! Must add a warning that this is enforced
  END IF

  ! An A0 1 card with multiple incident directions is only acceptable if all
  ! other A0 cards already belonging to this analysis, only have one incident direction,
  ! otherwise this A0 is treated as a new source:
  IF (A0data(A0_counter)%NThei*A0data(A0_counter)%NPhii.GT.1) THEN
    DO itemp1 = 1,A0_counter-1
      IF ((analysis_counter.EQ.A0data(itemp1)%analysis_no).AND. &
        (A0data(itemp1)%NThei*A0data(itemp1)%NPhii.GT.1)) THEN ! an A0 card with multiple directions
                                                               ! already exists in the current analysis
        IF (A0data(A0_counter)%anfl.NE.0) THEN
          A0data(A0_counter)%anfl = 0
          CALL ERROR_FEMFEKO(0,4002,int1=A0_counter)
         END IF

      END IF
    END DO
  END IF

  ! Assign the analysis number:
  SELECT CASE (A0data(A0_counter)%anfl)
  CASE (0) ! Start of new analysis
    analysis_counter = analysis_counter + 1
  CASE (1) ! Part of current analysis
  END SELECT
  A0data(A0_counter)%analysis_no = analysis_counter

END SUBROUTINE A0_CARD

!*******************************************************************************

SUBROUTINE A8_CARD
!***************************************************
! The A8 card specifies the position, dimensions and 
! value of an ideal current probe within the CBAA
! FEM volume.
!***************************************************
! Format I1 I2 I3 I4 I5 R1 R2 R3 R4 R5 R6 R7
!***************************************************
! I1 = 0,1 : New/Additional excitation.
! I2,I3  : Not presently used, but should have zero "place-holders"
! I4 = 1,2,3: Direction of current probe (+x,+y or +z direction respectively)
! I5     : Not presently used, but should have zero "place-holders"
! R1     : Magnitude of current on probe.
! R2     : Phase of current on probe.
! R3,R4,R5 : (x,y,z) position of origin of current probe (m). 
! R6     : Probe radius (m).
! R7     : Probe length and sense (+ or -) 
!***************************************************

  A8_counter = A8_counter+1
  CALL READ_LINE(string,A8data(A8_counter)%anfl,itemp2,itemp3,A8data(A8_counter)%prdir,  &
         itemp5,A8data(A8_counter)%prabs,A8data(A8_counter)%prphase,                     &
         A8data(A8_counter)%prx0,A8data(A8_counter)%pry0,                                &
         A8data(A8_counter)%prz0,A8data(A8_counter)%prrad,A8data(A8_counter)%prlen)

  SELECT CASE(A8data(A8_counter)%anfl)
  CASE(0) ! Valid data
    CONTINUE
  CASE(1) ! Enforce new excitation if this is the first analysis
    IF (analysis_counter.EQ.0) THEN
      A8data(A8_counter)%anfl = 0
      CALL ERROR_FEMFEKO(1,4005,int1=A8_counter)
    END IF
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4003,int1=A8data(A8_counter)%anfl,int2=A8_counter)
  END SELECT

  SELECT CASE(A8data(A8_counter)%prdir)
  CASE(1:3)
    CONTINUE ! Valid data
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4004,int1=A8data(A8_counter)%prdir,int2=A8_counter)
  END SELECT

  ! Assign the analysis number:
  SELECT CASE (A8data(A8_counter)%anfl)
  CASE (0) ! Start of new analysis
    analysis_counter = analysis_counter + 1
  CASE (1) ! Part of current analysis
  END SELECT
  A8data(A8_counter)%analysis_no = analysis_counter

END SUBROUTINE A8_CARD


!*******************************************************************************


! Extended DBD 20 Mar 01 and 2 Apr 01
SUBROUTINE A9_CARD
!***************************************************
! The A9 card specifies the relative excitation at 
! a specified port, the outward directed normal
! to the port and the "sense" of the port via a 
! tangential vector in the direction of the positive
! sense of the mode(s). 
!***************************************************
! Format I1 I2 I3 (I4) (I5) R1 R2 R3 R4 R5 R6 R7 R8
!***************************************************
! I1 = 0 : New excitation.
!    = 1 : Additional excitation.
! I2     : Port label that this card relates to.
! I3 = 1 : Port type=TE (mode defined on GW card)
! I4,I5  : Not presently used, but should have zero "place-holders"
! R1,R2  : Relative port excitation (Re,Im).
! DBD additions:
! R3, R4, R5
!        : Outward-directed normal at port. 
!        : x, y and z components. (Note: need not be unit normal). 
! R6, R7, R8
!        : Vector "sense" of the port (to resolve possible phase
!        : ambiguities). 
!        : x, y and z components of tangential vector to define 
!        : direction of postitive mode. 
! End DBD additions
!             
!***************************************************

  A9_counter = A9_counter+1
! DBD addition - normal and tangential added.
  CALL READ_LINE(string,A9data(A9_counter)%anfl,A9data(A9_counter)%port_num,        &
                 A9data(A9_counter)%port_type,itemp4,itemp5,                        &
                 A9data(A9_counter)%re_excitation,A9data(A9_counter)%im_excitation, &
                 A9data(A9_counter)%normal(1),A9data(A9_counter)%normal(2),         &
                 A9data(A9_counter)%normal(3),                                      &
                 A9data(A9_counter)%tangent(1),A9data(A9_counter)%tangent(2),       &
                 A9data(A9_counter)%tangent(3) ) 

  SELECT CASE(A9data(A9_counter)%anfl)
  CASE(0) ! Valid data
    CONTINUE
  CASE(1) ! Enforce new excitation if this is the first analysis
    IF (analysis_counter.EQ.0) THEN
      A9data(A9_counter)%anfl = 0
      CALL ERROR_FEMFEKO(0,4049,int1=A9_counter)
    END IF
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4048,int1=A9data(A9_counter)%anfl,int2=A9_counter)
  END SELECT

  ! Check that the port number is valid:
  IF ((A9data(A9_counter)%port_num.LT.1).OR.(A9data(A9_counter)%port_num.GT.NUM_ports)) THEN
    CALL ERROR_FEMFEKO(1,4053,int1=A9_counter,int2=A9data(A9_counter)%port_num)
  END IF

  ! Check that this analysis does not already contain an A9 card for this port:
  IF (A9data(A9_counter)%anfl.EQ.1) THEN
    DO itemp = 1,A9_counter-1
      IF ((A9data(itemp)%analysis_no.EQ.analysis_counter).AND. &
        (A9data(itemp)%port_num.EQ.A9data(A9_counter)%port_num)) THEN
        CALL ERROR_FEMFEKO(1,4054,int1=A9data(A9_counter)%port_num,int2=A9_counter)
      END IF  
    END DO
  END IF

  ! Assign the analysis number:
  SELECT CASE (A9data(A9_counter)%anfl)
  CASE (0) ! Start of new analysis
    analysis_counter = analysis_counter + 1
  CASE (1) ! Part of current analysis
  END SELECT
  A9data(A9_counter)%analysis_no = analysis_counter

END SUBROUTINE A9_CARD
!*******************************************************************************
! End extensions DBD 20 Mar 01


SUBROUTINE AB_CARD
!***************************************************
! The AB card: FEMFEKO Absorbing Boundary condition card. 
! Required to define an ABC. Parameters are similar to 
! A9 card. Note that normal must be unitary!
!***************************************************
! Format I1 I2 I3 (I4) (I5) R1 (R2) R3 R4 R5
!***************************************************
! I1      : unused
! I2 =n   : Integer n is a unique label for each ABC.
! I3 =1   : ABC type: 1st order. 
! I4,I5   : Not presently used, but should have zero "place-holders"
! R1      : Y_c (intrinsic impedance of infinite medium)
! R2      : unused
! R3, R4, R5
!        : Outward-directed unit normal at ABC. 
!        : x, y and z components.
!
!***************************************************
  INTEGER(I4B) BC_num, icount 
  LOGICAL(LGT) BC_card_found 

  AB_counter = AB_counter+1
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5,rtemp1,rtemp2,rtemp3,rtemp4,rtemp5) 

  ! Check for unit normal
  IF(ABS((rtemp3**2+rtemp4**2+rtemp5**2)-1.0_SP).GE.EPS) THEN
    CALL ERROR_FEMFEKO(1,4085,int1=itemp2)
  END IF

  ! Find the relevant BC card, or warn that is has not yet been defined or has already been defined.
  BC_card_found = .false.
  DO icount = 1,BC_counter
    IF (itemp2.EQ.BCs(icount)%ABC_num) THEN 
	  IF (BC_card_found) THEN
        CALL ERROR_FEMFEKO(1,4086,int1=itemp2)
	  END IF
	  BC_card_found = .true.
	  BC_num = icount
	END IF
  END DO 
  IF(.NOT.BC_card_found) THEN 
    CALL ERROR_FEMFEKO(1,4087,int1=itemp2)
  END IF  

  ABCs(AB_counter)%BC_num    = BC_num
  ABCs(AB_counter)%label     = itemp2
  ABCs(AB_counter)%type      = itemp3
  ABCs(AB_counter)%Yc        = rtemp1
  ABCs(AB_counter)%normal(1) = rtemp3
  ABCs(AB_counter)%normal(2) = rtemp4
  ABCs(AB_counter)%normal(3) = rtemp5
  
  IF(ABCs(AB_counter)%type.NE.1) THEN
    CALL ERROR_FEMFEKO(1,4095,int1=ABCs(AB_counter)%label)
  END IF        


END SUBROUTINE AB_CARD

!*******************************************************************************




SUBROUTINE AD_CARD
  USE adaptive_fem
  IMPLICIT NONE
!***************************************************
! The presence of the AD card indicates that
! adaptive FEM related actions must be performed,
! as specified by the card parameters.
!***************************************************
! Format I1 I2 I3 I4 I5 R1 R2 R3 R4 R5 R6 R7
!***************************************************
! I1 = 0: No action.
!    = 1: Write the element orders to a file.
!    = 2: Read the element orders from a file.
!    = 3: Read & Write the element orders from a file.
! I2 = 0: No error analysis performed.
!    = 1: Explicit residual error analysis.
!    = 2: Random selection of elements.
!    = 3: Element residual method. 0.5:0.5 weighing, L^2 error norm.
! I3 = 0: No upgrades. (Upgrade to CT/LN.)
!    = 1: Upgrade to LT/LN.
!    = 2: Upgrade to LT/QN.
!    = 3: Upgrade to QT/QN.
!    = 4: Graded upgrades, even distribution of LT/LN and 
!       : all higher order elements.
! I4    : Place holder.
! I5    : Place holder.
! R1    : R1 \in [0,1] Relative weighting of element to face contributions
!         to the residual estimator (R1*element + (1-R1)*face).
! R2    : R2 \in [0,1] Fraction of elements of which order must be upgraded.
!***************************************************

  ! Set flag that this is a adaptive run:
  ADAPTIVE_ANALYSIS = .TRUE.

  CALL READ_LINE(string, ADdata%file_orders, ADdata%analysis_type, &
                 ADdata%upgrade_type, itemp4, itemp5, ADdata%relative_res, &
				 ADdata%frac_upgrade)

  ! Check that <ADdata%file_orders> is in range:
  SELECT CASE (ADdata%file_orders)
  CASE (0:2)
    CONTINUE
  CASE DEFAULT
    STOP 'IE: Out-of-range <ADdata%file_orders> in SUBROUTINE AD_CARD.'
  END SELECT

  ! Check that <ADdata%analysis_type> is in range:
  SELECT CASE (ADdata%analysis_type)
  CASE (0:3)
    CONTINUE
  CASE DEFAULT
    STOP 'IE: Out-of-range <ADdata%analysis_type> in SUBROUTINE AD_CARD.'
  END SELECT

  ! Check that <ADdata%upgrade_type> is in range:
  SELECT CASE (ADdata%upgrade_type)
  CASE (0:4)
    CONTINUE
  CASE DEFAULT
    STOP 'IE: Out-of-range <ADdata%upgrade_type> in SUBROUTINE AD_CARD.'
  END SELECT

  ! Check that <ADdata%relative_res> is in range:
  IF ((ADdata%relative_res.LT.0.0).OR.(ADdata%relative_res.GT.1.0)) THEN
    STOP 'IE: Out-of-range <ADdata%relative_res> in SUBROUTINE AD_CARD.'
  END IF

  ! Check that <ADdata%frac_upgrade> is in range:
  IF ((ADdata%frac_upgrade.LT.0.0).OR.(ADdata%frac_upgrade.GT.1.0)) THEN
    STOP 'IE: Out-of-range <ADdata%frac_upgrade> in SUBROUTINE AD_CARD.'
  END IF

END SUBROUTINE AD_CARD
!*******************************************************************************


SUBROUTINE AX_CARD
!***************************************************
! The AX card specifies the position, dimensions and 
! value of an coax feef aperture on the PEC, cavity 
! boundary within the CBAA FEM volume.
!***************************************************
! Format I1 I2 I3 I4 I5 R1 R2 R3 R4 R5 R6 R7
!***************************************************
! I1 = 0,1 : New/Additional excitation.
! I2,I3  : Not presently used, but should have zero "place-holders"
! I4 = 1,2,3: Direction of current probe (+x,+y or +z direction respectively)
! I5     : Ratio of outer radius vs. inner radius, i.e. outer/inner.
! R1     : Magnitude of current on probe.
! R2     : Phase of current on probe (in degrees).
! R3,R4,R5 : (x,y,z) position of coax aperture center (m). 
! R6     : Outer coax radius (m).
! R7     : Center conductor length and sense (+ or -) 
!***************************************************

  analysis_counter = 1
  AX_counter = AX_counter+1
  CALL READ_LINE(string,AXdata(AX_counter)%anfl,itemp2,itemp3,AXdata(AX_counter)%coaxdir,  &
         itemp5, AXdata(AX_counter)%coax_Iabs, AXdata(AX_counter)%coax_Iphase,             &
         AXdata(AX_counter)%coaxcentre(1), AXdata(AX_counter)%coaxcentre(2),               &
         AXdata(AX_counter)%coaxcentre(3), AXdata(AX_counter)%coax_b,                      &
         AXdata(AX_counter)%coaxlen)

  ! Check that the direction is valid:
  SELECT CASE(AXdata(AX_counter)%coaxdir)
  CASE(1:3)
    CONTINUE ! Valid data
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4061,int1=AX_counter)
  END SELECT

  ! The length indicates sense, therefore must not be =0:
  IF (AXdata(AX_counter)%coaxlen.EQ.0.0) THEN
    CALL ERROR_FEMFEKO(1,4062,int1=AX_counter)
  END IF

  ! Calculate the normal vector, pointing away from the FEM volume:
  AXdata(AX_counter)%normal = 0.0 ! vector assignment
  AXdata(AX_counter)%normal(AXdata(AX_counter)%coaxdir) = 1.0
  IF (AXdata(AX_counter)%coaxlen.GT.0.0) THEN
    AXdata(AX_counter)%normal = - AXdata(AX_counter)%normal
  END IF

  ! Calculate coax end point:
  AXdata(AX_counter)%coaxend = AXdata(AX_counter)%coaxcentre
  AXdata(AX_counter)%coaxend( AXdata(AX_counter)%coaxdir ) = &
    AXdata(AX_counter)%coaxend( AXdata(AX_counter)%coaxdir ) + AXdata(AX_counter)%coaxlen

  ! Calculate inner radius:
  AXdata(AX_counter)%coax_a = AXdata(AX_counter)%coax_b/REAL(itemp5)

  ! Specify coax eps_r (a fixed internal program setting):
  AXdata(AX_counter)%coax_eps = 1.0
  AXdata(AX_counter)%coax_mu  = 1.0

END SUBROUTINE AX_CARD
!*******************************************************************************


SUBROUTINE BC_CARD
!***************************************************
! The BC card: FEMFEKO Boundary Condition card. Optional.
!***************************************************
! Format I1 I2 (I3) (I4) (I5) R1 R2 R3 (R4)
!***************************************************
! I1 = 0  : Free - i.e. not set.
!    = 1  : Port. Excitation value defined on A9 card.
!    = 2  : PEC BC
!    = 3  : Reserved for PMC BC.
!    = 4  : Reserved for Perodic BC.
!    = 5  : Absorbing boundary condition 
! I2      : Port or ABC label (unique for all ports, and for all ABC, 
!           but a port and ABC may have the same label).
! I3, I4 & I5 not presently used, but should have zero "place-holders".
! R1, R2, R3: coordinates of point 1 of quadrilateral at which BC is to be
! enforced followed by three further lines, with I1-5 unused but with zero
! placeholders, with subsequent R1, R2 and R3 values specifying coordinates
! of corners 2, 3, and 4 on lines 2, 3 and 4 respectively.
!
! Further parameters of ABC are set on the AB card.
!
! Note that both port and ABC boundary conditions may only
! be set on the boundary (closure) of the FEM mesh at present. 
!
!***************************************************
  BC_counter = BC_counter+1
  CALL READ_LINE(string,itemp,itemp1,itemp2,itemp3,itemp4, & 
    xyztemp1(1),xyztemp1(2),xyztemp1(3))
  IF (itemp.EQ.1) THEN ! this is a port
    BCs(BC_counter)%port_num = itemp1 ! assign port number
  END IF
! Added DBD 26 March 2003
  IF (itemp.EQ.5) THEN ! this is an ABC
	! Assign default admittance of free space
    BCs(BC_counter)%ABC_num = itemp1 ! assign ABC label
  END IF
! End added DBD 26 March 2003
  READ (MESHIN,*)itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp2
  READ (MESHIN,*)itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp3
  READ (MESHIN,*)itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp4
   BCs(BC_counter)%type = itemp
  SELECT CASE(itemp)
  CASE(0)    
    CONTINUE ! Valid BC specification - free, ie. not set.
  CASE(1) ! Valid BC specification - Port
    CONTINUE
  CASE(2)    
    CONTINUE ! Valid BC specification - PEC 
  CASE(5)    
    CONTINUE ! Valid BC specification - PEC 
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4007,int1=BC_counter)
  END SELECT
  BCs(BC_counter)%corner(1,1:3) = xyztemp1(1:3) ! These are the coordinates
  BCs(BC_counter)%corner(2,1:3) = xyztemp2(1:3) ! of the four corners of the 
  BCs(BC_counter)%corner(3,1:3) = xyztemp3(1:3) ! quadrilateral.
  BCs(BC_counter)%corner(4,1:3) = xyztemp4(1:3) 
END SUBROUTINE BC_CARD

!*******************************************************************************

SUBROUTINE BO_CARD
!***************************************************
! The BO card: FEKO ground plane card. 
!***************************************************
! Format I1
!***************************************************
! I1 = 0 : No ground plane using reflection coefficient approximation (RCA) used.
!    = 1 : Use ground plane with specified parameters R1..R4.
!    = 2 : PEC ground plane at z=0
!    = 3 : PEC ground plane at z=0
! NB! Only option 2 relevant for FEM at present! 
! Note that FEKO can also use R1..R4 parameters, but not FEMFEKO.
!***************************************************
  CALL READ_LINE(string,itemp1)
  SELECT CASE(itemp1)
  CASE(2) 
    PEC_ground = .true.
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4008,int1=itemp1)
  END SELECT
  
END SUBROUTINE BO_CARD

!*******************************************************************************

! DBD addition 28 Mar 03
SUBROUTINE DP_CARD
!***************************************************
! The DP card. A partial implementation of the full
! FEKO DP card. 
!***************************************************
!Format S1 (S2) (S3) (S4) (S5) R1 R2 R3
!***************************************************
! S1  : Name (max 5 characters)
! R1  : x coordinate in m
! R2  : y coordinate in m
! R3  : z coordinate in m
!***************************************************
  DP_counter =   DP_counter + 1
  CALL READ_STRING_LINE(string,stemp1,stemp2,stemp3,stemp4,stemp5,rtemp1,rtemp2,rtemp3)
  DPoints(DP_counter)%name = stemp1
  DPoints(DP_counter)%coords(1) = rtemp1
  DPoints(DP_counter)%coords(2) = rtemp2
  DPoints(DP_counter)%coords(3) = rtemp3
END SUBROUTINE DP_CARD
! End DBD addition 28 March 03
!*******************************************************************************

SUBROUTINE EA_CARD
!***************************************************
! The EA card: The FEMFEKO Eigenanalysis ARPARCK 
! control card. Required if sparse format used for 
! eigenanalysis.
!***************************************************
! Format I1 (I2) (I3) (I4) (I5) R1
!***************************************************
! I1 = nev: Number of eigenvalues requested.
! I2 = ncv: Number of Arnoldi vectors used. 
! I3 - Spectrum control:
!     = 1: 'SM' (Smallest magnitude)
!     = 2: 'LM' (Largest magnitude)
!     = 3: 'BE' (Both ends)
!     = 4: 'SA' (Smallest amplitude - not used)
!     = 5: 'LA' (Largest amplitude - not used)
! I4 = maxitr: Maximum number of iterations permitted.
! I5 = mode: ARPACK mode (2 to 5)
! R1 = sigma (Spectrum shift)
!***************************************************

  CALL READ_LINE (string,nev,ncv,itemp,maxitr,mode,sigma)
  SELECT CASE(itemp)
  CASE(1)
    which = 'SM'
  CASE(2)
    which = 'LM'
  CASE(3)
    which = 'BE'
  CASE(4)
    which = 'SA'
  CASE(5)
    which = 'LA'
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4009)
  END SELECT
  
  IF (maxitr.LT.500) THEN
    maxitr = 500 
    CALL ERROR_FEMFEKO(0,4010)
  END IF
 
  IF ((mode.GT.5).or.(mode.LT.2)) THEN
    CALL ERROR_FEMFEKO(0,4011) 
  END IF
    
  IF (sigma.LT.0) THEN
    CALL ERROR_FEMFEKO(1,4012)
  END IF 

 !  Check also that if SPARSE and EIGENANALYSIS, that EA has been defined. 
  EA_card_defined = .TRUE.  ! this is tested at the end of meshin to ensure that the card was indeed defined.

END SUBROUTINE EA_CARD

SUBROUTINE EI_CARD
!***************************************************
! The EI card: The FEMFEKO EIgenanalysis control card. Optional.
!***************************************************
! Format I1 (I2) (I3) (I4) (I5) R1
!***************************************************
! I1 < 0 : No eigenmodes computed.
!    = 0 : First eigenmode to be computed, as follows: 
!          Either the eigenmode with eigenvalue closest to R1 will be the first
!          computed, or if I3 is set, the first non-"trivial" (non-zero) 
!          eigenmode will be the first computed.
!    > 0 : The eigenmode corresponding to eigenvalue I1 will the first computed.
! I2 = x : x eigenmodes, starting as specified in I1, will be computed. 
! I3 = 1 : A theoretical prediction is to be used for the number of "trivial"
!          (zero-approximation) eigenvalues.    
!    <> 1: Either R1 is used (if I1=0), otherwise no effect.
! I4 & I5 not presently used, but should have zero "place-holders"
! R1     : An estimate of the eigenvalue of the first eigenmode required.
!          (only if I1=0 and I3 <>1, otherwise any value). 
!***************************************************
  USE eigen_analysis_data

  analysis_counter = 1
  CALL READ_LINE (string,itemp1,itemp2,itemp3,itemp4,itemp5,rtemp1)
  SELECT CASE(itemp1)
  CASE(:-1) 
    COMPUTE_EIGENMODES = .false.
  CASE(0:) 
    COMPUTE_EIGENMODES = .true.
    first_eigenmode = itemp1 ! Corrected in EIGEN_SYSMAT if theoretical
                             ! estimate used.
  END SELECT
  
  num_eigenmodes = itemp2
  IF(itemp3.EQ.1) THEN 
    PREDICT_SPURIOUS_EIGENMODES = .TRUE.
  ELSE 
    PREDICT_SPURIOUS_EIGENMODES = .FALSE.
  END IF
  IF (first_eigenmode.EQ.0.AND..NOT.PREDICT_SPURIOUS_EIGENMODES) THEN
    first_eigen_approx = rtemp1
  END IF


END SUBROUTINE EI_CARD

!*******************************************************************************

SUBROUTINE FA_CARD
!***************************************************
! The FA card: FEMFEKO Field Averaging card. Optional.
!***************************************************
! Format I1 I2 I3 I4 I5 R1 R2 R3 R4 R5 R6
!***************************************************
! I1 = 1 : Field in FE region computed without any averaging. Default
!    = 2 : Field in FE region computed using simple averaging.
! I2, I3, I4 & I5 not presently used, but should have zero "place-holders"
! R1     : the tolerance to use for averaging.
! R2, R3, R4, R5, R6 not used. 
!***************************************************

  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5, & 
                 eps_field_averaging)
  SELECT CASE(itemp1)
  CASE(1)
    FIELD_AVERAGING = .false.
  CASE(2) 
    FIELD_AVERAGING = .true.
  CASE DEFAULT    
    CALL ERROR_FEMFEKO(1,4013,int1=itemp1)
  END SELECT
  
END SUBROUTINE FA_CARD

!*******************************************************************************

SUBROUTINE FE_CARD
  USE scattering_analysis_data
!***************************************************
! The FE card: FEKO near FiEld  card: Optional for FEA. Same format as FEKO 
! card, but with new numeric option for FELTYP.
! [DBD 02 Apr 2003] 
! Note that for time domain analysis, the output 
! format still reports "phase" but this should 
! be zero and no meaning should be attached to this. 
!***************************************************
! Format I1 I2 I3 I4 I5 R1 R2 R3 R4 R5 R6
!***************************************************
! I1 = +8 : Total field in FE region computed.
!      -8 : Scattered field in FE region computed (only limited availability).
!
! I2,I3&I4: number of points in x,y and z direction respectively.
! I5      : Coordinate system used (0 for Cartesian, only one presently
!           implemented in FEMFEKO.
!
! R1,R2&R3: x,y and z coordinates respectively of first point.
! R4,R5,R6: Delta x, y and z respectively. 
!
!***************************************************

  FE_counter = FE_counter+1
  CALL READ_LINE(string,FEdata(FE_counter)%feltyp,FEdata(FE_counter)%n_x,         &
         FEdata(FE_counter)%n_y,FEdata(FE_counter)%n_z,FEdata(FE_counter)%felkor, &
         FEdata(FE_counter)%x0,FEdata(FE_counter)%y0,FEdata(FE_counter)%z0,       &
         FEdata(FE_counter)%delta_x,FEdata(FE_counter)%delta_y,                   &
         FEdata(FE_counter)%delta_z)
  SELECT CASE(FEdata(FE_counter)%feltyp)
  CASE(-8)
    IF(.NOT.(TD_ANALYSIS.AND.SCAT_FIELD)) CALL ERROR_FEMFEKO(1,4016,int1=FEdata(FE_counter)%feltyp,int2=FE_counter)
    CONTINUE ! Valid data.
  CASE(8)
    CONTINUE ! Valid data.
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4016,int1=FEdata(FE_counter)%feltyp,int2=FE_counter)
  END SELECT

  SELECT CASE(FEdata(FE_counter)%felkor)
  CASE(0)
    CONTINUE ! Cartesian coordinates
  CASE DEFAULT  
      CALL ERROR_FEMFEKO(1,4017,int1=FEdata(FE_counter)%felkor,int2=FE_counter)
  END SELECT

  ! Assign the analysis number:
  IF (analysis_counter.GT.0) THEN
    FEdata(FE_counter)%analysis_no = analysis_counter

  ELSE
    CALL ERROR_FEMFEKO(1,4018,int1=FE_counter)
 END IF

END SUBROUTINE FE_CARD

!*******************************************************************************

SUBROUTINE FF_CARD
!***************************************************
! The FF card: FEKO Far Field  card: Optional for FEA. Same format as FEKO 
! card.
!***************************************************
! Format I1 I2 I3 I4 (I5) R1 R2 R3 R4
!***************************************************
! I1 = 0 : No far field calculation done
!    = 1 : Far field calculated using following parameters.
!    = 2 : Far field calculated in same direction as incident plane wave.
! I2     : Number of observation points in theta direction.
!          Defaults to 1 as in FEKO.
! I3     : Number of observation points in phi direction.
!          Defaults to 1 as in FEKO
! I4 = 0 : Directivity calculated.
!    = 1 : Gain calculated.
! I5 not presently used, but should have zero "place-holders"
!
! R1     : Theta_start
! R2     : Phi_start
! R3     : Delta_theta
! R4     : Delta_phi
!***************************************************

  FF_counter = FF_counter+1
  CALL READ_LINE(string,FFdata(FF_counter)%FFreq,FFdata(FF_counter)%NTheta,       &
         FFdata(FF_counter)%NPhi,FFdata(FF_counter)%rige,itemp5,                  &
         FFdata(FF_counter)%Theta0,FFdata(FF_counter)%Phi0,                       &
         FFdata(FF_counter)%DTheta,FFdata(FF_counter)%DPhi)
  SELECT CASE(FFdata(FF_counter)%FFreq)
  CASE(0:2)
    CONTINUE ! Valid data.
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4019,int1=FF_counter)
  END SELECT
  IF (FFdata(FF_counter)%NTheta.EQ.0) THEN 
    FFdata(FF_counter)%NTheta = 1 ! FEKO default
  END IF
  IF (FFdata(FF_counter)%NPhi.EQ.0) THEN 
    FFdata(FF_counter)%NPhi = 1     ! Ditto
  END IF
  SELECT CASE(FFdata(FF_counter)%rige)
  CASE(0:1)
    CONTINUE ! Valid data.
  CASE DEFAULT  
    CALL ERROR_FEMFEKO(1,4020,int1=FF_counter)
  END SELECT
  
  ! Assign the analysis number:
  IF (analysis_counter.GT.0) THEN
    FFdata(FF_counter)%analysis_no = analysis_counter
  ELSE
    CALL ERROR_FEMFEKO(1,4021,int1=FF_counter)
  END IF

END SUBROUTINE FF_CARD

!*******************************************************************************


SUBROUTINE FM_CARD
!***************************************************
! The FM card: FEMFEKO card to control the 
! Finite element Matrix storage and solution. Optional.
! See also EA card for eigenanalysis options.
!***************************************************
! Format I1 I2 I3 (I4) (I5) R1 R2 R3 R4
!***************************************************
! I1 = Renumbering undertaken to minimize matrix bandwidth (obsolete)
!    <>1: Not thus renumbered (current default).  
! I2 = 0: Matrix stored as sparse matrix (current default).
!    = 1: Matrix stored as full format.
!    = 2: Matrix stored in banded format (obsolete).
! I3      Solution scheme for Ax=b (defaults: LU for full, Bi-CG for sparse):
!    = 0: LU factorization (only available for sparse matrices for TD analysis at present)
!    = 1: Iterative solver - Bi-conjugate gradient.
!    = 2: Iterative solver - Conjugate gradient.
!    = 3: QMR
!    = 4: GMRES
!         NB! Does not effect eigenanalysis - see EA card for this.
!    = 5: Experimental preconditioner using mixed potentials. 
! I4 = 0: Matrix elements scaled by edge lengths (current default).
!    <>0: Not thus scaled.
! I5      On-screen solution progress monitor for iterative solver.
!    = 0: No on-screen reporting (current default).
!      1: On-screen reporting of iteration progress.
! R1    : The residual norm used by iterative linear algebra routines.
!         as a stopping criteria. If iterative solvers not used, 
!         an arbitrary value (eg 0.0) is still required.
! R2    : Factor multiplying degrees of freedom to give maximum number
!         of iterations. (Rounded off to an integer subsequently).
!         A value is required: 2.0 is adequate for most purposes.
! R3    : Flag for selecting pre-conditioner (>0) or no preconditioner
!       : <= 0 (default)
!         A value (eg 0.0) is required even if pre-conditioning not required. 
! R4    : Number of cycles before restart (GMRES). Specify as real, but 
!         converted to integer in code. 
!         A value (eg 0.0) is required even if GMRES is not used. 
!***************************************************
! DBD addition 25 April 2001
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5,& 
                 rtemp1,rtemp2,rtemp3,rtemp4)
! End DBD addition 25 April 2001
  SELECT CASE(itemp1)
  CASE(1)
    BANDRENUM = .true.
  CASE DEFAULT 
    BANDRENUM = .false.
  END SELECT
  
  SELECT CASE(itemp2)
  CASE(0)
    SPARSE = .true.
  CASE(1)
    SPARSE = .false.
  CASE DEFAULT 
    CALL ERROR_FEMFEKO(1,4022)
  END SELECT

  SELECT CASE(itemp3)
! changed by MMB, 16 Apr 2001
! was 
!  CASE(0:2)
! and then 
! CASE(0:4)
! end MMB change 
! and is now
  CASE(0:5)
! and DBD extension
    SOLVER_TYPE = itemp3
  CASE DEFAULT 
    CALL ERROR_FEMFEKO(1,4023)
  END SELECT

  SELECT CASE(itemp4)
  CASE(0)
    SCALE_BY_EDGE_LENGTH = .true.
  CASE DEFAULT 
    SCALE_BY_EDGE_LENGTH = .false.
  END SELECT

  SELECT CASE(itemp5)
  CASE(1)
    ON_SCREEN_REPORTING = .true.
  CASE DEFAULT 
    ON_SCREEN_REPORTING = .false.
  END SELECT
 
  residual_norm = rtemp1
  max_iter_factor = rtemp2
! DBD addition 25 April 2001
  IF(rtemp3.GT.0.0_SP) THEN 
    USE_PRE_CONDITIONER = .TRUE.
  ELSE 
    USE_PRE_CONDITIONER = .FALSE.
  END IF
! End DBD addition 25 April 2001

! DBD addition 27 April 2001
  restart_GMRES = INT(rtemp4)
! End DBD addition 27 April 2001


END SUBROUTINE FM_CARD

!*******************************************************************************


! Changed DBD 08 March 2002.
SUBROUTINE GW_CARD
!***************************************************
! The GW card: FEMFEKO card to control various
! parameters of the Guided Wave analysis module.
! Optional.
!***************************************************
! Format I1 (I2) (I3) (I4) (I5)
!***************************************************
! I1    : Controls type of analysis.
!       : 0 (Default) Full n-port analysis undertaken.
!       : 1 Excitation for each port used.
! I2    : Controls numerical integration accuracy over port.
!       : 0 Low-order scheme used.
!       : 1 (Default) Moderate-order scheme used.
!       : 2 High-order scheme used.
!       : 3 Very high-order scheme used.
!       : 4 Extremely high-order scheme used.
! I3 and I4:
!       Multi-mode analysis; number of modes in width (m) and height (n)

!       respectively. Set to 1 and 0 for single mode TE10 analysis.
! I5    : Presently unused, but should have a zero "place holder".
!***************************************************
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5)
  SELECT CASE(itemp1)
  CASE(0)
    GW_COMPUTE_S_PARAMS = .TRUE.
  CASE(1)
    GW_COMPUTE_S_PARAMS = .FALSE.
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4037,int1=itemp1)
  END SELECT
  SELECT CASE(itemp2)
  CASE(-1)
    gauss_points = 1
  CASE(0)
    gauss_points = 3
  CASE(2)
    gauss_points = 6
  CASE(3)
    gauss_points = 13
  CASE(4)
    gauss_points = 25
  CASE DEFAULT ! Includes case 1
    gauss_points = 4
  END SELECT
  num_wg_modes_m = itemp3
  num_wg_modes_n = itemp4
END SUBROUTINE GW_CARD
!*******************************************************************************
! End DBD changes 08 March 2002


! DBD addition 20 May 2003 - presently removed, not needed. 
!SUBROUTINE KU_CARD
!***************************************************
! The KU card. A partial implementation of the full
! FEKO KU card. Not all parameters are required. 
! Used with the TD card to speficy an interior surface
! Sprime containing all inhomogenities and 
! within which only the scattered fields are computed.
! Only one permitted.
!
!***************************************************
!Format S1 S2 S3 (S4) (S5) R1 R2 R3 R4 R5 
!***************************************************
! S1-S3: See KU Card in FEKO manual 
! S4: unimplemented, but should have zero "place holder".
!     Normal must be outward pointing in FEMFEKO.
! S5: unused, but should have zero "place holder".
! R1-R4: See KU Card in FEKO manual. Not used at present.
! R5: unimplemented, but should have zero "place holder"
!***************************************************
!  ! Test if a QU card has already been read in.
!  IF (KU_card_defined) THEN
!    CALL ERROR_FEMFEKO(1,4096)
!  END IF
!  ! Set the relevant flag. 
!  KU_card_defined = .true.
!  CALL READ_STRING_LINE(string,stemp1,stemp2)
!  Spherical_Boundary%S1  = stemp1
!  Spherical_Boundary%S2  = stemp2
!  Spherical_Boundary%S3  = stemp2
!END SUBROUTINE KU_CARD
! End DBD addition 20 May 2003
!*******************************************************************************


SUBROUTINE MA_CARD
!***************************************************
! The MA card: FEMFEKO card for the MAterial parameters.
! Optional (all materials treated as free space if 
! not explicitly specified. A warning is printed about
! any such implicitly specified materials). 
!***************************************************
! Format I1 I2 (I3) (I4) (I5) R1 R2 R3 R4
!***************************************************
! I1 = x: All elements with label (layer in FEMAP) x will be assigned
!         the following material parameters. Note that I1 >=0.  
!
! I2 = x: x<0 - Default. Material not defined.
!         x=0 - Element consists of PEC. All d.o.f.'s set to zero
!               Reserved for use, but not yet used in code.
!         x=1 - Isotropic, possibly lossy, dielectric or magnetic materials
!         x=2 - Anisotropic, possibly losssy, dielectric or magnetic materials 
!               with DIAGONAL permittivity or permeability tensors
!         x=3 - PML. Used internally by code, not user-accessible.
! I3,I4,I5: unused, but should have zero "place holders".
! R1, R2, R3, R4: eps_r (Re & Im); mu_r(Re & Im)
! For anisotropic materials: another 8 lines, with the   
! full permittivity and permeability matrix defined with eps_r and mu_r 
! in the sequence
! xx (on same line as MA card) xy xz yx yy yz zx zy zz 
! per line
! 
! Note that the material layer numbers do not have to be numbered from
! 0, and also can "skip" numbers, but should be kept as small as possible. 
! The limit to the maximum number of definable materials is dynamically
! calculated. 
! Note that the materials need not be listed in consecutive order
! (i.e. materials can be defined in the sequence 2,5,3, for example.
! but this obviously can waste memory, since space for the "unused"
! material labels is also allocated - such materials are automatically
! assigned the properties of free space.
!
! Added DBD 25 March 2003
! For time domain analysis, the real part of eps_r or mu_r 
! is the relative permittivity or permeability respectively; 
! the "imaginary" part is used to store the electric or magnetic conductivity
! sigma or sigma_m respectively. Note that this is only for coding convenience;
! this is NOT the imaginary part of the permittivity which would be sigma/(j omega)
!***************************************************

  MA_counter = MA_counter + 1
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5,&
    parms(1),parms(2),parms(3),parms(4))
  IF(itemp1.GT.MAX_MATERIALS) THEN ! Error check - this shouldn't happen.
    STOP 'IE: MAX_MATERIALS: somehow exceeded internally.'
  ELSE IF (itemp1.GT.highest_material) THEN 
     highest_material = itemp1 
  END IF
  IF (itemp1.LT.0) THEN
    CALL ERROR_FEMFEKO(1,4025,int1=itemp1)
  END IF
  IF (MA_counter.EQ.1) THEN
     lowest_material = itemp1  
  ELSE IF (itemp1.LT.lowest_material) THEN
     lowest_material = itemp1
  END IF    
  
  MAdata(MA_counter)%label         = itemp1
  MAdata(MA_counter)%material_type = itemp2
  MAdata(MA_counter)%analysis_no   = analysis_counter + 1  ! applicable to all 
                                                           ! subsequent analysis


  SELECT CASE(MAdata(MA_counter)%material_type) 
  CASE(0) ! PEC 
  CASE(1) ! Isotropic material
    MAdata(MA_counter)%eps_r = CMPLX(parms(1),parms(2))
    MAdata(MA_counter)%mu_r  = CMPLX(parms(3),parms(4))
  CASE(2) ! Anisotropic material
          ! Defined in sequence xx xy xz yx yy yz zx zy zz on MA cards
    MAdata(MA_counter)%eps_r_tensor(1,1) = CMPLX(parms(1),parms(2))
    MAdata(MA_counter)%mu_r_tensor(1,1)  = CMPLX(parms(3),parms(4))
    READ(MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,parms
    MAdata(MA_counter)%eps_r_tensor(1,2) = CMPLX(parms(1),parms(2))
    MAdata(MA_counter)%mu_r_tensor(1,2)  = CMPLX(parms(3),parms(4))
    READ(MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,parms
    MAdata(MA_counter)%eps_r_tensor(1,3) = CMPLX(parms(1),parms(2))
    MAdata(MA_counter)%mu_r_tensor(1,3)  = CMPLX(parms(3),parms(4))
    DO i_row=2,3 ! Read other six entries
      DO i_col=1,3
        READ(MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,parms             
        MAdata(MA_counter)%eps_r_tensor(i_row,i_col) = CMPLX(parms(1),parms(2))
        MAdata(MA_counter)%mu_r_tensor(i_row,i_col)  = CMPLX(parms(3),parms(4))
      END DO
    END DO
  CASE DEFAULT ! An error condition
    CALL ERROR_FEMFEKO(1,4026,int1=itemp1)
  END SELECT

END SUBROUTINE MA_CARD

!*******************************************************************************


SUBROUTINE PL_CARD
!***************************************************
! The PL card: FEMFEKO Perfectly matched Layer card. 
! Required to define a PML. Note that this card is 
! independent of the AB card - a PML region may 
! be terminated by an AB card, or left unterminated (in
! which case a PEC termination is usually assumed). 
! A rectangular external boundary is ASSUMED!!
!
! The PML characteristics are the same for x,y and z 
! absorption. However, absorption in each of these 
! directions may be selected or not. In overlap regions,
! the usual approach of absoption in both (or all three)
! directions is used.
!
! For free space scattering, absorption in all directions
! would usually be chosen and this option is primarily for testing. 
!
! Two PML's are available. The first (I5=1) is based on an 
! approximate hyperbolic approximation
! of the ideal PML. The approximation improves as frequency
! increases. The second (I5=6) is a numerically exact implementation
! as proposed by Jiao and Jin, IEEE T-AP, Feb 2003. 
!
!***************************************************
! Format I1 I2 I3 (I4) (I5) R1 (R2) R3 R4 R5
!***************************************************
! I1      : polynomial order of PML
! I2      =  1: Apply to + and - \hat{x} travelling waves. 
!         <> 1: Do not apply. 
! I3      =  1: Apply to + and - \hat{y} travelling waves. 
!         <> 1: Do not apply. !
! I4      =  1: Apply to + and - \hat{z} travelling waves. 
!         <> 1: Do not apply. 
!           Note: I2,I3 and I4 may be combined as desired; 
!           usually, all will be 1 for PML all around. 
! I5      : Full (=1) or approximate (<>1) treatment.
! R1      : Thickness of PML.
! R2      : Maximum conductivity.
!
!***************************************************
  ! Test is a PL card has already been read in.
  IF (PL_card_defined) THEN
    CALL ERROR_FEMFEKO(1,4097)
  END IF
  ! Set the relevant flag. 
  PL_card_defined = .TRUE.
  PML_PRESENT     = .TRUE. ! Actually already set.
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5,rtemp1,rtemp2) 

  PML%m= itemp1
  IF(itemp2.EQ.1) PML%absorb_x = .TRUE.
  IF(itemp3.EQ.1) PML%absorb_y = .TRUE.
  IF(itemp4.EQ.1) PML%absorb_z = .TRUE.
  IF(itemp5.EQ.1) PML%full     = .TRUE.

  PML%thickness = rtemp1
  PML%sigma_max = rtemp2

END SUBROUTINE PL_CARD
!*******************************************************************************



! DBD addition 28 Mar 03
SUBROUTINE QU_CARD
!***************************************************
! The QU card. A partial implementation of the full
! FEKO QU card. Not all parameters are required. 
! Used with the TD, BC and AB cards to specify
! a rectangular box. Only one permitted.
! 
! S1 should be closer to the origin than S2. 
!
!***************************************************
!Format I1 I2 
!***************************************************
! I1    : S1 See QU Card in FEKO manual 
! I2    : S2 Name of opposite corners.
!***************************************************
  ! Test if a QU card has already been read in.
  IF (QU_card_defined) THEN
    CALL ERROR_FEMFEKO(1,4088)
  END IF
  ! Set the relevant flag. 
  QU_card_defined = .true.
  CALL READ_STRING_LINE(string,stemp1,stemp2)
  ABC_Box%S1  = stemp1
  ABC_Box%S2  = stemp2
END SUBROUTINE QU_CARD
! End DBD addition 28 March 03
!*******************************************************************************


! DBD addition 24 Mar 03
! Extended DBD 8 May 03
SUBROUTINE TD_CARD
  USE scattering_analysis_data
!***************************************************
! The TD card: The FEMFEKO Time Domain analysis card.
! Optional. Frequeny domain analysis performed 
! unless this card is defined. 
! Only one of these cards may be defined. 
! Note that amplitude of plane wave is defined in A0
! card as usual, as well as direction of polarization and 
! propagation.
! It is presently assumed that accompanying this 
! card are six ABC's numbered 1-6 and forming 
! a rectangular box aligned with the x-y-z axes. 
! This must be defined with an appropriate QU card. 
!
! If the scattered field formulation is used, an 
! additional internal spherical surface Sprime
! must be defined using the KU card. This must 
! be the same as one defined in FEMAP. This cannot be checked.
!
!***************************************************
!Format I1 I2 (I3 I4 I5) R1 R2 R3 R4 R5
!***************************************************
!         Analysis type
! I1 = 1: Total field solver. Field introduced at boundary.
!    = 2: Scattered field within interior surface Sprime. 
!
!         Source type 
! I2 = 1: Differentiated Gaussian source.
! I3, I4&I5: Unused, but zero placeholders required. 
! R1    : Number of time steps (converted to integer)
!         Done to allow a very large number of time steps
!         to be defined. 
! R1    : delta t
! R2    : Newmark parameter beta
!         For Gaussians sources:
! R3    : sigma (in seconds)
! R4    : time offset (in seconds)
!
!***************************************************
  ! Test is a TD card has already been read in.
  IF (TD_card_defined) THEN
    CALL ERROR_FEMFEKO(1,4083)
  END IF
  ! Set the relevant flag. 
  TD_card_defined = .true.
  CALL READ_LINE(string,itemp1,itemp2,itemp3,itemp4,itemp5,& 
                 rtemp1,rtemp2,rtemp3,rtemp4,rtemp5)
  SELECT CASE(itemp1)
  CASE(1) 
    SCAT_FIELD = .FALSE.
  CASE(2) 
    SCAT_FIELD = .TRUE.
  CASE DEFAULT
    CALL ERROR_FEMFEKO(1,4084)
  END SELECT
  TD_source_type = itemp2
  num_timesteps  = REAL(rtemp1)
  delta_t        = rtemp2
  Newmark_beta   = rtemp3
  sigma_pulse    = rtemp4
  offset_pulse   = rtemp5

END SUBROUTINE TD_CARD
! End DBD addition 24 March 03 and extensions 08 May 03
!*******************************************************************************

SUBROUTINE TE_CARD
!***************************************************
! The TE card : FEMFEKO card for each TEtrahedral element
! Required for each tetrahedral element.
!***************************************************
! Format I1 (I2) (I3) (I4) (I5) R1 R2 R3 
!***************************************************
! I1 = x: Label for this element
! I2-5: Not used, but should have zero "placeholders"
! R1, R2, R3: coordinates of node 1.
! followed by three further lines, with I1-5 unused but with placeholders,
! with subsequent R1, R2 and R3 values specifying coordinates of 
! nodes 2, 3, and 4 on lines 2, 3 and 4 respectively.
!***************************************************
   CALL READ_LINE(string,itemp,itemp1,itemp2,itemp3,itemp4, & 
     xyztemp1(1),xyztemp1(2),xyztemp1(3))
   READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp2
   READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp3
   READ (MESHIN,*) itemp1,itemp2,itemp3,itemp4,itemp5,xyztemp4

   !***************************************************
   !Debugging section : check that data has been read correctly...
   !print *, string1, element_counter, itemp,xyztemp1
   !print *, xyztemp2
   !print *, xyztemp3
   !print *, xyztemp4
   !***************************************************

   xyz_pos(1,1:3) = xyztemp1(1:3); !These are the vertex coordinates: 
   xyz_pos(2,1:3) = xyztemp2(1:3); !three in the case of a triangle,  
   xyz_pos(3,1:3) = xyztemp3(1:3); !and four in the case of a tetrahedron  
   xyz_pos(4,1:3) = xyztemp4(1:3); 
   !debugging********************************************
    !WRITE(FILEOUT,*) '***TE card ***************************'
    !WRITE(FILEOUT,*) 'xyz_pos : '
    !DO i = 1,4 
    ! WRITE(FILEOUT,'(T10,F19.4, T30, F19.4, T50, F19.4)') xyz_pos(i,:)
    !END DO
   !******************************************************
   element_counter = element_counter+1 ! counts the number of TE cards, or elements
   elements(element_counter)%material = itemp
   ! Following is checked later to ensure that materials (labels) referenced
   ! here are indeed defined.
   IF(itemp.GT.highest_material_referenced) THEN
     highest_material_referenced = itemp
   END IF
   VERTEX_LOOP: DO i_vert = 1,NUM_VERTICES 
     new_node = .TRUE.  !node has not been entered in vertices yet

     IF (node_counter.GT.2) THEN !The first four nodes must be new
       NODE_LOOP: DO i_node = 1,node_counter-1
         IF((ABS(xyz_pos(i_vert,1)-vertices(i_node)%coord(1)).LT.EPS).AND.  & 
            (ABS(xyz_pos(i_vert,2)-vertices(i_node)%coord(2)).LT.EPS).AND.  & 
            (ABS(xyz_pos(i_vert,3)-vertices(i_node)%coord(3)).LT.EPS) ) THEN
           new_node = .FALSE.  !this node has been entered in the list of nodes already
           ipos = i_node       !this is the position of the node in vertices   
           EXIT NODE_LOOP
         END IF
       END DO NODE_LOOP 
     END IF
     IF (new_node)  THEN
       elements(element_counter)%nodes(i_vert) = node_counter
       vertices(node_counter)%coord = xyz_pos(i_vert,1:3)! enter node in list
       node_counter = node_counter +1
     ELSE 
       elements(element_counter)%nodes(i_vert) = ipos 
                                ! The node is already in node list; 
                                ! enter the number of the node in element list.
     END IF
  END DO VERTEX_LOOP
  num_nodes = node_counter
END SUBROUTINE TE_CARD


! Obsolete MIPS workaround removed 28 March 2003 DBD.
! Extended DBD 20 March to read 8 real parameters.
SUBROUTINE READ_LINE(string,int1,int2,int3,int4,int5,&
 real1,real2,real3,real4,real5,real6,real7,real8)
!***************************************************
! This routine was required due to a compiler problem
! under MIPS 7.1 under Irix 6.3, which was unable to 
! correctly compile ADVANCE='NO' on the read statements, 
! forcing a clumsy workaround on subsequent reads.
! (Option removed DBD 28 March 2003). 
! This routine  essentially groups much of the workaround. 
! Note that on multi-line cards (eg TE), a blank string may
! be read on subsequent lines using the work-around: this is irrelevant.
!***************************************************

  CHARACTER(15), INTENT(out) :: string

  INTEGER(I4B), INTENT(out), OPTIONAL  :: int1,int2,int3,int4,int5
      
  REAL(SP), INTENT(out), OPTIONAL  :: real1,real2,real3,real4,real5,real6,real7,real8

  string  = ' '  ! This paramater is now redundant, but would require too many changes to 
                 ! code to fix easily.  

  IF(PRESENT(real8)) THEN
     READ(MESHIN,*)int1,int2,int3,int4,int5,real1,real2,real3,real4,real5,real6,real7,real8
  ELSE IF(PRESENT(real7)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2,real3,real4,real5,real6,real7
  ELSE IF(PRESENT(real6)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2,real3,real4,real5,real6
  ELSE IF(PRESENT(real5)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2,real3,real4,real5
  ELSE IF(PRESENT(real4)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2,real3,real4
  ELSE IF(PRESENT(real3)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2,real3
  ELSE IF(PRESENT(real2)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1,real2
  ELSE IF(PRESENT(real1)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5,real1
  ELSE IF(PRESENT(int5)) THEN
     READ(MESHIN,*) int1,int2,int3,int4,int5
  ELSE IF(PRESENT(int4)) THEN
     READ(MESHIN,*) int1,int2,int3,int4
  ELSE IF(PRESENT(int3)) THEN
     READ(MESHIN,*) int1,int2,int3
  ELSE IF(PRESENT(int2)) THEN
     READ(MESHIN,*) int1,int2
  ELSE IF(PRESENT(int1)) THEN 
     READ(MESHIN,*) int1
  ELSE 
     CONTINUE ! No action required.
  END IF
END SUBROUTINE READ_LINE

! Extended DBD 20 March to read 8 real parameters.
SUBROUTINE READ_STRING_LINE(string,stemp1,stemp2,stemp3,stemp4,stemp5,&
 real1,real2,real3,real4,real5,real6,real7,real8)
!***************************************************
! This routine was required due to a compiler problem
! under MIPS 7.1 under Irix 6.3, which was unable to 
! correctly compile ADVANCE='NO' on the read statements, 
! forcing a clumsy workaround on subsequent reads.
! This routine  essentially groups much of the workaround. 
! Note that on multi-line cards (eg TE), a blank string may
! be read on subsequent lines using the work-around: this is irrelevant.
!***************************************************

  CHARACTER(15), INTENT(out) :: string

  CHARACTER(5), INTENT(out), OPTIONAL  :: stemp1,stemp2,stemp3,stemp4,stemp5
      
  REAL(SP), INTENT(out), OPTIONAL  :: real1,real2,real3,real4,real5,real6,real7,real8

  string  = ' '  ! This paramater is now redundant, but would require too many changes to 
                 ! code to fix easily.  

  IF(PRESENT(real8)) THEN
     READ(MESHIN,*)stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3,real4,real5,real6,real7,real8
  ELSE IF(PRESENT(real7)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3,real4,real5,real6,real7
  ELSE IF(PRESENT(real6)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3,real4,real5,real6
  ELSE IF(PRESENT(real5)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3,real4,real5
  ELSE IF(PRESENT(real4)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3,real4
  ELSE IF(PRESENT(real3)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2,real3
  ELSE IF(PRESENT(real2)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1,real2
  ELSE IF(PRESENT(real1)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5,real1
  ELSE IF(PRESENT(stemp5)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4,stemp5
  ELSE IF(PRESENT(stemp4)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3,stemp4
  ELSE IF(PRESENT(stemp3)) THEN
     READ(MESHIN,*) stemp1,stemp2,stemp3
  ELSE IF(PRESENT(stemp2)) THEN
     READ(MESHIN,*) stemp1,stemp2
  ELSE IF(PRESENT(stemp1)) THEN 
     READ(MESHIN,*) stemp1
  ELSE 
     CONTINUE ! No action required.
  END IF
END SUBROUTINE READ_STRING_LINE


END SUBROUTINE MESHIN_FEK_OLD
