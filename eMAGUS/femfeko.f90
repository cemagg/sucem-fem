! Status on 26 Aug 04: higher order elements not working properly with scattered field formulation. 
! Interpolation (tested using ONLY prescribed dof) seems to be OK, but 
! overall field solution is poor, esp. near target, with scattered field formulation.

! 27 AUG - STILL UNABLE to find bug. May have to re-test interpolation scheme. 
! Only three places where error can be: in the computation of prescribed dof,
! in the computation of the RHS vector, and in post-processing, 
! assuming that all the data structures for prescribed 2nd type dofs are correct. 
! (Maybe re-check this???!)

! WARNING:  
! An apparent error persists in the scattered field formulation section, when applied to cubical PEC scatterers. 
! DBD Feb 2005.

! Last changed:
! 24 March 03. Work started on time domain version. 
! 14 Dec 02. Additional iterative solver added. DBD. 
! 1 May 02
! Line 823ff: documentation added in routine FEMFEKO_OUTPUT_FE_RESULTS
! to clarify output.
! Line 1020ff: support added for GW fields in non-S-parameter case. 

MODULE FEMFEKO
!!$  USE write_matrix
  USE adaptive_fem
  USE problem_info
  USE boundary_conditions
  USE CBAA_data
  USE coax_feed
  USE eigen_analysis_data
  USE cbaa_sys, ONLY: CBAA_BOUNDARY_SEARCH, &
       CBAA_MAKE_COAX_BVECTOR, &
       CBAA_MATRIX_ALLOCATE, &
       CBAA_COUNT_APERTURE_DOFS, &
       CBAA_PREPROCESSING, &
       COAX_BOUNDARY_SEARCH, &
       CBAA_SYSMAT
  USE gw_sys, ONLY: GW_MATRIX_ALLOCATE, &
       GW_PREPROCESSING, &
       PORT_BOUNDARY_SEARCH, &
       GW_SYSMAT 
  USE feminterface, ONLY: ASSIGN_FACE_AND_EDGE_ORDERS, &
       CONTROL_ADAPTIVE, &
       DIRECT_SOLVE, &
       EIGEN_SYSMAT, &
       FD_SCAT_PRE_X_VEC,FD_SCAT_MATRIX_ALLOCATE, & 
       FD_SCAT_BVECTOR_SCAT, FD_SCAT_BVECTOR_TOT, FD_SCAT_INC_SCAT, & 
       INPUT_ELEMENT_ORDERS, &
       ITER_SOLVE, & 
       MESH_INFO_WRITE, &
       NUMBER_DOF, &
       REORDER,&
       TD_SYSMAT,TD_TIMESTEP,TD_MATRIX_ALLOCATE
   USE far_field_data
   USE frequency_data
   USE geometry
   USE gw_data
   USE inc_field
   USE math_tools, ONLY: TIME_DIFFERENCE
   USE matrix
   USE near_field_data
   USE nrtype
   USE output_error
   USE output_data
   USE probe_feed
   USE quad_tables
   USE unit_numbers
   USE parseinput, ONLY: mesh_format => mesh_type
   USE scattering_analysis_data
   USE material_properties
   IMPLICIT NONE
!*******************************************************************************
! PROGRAM DESCRIPTION
!*******************************************************************************
! This is a 3D Finite Element Method program for the analysis of scattering 
! from strutures contained within the mesh volume. 
! It performs eigenanalysis of unloaded or lossless dielectric loaded 
! metal cavities.
! Lastest extensions include lossy anisotropic materials for perfectly matched
! absorber implementations, sources for scattering analyses of periodic
! surfaces (presently normal incidence only), an interface to FAM
! for post-processing, and the inclusion of perfect electric conductor objects
! within the mesh.
!
! The elements used are tetrahedral elements, H_0(curl) and 
! H_1(curl) (CT/LN and LT/QN respectively, also known as order
! 1/2 and 1 1/2). 
!*******************************************************************************
! AUTHOR
!*******************************************************************************
! This code was originally written by Dr. David B Davidson.
!   Dept Electrical and Electronic Engineering
!   University of Stellenbosch
!   Stellenbosch, South Africa
! during the author's sabbatical visit to Trinity College, Cambridge, England
! Jan - July 1997. 
! It was subsequently extended in Stellenbosch, starting November 1998, by 
! David B Davidson (DBD) and Riana H Geschke in colloboration with Frans JC Meyer
! who developed the pre- and post-processing environments, under 
! contracts with EM  Software and Systems Pty Ltd, 
! Technopark, Stellenbosch, partially supported by THRIP funding from 1997-1999.
! DBD continued to add functionality, including waveguide analysis under contract to EMSS 
! (with SPII support) during the period 2000-2002, 
! and during a sabbatical at TU Delft (2003), new time-domain scattering capabilities were added 
! as academic research. 
! 
! During 2000-2003, major additional functionality was added by Matthys M Botha 
! during his PhD work at the Univ of Stellenbosch, in particular on cavity backed antennas and p-adaptation. 
! Some of this work was done under contract to EMSS, much purely
! as part of his doctoral research. 
!*******************************************************************************
! REVISION HISTORY
!******************************************************************************
!
! This program was originally called EMFEMUS, the last version of which was 
! Version 2.3, released internally only in February 1999. See below for the
! development history of EMFEMUS.
! 
! Version 0.90: Started March 1999 by DBD and RHG. Major new pre-processing
!               facilities added to interface with FEKO and FEMAP. 
!               Development version frozen May 13 1999. 
!               Major extensions implemented: 
!               o Reads in the .pre file created by PreFEM from a FEMAP .neu
!                 file.
!               o Degrees of freedom can be renumbered to minimize bandwidth.
!
! Version 0.91: Started 26 May 1999 by DBD and RHG.
!               Features implemented:
!               o Banded matrix solvers.
!               o Eigenvector computation
!               o Post-processing in FEKO format.
!               Development version frozen 31 Aug 1999.
! 
! Version 0.92: Nov 26 1999. Minor corrections of 0.91 during re-compilation
!               on Irix 6.5.4 using MIPS 7.2.1 f90 compiler. 
!
! Version 0.93: Started Nov 26 99 by DBD. Finished early February 2000. 
!               Intermediate development version, frozen for futher work by
!               FJCM.
!
! Version 0.94: Started early Feb 00 by DBD. Finished 7 March 2000. Mainly
!               internal restructuring, error handling, use of records,  
!               and facility to predict number of spurious eigenmodes.
!
! Version 0.95: First version featuring following new capabilities:
!               o Waveguide analysis.
!               o Cavity-backed aperture antennas. 
!               o Sparse matrix solution for eigenanalysis.
!               Finished 22 June 2000, DBD/RHG/MMB.
!
! Version 0.96: Sparse matrix routines implemented for waveguide and
!               cavity-backed aperture analysis as well. Also includes
!               LT/QN element support for guided wave analysis. 
!               Present release finished 20 July 2000, DBD.
!
! Version 0.97r1: Upgrade of FEKO-like functionality. LT/QN also incorporated
!                 in CBAA analysis. Error handling changed. 31 Jan 2001. MMB/DBD
!               
! Version 0.97r2: GW fully incorporated in the general FEKO-like processing.
!                 <cbaa_pp.f90> removed from makefile.
!                 FJCM updates for DF compiler included, thus the code can
!                 now compile on Win(PC)-Digital Fortran and on SGI-MIPS compiler.
!                 09 Feb 2001. MMB/FJCM
! Version 0.97r3: Release 1 with cubature and new A&V LT/QN elements 
!                 (not implemented for CBAA). 17 Feb 2001. DBD. 
! Version 0.97r4: Release 2 merged with 3. 
!                 19 Feb 2001. DBD. 
! Version 0.97r5: Port definitions generalized - no longer constrained to lie
!                 along z-axis. 7 April 2001. DBD.
! Version 0.97r6: Version 5 with some errors corrected. Apr 2001 DBD.
! Version 0.97r7: Version 5 with LT/QN elements implemented for FMM in CBAA; 
!                 new iterative solvers (QMR and GMRES) and diagonal 
!                 preconditioner. Apr 2001. MMB 
! 
! Version 0.97r8: Version 6 and 7 merged. 24 April 2001. DBD/MMB.
! Version 0.97r9: Corrections to r8, 9 May 2001.  
! Version 0.97r10: Continuing corrections. 22 June 2001. 
!                  FSS removed. MMB 17 Sept 2001. 
! Version 0.98r1: Removed nearly all dead code. Started to restructure 
!                 FEMFEKO such that all routines are located within modules
!                 with their associated data if applicable. Much restructuring 
!                 still to be done before the whole code will conform to this 
!                 idea. MMB Oct 2001.
! Version 0.98r2: Changed FEMFEKO such that all elements do no longer need to
!                 be of the same hierarchal order. MMB 14 Oct 2001.
! Version 0.98r3: Added rigorous, dominant mode, coax port formulation - but
!                 it may still contain bugs. The Gong&Volakis coax model is also
!                 suspected to be buggy. Moved some more routines to
!                 modules <geometry> and <basis_function>. Rewrote B_MAKE_HIERARCHAL
!                 so that it always uses qudrature. Old verson now unused, but
!                 still available as B_MAKE_HIERARCHAL_ANALYTIC. MMB 2002-02-20.
! Version 0.98r4: Moved cubature tables to the <quad_tab> module, and re-wrote
!                 them in the same format as the 2D quadrature routines. 
!                 Added new flags into code to permit polynomial-complete elements
!                 to be incorporated. DBD 2002-02-21.
! Version 0.98r5: Implementation started on polynomial complete 
!                 elements, and also the 2nd order elements proposed
!                 by Webb [IEEE T-AP Aug 99 p 1244ff]. DBD 2002-02-26. 
! Version 0.98r6: 1st order complete elements [Webb99] implemented. DBD. 
!                 S and T matrix routines simplified, only cubature to be 
!                 supported from now on. Support added for up to QT/QN elements.
!                 (Not fully implemented however).
!                 DBD 04 March 2002.
! Version 0.98r7: Additional QT/QN elements [Webb99] fully implemented.
!                 DBD 07 March 2002.
! Version 0.98r8: Corrections to U_MAKE_HIERARCHAL, MAX_ORDER, MIXED_ORDER,
!                 ASSIGN_FACE_AND_EDGE_ORDERS. GW_PREPROCESSING, NUMBER_DOF.
!                 Added explicit, residual estimator (only GW at present). 
!                 Added AD card for full control of the adaptation tasks
!                 from the *.fek file. Additional library file now needed
!                 for sorting in the error estimation process. The library file 
!                 contains the relevant routines from the slatec-library, 
!                 available from www.netlib.org. MMB 2002-03-19.
! Version 0.98r9: Incorporated corrections recieved via email from DBD. Extended
!                 the functionality of the AD card (see SUBROUTINE CONTROL_ADAPTIVE
!                 for all its current capabilities). MMB 2002-04-22.
! Version 0.98r10: Minor extensions by DBD to support GW fields in non-S-parameter case. 
! Version 0.98r11: Added implicit ERM error etimator for CBAA/GW/COAX. Extended 
!                  explicit residual indicator to work for CBAA/GW/COAX. Complete
!                  revamp of CBAA BI routines - they are now essentially basis
!                  function independent, just add VBF definitions in
!                  POINT_EVALUATE_FACE_FUNCTIONS. CBAA BI upgraded to Webb99,
!                  2nd order complete. Temporary, additional real parameter in FC
!                  card. (MMB Note: Will change again soon. 2002-05-24.)
!                  DBD note: MMB developed some further PhD specific code which is not contained
!                  in this development line. 
!
! Version 0.98r12: Experimental mixed potential preconditioner 
!                  added. DBD 2002-12-14.
!
! Version 0.99r1:  First version supporting time-domain analysis. Initial development
!                  started during DBD's sabbatical at TU Delft, Feb-July 2003. 
!                  Finished 9 April 2003. First order ABC's are 
!                  incorporated via the BC card, AB and QU card. 
!                  The formulation used is the scattered field formulation
!                  in Chapter 12, J. Jin, "The FEM in EM", 2nd edn, Wiley 2002.
! Version 0.99r2:  Support added for sparse matrix routines 
!                  using direct solvers (LU factorization)
!                  from the Compaq Extended Math Library  
!                  routines, using skyline storage. DBD. 28 April 2003.
! Version 0.99r3:  Support added for sparse matrix routines 
!                  using iterative solvers using existing FEMFEKO routines 
!                  modified for REAL*8 calculations. Finished 1 May 2003.
!
! Version 0.99r4:  Scattered/total field formulation added (to LT/QN order, but only tested to CT/LN order). 
!                  DBD Finished 6 June 2003. Formulation assumes an interior, inhomogeneous
!                  scatterer and an exterior, homogeneous region (HOMOG_MEDIUM, presently 1). 
!                  Incident field enters via intergral over this interior fictitious boundary 
!                  and also volume integral over inhomogeneous region. Preliminary testing only. 
! Version 0.99r5:  (Approximate) pefectly matched layer added. Uses the approach
!                  of V Mathis ("An Anistropic PML-Absorbing Medium in FETD method for Maxwell's Equations", 
!                  IEEE AP-S 1997. Requires uses of scattered field formulation. 
!                  DBD Finished 12 June 2003. Only preliminary testing performed - with zero sigma's,
!                  code gives same result as one without PML. 
! Version 0.99r6.  Full implementation of PML, based on D.Jiao et al, 
!                  "Time-domain finite-element simulation of three-dimensional scattering and radiation
!                  problems using perfectly matched layers", IEEE T-AP, Feb 2003. 
!                  DBD. Coding finished 16 July 2003. Preliminary testing only completed.
!
! Version 1.00r1.  Research FEMFEKO implementation, with input upgraded to FEMFEKO StandAlone and current WinFEKO 
!                  compliance. (WinFEKO does not support all FEMFEKO cards, however). 
!                  DBD. 15 December 2003. Not extensively tested!
! 
! Version 1.00r2.  Frequency domain scattering analysis option added, analogous to TD option added ver. 0.99, along with
!                  1st order ABC. Total field formulation only; incident field injected at boundary. 
!                  Work finished 18 Dec 2003 DBD.  
! Version 1.00r3.  Continuation of work using 2nd order ABC, with both total (1st order only) and scattered field
!                  (1st and 2nd order) formulations. Only CT/LN order elements fully implemented. Finished 16 Aug 2004 DBD..
! Version 1.00r4.  LT/LN and LT/QN implementation added and tested. Finished 01 Sept 2004 DBD.
! Version 1.00r5.  Unsucessful attempt to implement full 2nd order ABC including surface div term. 30 Sept 2004. DBD.
! Version 1.00r6.  Valid numerical approximation of surface div term implemented. 23 Nov 2004. DBD.
! Version 1.00r7.  Additional scatterers supported. 30 Nov 2004. DBD. Problem with scattering from cube
!                  using scattered fielf formulation noted; problem not yet identified. DBD 27 Jan 2005.
! Version 1.00r8:  Support for non-spherical ABC's and dielectric scatterers (total field formulation only) added.
!                  DBD 28 Jan 2005.
! Version 1.00r9:  Ported to SGI. DBD 8 Feb 2005.
! Version 1.10r1:  Experimental 2nd order curvilinear element implementation, generating mid-point nodes 
!                  using an ad-hoc linear interpolation. DBD 20 July 2005.
! Version 1.10r2:  Further work on curvilinear elements as in 1.10, permitting an option to generate a curvilinear approximation
!                  in the specific case of a spherical scatterer. DBD 27 July 2005.
! Version 1.10r3:  Added treatment for curvilinear elements to post-processing routines.
!
!*******************************************************************************
! DEVELOPMENT HISTORY OF EMFEMUS (predecessor to FEMFEKO)
!******************************************************************************
! 
! Development started Feb 1997.
! Version 1.00: frozen on 4 April 1997.
!               Eigenvalue analysis of unloaded metal cavities.
! Version 1.01: very minor changes; 8 April 1997.             
!
! Version 1.10: Considerable resturcturing to improve:
!               I/O; dynamic memory allocation; material handling.
!               28-30 April 1997.
! 
! Version 2.00: Frozen on 3 June 1997 as Version 2.00.
!               Extension to include periodic surface analysis. 
!               Also includes plane wave sources and lossy anisotropic 
!               materials (diagonal permittivity and permeability 
!               matrices only) to implement a perfectly matched 
!               absorber.
!               (See "Comments and extensions of `A note on the application
!               of edge-elements for modelling three-dimensional 
!               inhomogenously-filled cavities'", DB Davidson
!               IEEE Trans. MTT, Sep 98, pp. 1344-6 for theoretical extensions).
!               Preliminary post-processing done via the generation of 
!               a FAM neutral file (fnf) output.
!               Cavity analysis validated for inhomogenous cavity.
!               
! Version 2.01: June 17-18 1997
!               Following minor facilities added following initial tests:
!               a). Option added to read in data, but not generate
!                   system or solve it. (Useful for testing large datasets).
!                   (Only implemented for FSS analysis).
!               b). Solid metal objects can now be located within mesh.
!                   These objects must be of finite thickness.
!               c). Options added to suppress some of the output listing.
!               d). Timing data added. 
! 
! Version 2.02: June 19 1997. Error in incident field corrected in FSS_SYSMAT.
!
! Version 2.10: June 20 - 22 1997
!               Banded matrix storage with iterative (BiCG solver) implemented.
!
! Version 2.11: June 29 1997. Error in BUILD_FSS_SYSMAT corrected.
!               Conjugate gradient iterative solver option added.
!
! Version 2.12: Oct 2 1998.
!               Some minor errors with "old" USE statements corrected.
!
! Version 2.20: Dec 13 1998. 
!               H_0 (curl) elements (1/2 order, CT/LN) re-implemented 
!               for eigenvalue analysis.
!               See "Hierarchal 2D and 3D Vector Finite 
!               Elements for Electromagnetic Wave Eigenvalue Problems",
!               DB Davidson, RH Hansmann, ACES 15th Annl. Review of Progress
!               in Applied Computational Electromagnetics, March 1999, pp.
!               518--521, for underlying theory.
!               
! Version 2.21: Dec 1998
!               H_1(curl) (1&1/2-th order, LT,QN) hierarchal elements 
!               added. Ibid. Not tested.
!
! Version 2.22: 29 Jan 1998. DBD.
!               Whole code changed to uniformly number nodes from 1 to 4,
!               in accordance with Savage and Peterson. (Apart from one small
!               part of the S&T matrix generation routines specific to 
!               Lee and Mittra, where the old numbering scheme is retained
!               and the necessary conversion done).
!               Eigenvalue part tested with renumbering; works. 
!               CAUTION: FSS part, along with post-processing and 
!               anisotropic media, NOT properly tested following this 
!               modification!
!               Also put through the Salford FORTRAN95 compiler; 
!               some run-time checking done (this option not very reliable).
!               "Frozen" after this port.
!               
! Version 2.23  Work on higher-order elements finished 10 Feb 1999
!               followed by debugging. DBD.
!               
! Version 2.3   Feb 18 1999. DBD. Production code; cleaned-up version of 2.23.
!                             
!*******************************************************************************
! EXECUTION ENVIRONMENT
!*******************************************************************************
! Developed in FORTRAN 90 on Silicon Graphics Power Indigo2, 02 and Octane
! workstations, running IRIX 6.2, 6.3 and 6.5 respectively, 
! and compiled using the MIPS f90 compiler in -n32 bit mode. 
!
! Also ported to the Salford FTN95 compiler. (Obsolete) routine SMAKE will not 
! run with the checking options activated; the error appears to be in the 
! compiler. Several other routines will also not run with checking options
! enabled.
!
! Now primarily developed using the Compaq Visual Fortran suite, with ports to SGi's (R12 000)
! and compiled under MIPS Fortran from time to time. 
!
!*******************************************************************************
! INPUT 
!*******************************************************************************
! Input data uses .pre files created by WinFEKO and FEMAP and an input.dat file.
!
!*******************************************************************************
! OUTPUT 
!*******************************************************************************
! Either: 
! Eigenvalue of  possibly loaded but lossless metal cavities.
! or
! Cavity backed aperture antenna analysis.
! or 
! Guided wave analysis.
! or
! or
! Frequency domain scattered or total field data
! Time domain field data.
!
! The following functionality is obsolete:
! Periodic structure analysis (not currently properly tested).
!
! Data is output is a format suitable for graphical
! post-processing. Development is presently in progress in the FEKO
! environment.An obsolete post-processing interface to FAM 
! has now been removed. 
!
!*******************************************************************************
   CHARACTER(FILENAMELENGTH) :: infilename

! Changed  DBD 20 Feb 02
   CHARACTER(RUN_LABEL_LINE_LENGTH), DIMENSION(RUN_LABEL_NUM_LINES) :: run_label 
   REAL(SP) :: timer2, timer3
   INTEGER(I4B) :: bw, s, k, kedge
   INTEGER(I4B) :: cardcount

   ! Main loop variables:
   INTEGER(I4B) :: jj,cardi,cardj,port_counter          ! counters
   REAL(SP) :: t_pw,p_pw,e_pwE,mag_pwE,phase_pw         ! specifying incident plane wave
   INTEGER(I4B) :: freqcount,                           &
                   analysis_count,MAX_analysis_loop,    &
                   theta_total_points,phi_total_points  ! counters
   INTEGER(I4B), DIMENSION(2) :: excitation_counters    ! (1) <- no. of current probes in this analysis
                                                        ! (2) <- no. of plane waves in this analysis
   INTEGER(I4B), DIMENSION(8) :: time1,time2,time_start,time_finish ! temp variables for timing all events
   REAL(SP) :: admin_time,total_time,reading_time,geoprocess_time,matbuild_time, &
               exbuild_time,solution_time,postcalc_time, &
               FE_time,FF_time,temp_time                      ! timers in seconds

   INTEGER(I4B) debug_temp,debug_counter,ios
   
CONTAINS

 SUBROUTINE run_femfeko
   IMPLICIT NONE

   CALL DATE_AND_TIME (values=time_start)

   infilename='input.dat'

   CALL INIT_TIMERS             ! Initialise timers
   CALL READ_FILES              ! Read input and mesh files
   CALL PROCESS_GEOMETRY        ! General processing of the geometry data
   CALL PROCESS_ANALYSIS_GEO    ! Other analysis-specific, geometry 
                                ! pre-processing
   IF (ADAPTIVE_ANALYSIS.AND. &
     ((ADdata%file_orders.EQ.2).OR.(ADdata%file_orders.EQ.3))) THEN ! read the element orders from a file     
     CALL INPUT_ELEMENT_ORDERS
	 CALL ASSIGN_FACE_AND_EDGE_ORDERS(.FALSE.)
   ELSE
     CALL ASSIGN_FACE_AND_EDGE_ORDERS(.TRUE.) ! newly initialise the element orders
   END IF
   !CALL NUMBER_DOF
   IF (FD_SCAT_ANALYSIS.AND.SCAT_FIELD) CALL NUMBER_PRE_DOF
   CALL MESH_INFO_WRITE
! Added DBD 2 April 2003
   CALL EDGE_FACE_OUTPUT
! End DBD 2 April 2003

!   CALL ALLOCATE_SYS_MATRICES   ! Allocate memory for the system equation

!   CALL material_properties_errcheck

   WRITE (*,'(A)') 'Finished processing geometry data'
   CALL DATE_AND_TIME (values=time2)
   geoprocess_time = geoprocess_time + TIME_DIFFERENCE(time1,time2)
   ! -----------------------------------------------------------------

   ! EIG case treated seperately from CBAA and GW because it is frequency 
   ! independent: (Only ONE analysis should be present.)
   
   IF (REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS) THEN
      CALL output_materials_info
      CALL EIGEN_SYSMAT(timer2,timer3)  
     CALL FEMFEKO_OUTPUT_FE_RESULTS(1) 
   END IF
   ! Ditto TD analysis case.
   IF (TD_ANALYSIS) THEN
       CALL output_materials_info
      IF (PML_PRESENT) CALL PML_ASSIGN_PROPERTIES
      CALL DATE_AND_TIME (values=time1)
      CALL TD_SYSMAT
      CALL DATE_AND_TIME (values=time2)
      matbuild_time = matbuild_time + TIME_DIFFERENCE(time1,time2)
      CALL DATE_AND_TIME (values=time1)
      CALL TD_TIMESTEP
      CALL DATE_AND_TIME (values=time2)
      solution_time = solution_time + TIME_DIFFERENCE(time1,time2)
   END IF

   IF (.NOT.GW_ANALYSIS) THEN 
      ! (Re)set these so that the loops will execute once only.
      num_wg_modes_m = 1
      num_wg_modes_n = 0
   END IF

   ! Main system of loops for CBAA and GW and Freq Domain scattering analysis: 
   ! (will not execute in EIG case because then no FR cards are allowed.)
   ! Find MAX_analysis_loop:
   

   PRINT*, "Before frequency loop"
   FREQUENCY_LOOP: DO freqcount = 0,FRdata%nfreq-1
      ! Calculate frequency,k0 and lambda0:
      PRINT*, "In Frequency loop"
      SELECT CASE(FRdata%freqf)
      CASE(0)
         frequency =   FRdata%freq0                   &
              + REAL(freqcount)*FRdata%dfreq
      CASE(1)
         frequency =   FRdata%freq0                   &
              * ( FRdata%dfreq**freqcount )
      END SELECT
      lambda0 = c_0/frequency
      k0 = 2.0*PI/lambda0

      ! All analysis' using this frequency&material (i.e system matrix) 
      ! will be processed:
      MAX_analysis_loop = NUM_analysis
      
!!! This next do loop loops over all the analyses to be done. In the past 
!!! multiple FR cards were supported, but now we can just start from 1 NM
      DO analysis_count=1, MAX_analysis_loop
         ! Loop over the waveguide modes - or execute loop once if not 
         ! guided wave analysis. 
         CALL MODE_LOOPS
      END DO

   END DO FREQUENCY_LOOP
   ! -----------------------------------------------------------

   IF (COUPLED_TD_ANALYSIS) THEN
!      call coupled_td_dostuff

      RETURN
!!$      STOP
      
   END IF
   ! Calculate element adaptations: (Added: 2002-03-11. MMB.)
   IF (ADAPTIVE_ANALYSIS) THEN

     ! Correct the the current GW mode values (the loops set them to one more
     ! than they should be):
     mode_m = num_wg_modes_m
     mode_n = num_wg_modes_n

     ! Do error stuff:
     CALL CONTROL_ADAPTIVE 

   END IF
   ! -----------------------------------------------------------

   ! Calculate the total time and admin time:
   CALL DATE_AND_TIME (values=time_finish)
   total_time = TIME_DIFFERENCE(time_start,time_finish)
   admin_time = total_time - (reading_time+geoprocess_time+matbuild_time+exbuild_time &
                              +solution_time+postcalc_time+FE_time+FF_time)

   ! Write breakdown of the runtime to the output file:
   WRITE(FILEOUT,'(//,20X,A)') 'PROGRAM TIMINGS IN SECONDS'
   WRITE(FILEOUT,'(/,1X,A,2X,F12.3)') 'Reading input file(s):                        ', reading_time
   WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Processing the geometry:                      ', geoprocess_time
   IF(TD_ANALYSIS) THEN
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Building and factoring (if relevant)'
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'the system matrix:                            ', matbuild_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Time stepping (overall)             :         ', solution_time
   ELSE
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Building the system matrix:                   ', matbuild_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Building the excitation vector:               ', exbuild_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Solving the matrix equation:                  ', solution_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Calculation of impedances/powers/losses:      ', postcalc_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Calculation of near fields:                   ', FE_time
     WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Calculation of far fields:                    ', FF_time
   END IF
   WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Other administration:                         ', admin_time
   WRITE(FILEOUT,'(1X,A)')            '------------------------------------------------------------'
   WRITE(FILEOUT,'(1X,A,2X,F12.3)')   'Total elapsed time:                           ', total_time

   CLOSE (UNIT=INFILE) 
   
   CALL CLEANUP_MATRIX_MEMORY
   CALL CLEANUP_MODEL_MEMORY
   CLOSE (FILEOUT) ! Opened in subroutine MESHIN_FEK

 END SUBROUTINE run_femfeko

!*******************************************************************************

! Changes to subroutine Mar 16, 20  & 23 DBD.

  SUBROUTINE A0_LOOPS
    A0_THETA_LOOP: DO cardi = 0,theta_total_points-1
       A0_PHI_LOOP: DO cardj = 0,phi_total_points-1

          IF (GW_ANALYSIS.AND.GW_COMPUTE_S_PARAMS) THEN 
             CALL S_PARAMS_CALC ! Calculate S parameters
          ELSE 
             
             ! DBD change 30 Aug 2004. Test purposes only. 
             CALL FEMFEKO_EXCITATION_VECTOR(analysis_count,cardi,cardj)
             IF (.NOT. EXECUTE) THEN ! This flag skips matrix solution 
                WRITE(FILEOUT,'(A)') 'EXECUTE FLAG OFF, no matrix solution'
                EXIT
             END IF
             RUN_SOLVER: IF(.NOT.TEST_INTERPOLATE_FIELD) THEN 
                ! Perhaps this test should be combined with above EXECUTE test
                ! write(fileout,*) b_vec_c 
                ! Solve the matrix equation:
                SELECT CASE (SOLVER_TYPE)
                CASE (0)
                   ! Extended DBD 16 Dec
                   IF (CBAA_ANALYSIS.OR.GW_ANALYSIS.OR.FD_SCAT_ANALYSIS) CALL DIRECT_SOLVE(1,temp_time) ! back substitution
                   !write (fileout,*)'Solution vector'
                   !write (fileout,*) x_vec_c
                   ! End extended DBD 16 Dec
!!! DBD change 06 Dec: was
!!!                   CASE (1:4)
!!! now 
                CASE (1:5)
                   ! Extended DBD 16 Dec
!!$                   CALL write_sparse_c(Asparse_c, col_ind, row_ind)
!!$                   CALL write_vector_c(b_vec_c, "out.b_vec_c")
                   IF (CBAA_ANALYSIS.OR.GW_ANALYSIS.OR.FD_SCAT_ANALYSIS) CALL ITER_SOLVE(temp_time)
!!$                   CALL write_vector_c(x_vec_c, "out.x_vec_c")
!!$                   STOP
                   ! End extended DBD 16 Dec
                END SELECT

                !Write(fileout,'(///)')
                !Write(fileout,*) 'x_vec_c', x_vec_c
             END IF RUN_SOLVER
             solution_time = solution_time + temp_time
             ! Post processing:
             CALL FEMFEKO_POSTPROCESS(analysis_count)
             ! Write results to the *.out file:
             CALL FEMFEKO_OUTPUT_FE_RESULTS(analysis_count)
             CALL FEMFEKO_OUTPUT_FF_RESULTS(analysis_count)
          END IF 
       END DO A0_PHI_LOOP
    END DO A0_THETA_LOOP
  END SUBROUTINE A0_LOOPS

  SUBROUTINE ALLOCATE_SYS_MATRICES
    ! Allocate memory for the system equation:
    IF (CBAA_ANALYSIS) THEN 
       CALL CBAA_COUNT_APERTURE_DOFS
       PRINT*, 'Allocating memory for CBAA'
       CALL CBAA_MATRIX_ALLOCATE
    END IF
    IF (GW_ANALYSIS)   CALL GW_MATRIX_ALLOCATE
    IF (FD_SCAT_ANALYSIS)   CALL FD_SCAT_MATRIX_ALLOCATE
    ! DBD addition 25 March 2003   
    IF (TD_ANALYSIS)   CALL TD_MATRIX_ALLOCATE
    ! End DBD addition 25 March 2003   
  END SUBROUTINE ALLOCATE_SYS_MATRICES


SUBROUTINE CLEANUP_MATRIX_MEMORY
!*******************************************************************************
! This routine deallocates the remainder of the dynamically allocated data 
! structures. Some may not have been allocated, hence the tests.
! Last changed 08 Aug 2004 DBD.
!*******************************************************************************

  IF (CBAA_ANALYSIS.OR.GW_ANALYSIS) THEN
    IF (.NOT.SPARSE) THEN
      IF (ALLOCATED(A_mat_c)) DEALLOCATE (A_mat_c) ! Storage for the system matrix
      IF (ALLOCATED(ipiv))    DEALLOCATE (ipiv)
      IF (TD_ANALYSIS) THEN
        IF (ALLOCATED(A_mat)) DEALLOCATE (A_mat) 
        IF (ALLOCATED(B_mat)) DEALLOCATE (B_mat) 
        IF (ALLOCATED(C_mat)) DEALLOCATE (C_mat)
	  END IF 
    ELSE
	  IF (TD_ANALYSIS) THEN
        IF (ALLOCATED(Asparse)) DEALLOCATE (Asparse) 
        IF (ALLOCATED(Bsparse)) DEALLOCATE (Bsparse) 
        IF (ALLOCATED(Csparse)) DEALLOCATE (Csparse) 
        IF (ALLOCATED(Asparse_skyline )) DEALLOCATE (Asparse_skyline ) 
        IF (ALLOCATED(IAUdiag)) DEALLOCATE (IAUdiag) 
        IF (ALLOCATED(iwrk)) DEALLOCATE (iwrk) 
        IF (ALLOCATED(rwrk)) DEALLOCATE (rwrk) 
	  ELSE 
        IF (ALLOCATED(Asparse_c)) DEALLOCATE (Asparse_c)
        IF (ALLOCATED(col_ind))   DEALLOCATE (col_ind)
        IF (ALLOCATED(row_ind))   DEALLOCATE (row_ind)
      END IF
    END IF  
    IF (ALLOCATED(x_vec_c)) DEALLOCATE (x_vec_c)
    IF (ALLOCATED(x_pre_vec_c)) DEALLOCATE (x_pre_vec_c)
    IF (ALLOCATED(b_vec_c)) DEALLOCATE (b_vec_c)
  END IF

  IF (CBAA_ANALYSIS) THEN
    IF (SPARSE) THEN
      IF (.NOT.CBAA_FMM_storage) THEN
        IF (ALLOCATED(CBAA_BE_mat)) DEALLOCATE(CBAA_BE_mat)
      ELSE
        IF (ALLOCATED(edge_grnum))         DEALLOCATE (edge_grnum)
        IF (ALLOCATED(gr_matindex))        DEALLOCATE (gr_matindex)
        IF (ALLOCATED(gredge_colind))      DEALLOCATE (gredge_colind)
        IF (ALLOCATED(gredge_rowind))      DEALLOCATE (gredge_rowind)
        IF (ALLOCATED(CBAA_BE_colind))     DEALLOCATE (CBAA_BE_colind)
        IF (ALLOCATED(CBAA_BE_rowind))     DEALLOCATE (CBAA_BE_rowind)
        IF (ALLOCATED(k_dirs))             DEALLOCATE (k_dirs)
        IF (ALLOCATED(CBAA_BE_val))        DEALLOCATE (CBAA_BE_val)
        IF (ALLOCATED(CBAA_groupmat))      DEALLOCATE (CBAA_groupmat)
        IF (ALLOCATED(CBAA_elemvec_term1)) DEALLOCATE (CBAA_elemvec_term1)
        IF (ALLOCATED(CBAA_elemvec_term2)) DEALLOCATE (CBAA_elemvec_term2)
        IF (CBAA_FMM_debug) THEN
		  IF (ALLOCATED(temp_BE_mat)) DEALLOCATE(temp_BE_mat)
        END IF
	  END IF
    END IF
  END IF   

  IF (REAL_EIGEN_ANALYSIS.OR.CMPLX_EIGEN_ANALYSIS) THEN
    IF (ALLOCATED(eigenvectors)) DEALLOCATE(eigenvectors)
    IF (ALLOCATED(eigenvalues))  DEALLOCATE(eigenvalues)
  END IF

  IF (ALLOCATED(renumbered_e1)) DEALLOCATE (renumbered_e1)
  IF (ALLOCATED(renumbered_e2)) DEALLOCATE (renumbered_e2)   
  IF (ALLOCATED(renumbered_e3)) DEALLOCATE (renumbered_e3)   
  IF (ALLOCATED(renumbered_f1)) DEALLOCATE (renumbered_f1)   
  IF (ALLOCATED(renumbered_f2)) DEALLOCATE (renumbered_f2) 
  IF (ALLOCATED(renumbered_f3)) DEALLOCATE (renumbered_f3) 
  IF (ALLOCATED(dof_type))      DEALLOCATE (dof_type)
 
  IF (ALLOCATED(renumbered_pre_e1)) DEALLOCATE (renumbered_pre_e1)
  IF (ALLOCATED(renumbered_pre_e2)) DEALLOCATE (renumbered_pre_e2)   
  IF (ALLOCATED(renumbered_pre_e3)) DEALLOCATE (renumbered_pre_e3)   
  IF (ALLOCATED(renumbered_pre_f1)) DEALLOCATE (renumbered_pre_f1)   
  IF (ALLOCATED(renumbered_pre_f2)) DEALLOCATE (renumbered_pre_f2) 
  IF (ALLOCATED(renumbered_pre_f3)) DEALLOCATE (renumbered_pre_f3) 
  IF (ALLOCATED(pre_dof_type))      DEALLOCATE (pre_dof_type)
 
END SUBROUTINE CLEANUP_MATRIX_MEMORY
!*******************************************************************************


SUBROUTINE CLEANUP_MODEL_MEMORY
  USE material_properties
!*******************************************************************************
! This routine deallocates the remainder of the dynamically allocated data 
! structures. Some may not have been allocated, hence the tests.
!*******************************************************************************

  IF (CBAA_ANALYSIS) THEN
    DEALLOCATE (ap_elnumbers)
    DEALLOCATE (which_local_edges)
    DEALLOCATE (which_local_face)
  END IF   

  ! Deallocate vertices, edges, faces and elements data:
  CALL DEALLOCATE_GEOMETRY(dealloc_vertices=.TRUE., dealloc_edges=.TRUE., &
                           dealloc_faces=.TRUE., dealloc_elements=.TRUE., &
                           dealloc_connect=.TRUE.)  
  
  ! Deallocate triangular and tetrahedral quadrature rules:
  CALL QUAD_LINE_RULES_CLEAN ! Added DBD 04 June 2003
  CALL QUAD_TRI_RULES_CLEAN
  CALL CUBE_TET_RULES_CLEAN

  CALL material_properties_deallocate

  DEALLOCATE (BCs)
  DEALLOCATE (ports)
  IF(ALLOCATED(FEdata))  DEALLOCATE (FEdata) ! FE card parameters
  IF(ALLOCATED(FFdata))  DEALLOCATE (FFdata) ! FF card parameters
  IF(ALLOCATED(A0data))  DEALLOCATE (A0data) ! A0 card parameters
  IF(ALLOCATED(A8data))  DEALLOCATE (A8data) ! A8 card parameters
  IF(ALLOCATED(A9data))  DEALLOCATE (A9data) ! A9 card parameters
 
END SUBROUTINE CLEANUP_MODEL_MEMORY
!*******************************************************************************


  SUBROUTINE FEMFEKO_ASSIGN_PORTS(this_analysis,port_counter)
!*******************************************************************************
! Assign the relevant values to the port excitations in the array <ports>, for
! the current analysis.
! Start DBD addition 23 March
! If the waveguide is recognized as a standard size, this is printed in the
! output file.
! Standard sizes from:
!   D M Pozar, "Microwave Engineering", 2nd edn, Wiley 1998, 
!   Appendix I, p.706
! End DBD addition 23 March
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: this_analysis       ! analysis number to be processed
    INTEGER(I4B), INTENT(IN), OPTIONAL ::  port_counter ! Used for full S 
                                                          ! parameter excitation.  

    LOGICAL (LGT),  SAVE ::  output_port_info = .TRUE.
                                     ! triggers or suppresses writing of port info.
                                     ! Always true on first pass, may be set
                                     ! false in routine and SAVEd for next pass.
    INTEGER(I4B), SAVE:: last_mode_m, last_mode_n
    INTEGER(I4B) :: i_port           ! counters
! DBD change 2 Apr 2001:
    CHARACTER(LEN=3) :: port_mode    ! used for writing out the mode type
    REAL(SP) :: a,b                  ! temp variables to identify the WG type

! Change of port definition DBD: 
    ! First initialize ALL excitations to zero, TE:
    ports%type = 1
    ports%excitation = (0.0,0.0)
! DBD addition 2 Apr 2001:
    IF (output_port_info) THEN ! First pass, port data needed. 
      last_mode_m = mode_m
      last_mode_n = mode_n
    END IF
    IF (last_mode_m.NE.mode_m.OR.last_mode_n.NE.mode_n) THEN
      ! Mode has changed, output port info again.
      output_port_info = .TRUE.
    END IF
! End DBD addition 2 Apr 2001
    ! Cycle through all A9 cards, using those corresponding to this aor earlier analysis:
    ! (If an excitation was redefined in this or previous analysis, it will simply be overwritten
    ! with its last definition falling within the group of considered A9 cards.)
    DO i_port = 1,NUM_A9_cards
      IF (A9data(i_port)%analysis_no.LE.this_analysis) THEN ! this is a relevant card
        ports(A9data(i_port)%port_num)%type       = A9data(i_port)%port_type
! DBD change 16 Mar 01
       IF (GW_COMPUTE_S_PARAMS) THEN
         IF(i_port.EQ.port_counter) THEN 
           ports(A9data(i_port)%port_num)%excitation = CMPLX(1.0_SP,0.0_SP)
         END IF 
       ELSE 
         ports(A9data(i_port)%port_num)%excitation = &
           CMPLX(A9data(i_port)%re_excitation,A9data(i_port)%im_excitation)
       END IF
! End DBD change 16 Mar
      END IF
    END DO

    ! Write the excitations to the output file, if required:
    IF (output_port_info) THEN  ! DBD change 19 Mar 01
      PORT_LOOP: DO i_port = 1,NUM_ports
    WRITE (FILEOUT,'(//,20X,A)') 'PORT INFORMATION'
    WRITE (FILEOUT,'(/,1X,A,I4)')  'Number of port:                   N = ', i_port
    WRITE (FILEOUT,'(1X,A,E12.5)') 'Frequency in Hz:               FREQ = ', frequency
    WRITE (FILEOUT,'(1X,A,E12.5)') 'Port dimension a in m:            a = ', ports(i_port)%a
    WRITE (FILEOUT,'(1X,A,E12.5)') 'Port dimension b in m:            b = ', ports(i_port)%b
! DBD change 20 Mar 01
    WRITE (FILEOUT,'(1X,A,E12.5)') 'Normal vec. away from mesh:     n_x = ', ports(i_port)%normal(1)
    WRITE (FILEOUT,'(1X,A,E12.5)') '                                n_y = ', ports(i_port)%normal(2)
    WRITE (FILEOUT,'(1X,A,E12.5)') '                                n_z = ', ports(i_port)%normal(3)
    WRITE (FILEOUT,'(1X,A,E12.5)') 'Tangential vector on port:      t_x = ', ports(i_port)%tangent(1)
    WRITE (FILEOUT,'(1X,A,E12.5)') '(This defines the positive      t_y = ', ports(i_port)%tangent(2)
    WRITE (FILEOUT,'(1X,A,E12.5)') 'sense for the eigenmode.)       t_z = ', ports(i_port)%tangent(3)
! End DBD change 20 Mar 01
    SELECT CASE (ports(i_port)%type)
    CASE (1)
          port_mode = 'TE'
    END SELECT
    WRITE (FILEOUT,'(1X,A,A,I1,I1)')     'Excited mode:                  MODE = ',& 
               port_mode,mode_m,mode_n
! DBD change 16 Mar 01 and 9 Apr

        IF (.NOT.GW_COMPUTE_S_PARAMS.OR.DEBUG_GW_PRE) THEN
          WRITE (FILEOUT,'(1X,A,F9.3,A7,F9.3,A2)')                                     &
                                       'Relative excitation:             E0 = ',     &
                                       REAL(ports(i_port)%excitation), ' + i*( ',   &
                                       AIMAG(ports(i_port)%excitation), ' )'
    END IF
! End DBD change 16 Mar and 9 Apr
    a = ports(i_port)%a
    b = ports(i_port)%b
    IF (ABS(a-16.51E-2).LT.EPS.AND.ABS(b-8.255E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = L-band WR-650'
    ELSE IF (ABS(a-7.214E-2).LT.EPS.AND.ABS(b-3.404E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = S-band WR-284'
    ELSE IF (ABS(a-3.485E-2).LT.EPS.AND.ABS(b-1.580E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = C(J)-band WR-137'
    ELSE IF (ABS(a-2.850E-2).LT.EPS.AND.ABS(b-1.262E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = W(H)-band WR-112'
    ELSE IF (ABS(a-2.286E-2).LT.EPS.AND.ABS(b-1.016E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = X-band WR-90'
    ELSE IF (ABS(a-1.580E-2).LT.EPS.AND.ABS(b-0.790E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = Ku(P)-band WR-62'
    ELSE IF (ABS(a-1.07E-2).LT.EPS.AND.ABS(b-0.43E-2).LT.EPS) THEN
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = K-band WR-42'
    ELSE 
          WRITE(FILEOUT,'(1X,A)')      'Standard waveguide type:       TYPE = unrecognized'
    END IF
      END DO PORT_LOOP ! DBD change 19 Mar 01
    END IF 

! DBD change 19 Mar 01 and 2 Apr and 9 Apr
    IF (GW_COMPUTE_S_PARAMS.AND..NOT.DEBUG_GW_PRE) THEN
      output_port_info = .FALSE. ! Don't output data again.
    END IF 
! DBD change 19 Mar 01 and 2 Apr and 9 Apr

  END SUBROUTINE FEMFEKO_ASSIGN_PORTS
!*******************************************************************************


  SUBROUTINE FEMFEKO_EXCITATION_VECTOR(this_analysis,theta_num,phi_num,port_counter)
    USE coax_feed
    USE gw_sys, ONLY: GW_MAKE_PORT_BVECTOR
    USE cbaa_sys, ONLY: CBAA_MAKE_COAX_A_AND_B, &
         CBAA_MAKE_CURRENTPROBE_BVECTOR,               &
         CBAA_MAKE_PLANEWAVE_BVECTOR
!*******************************************************************************
! Sets up the excitation vector in the CBAA/GW/FD scattering cases.
! Last changed: DBD 07 Aug 2004
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: this_analysis           ! analysis number to be processed
    INTEGER(I4B), INTENT(IN) :: theta_num,phi_num       ! A0 direction counters
    INTEGER(I4B), INTENT(IN), OPTIONAL :: port_counter  ! WG port counter
    INTEGER(I4B) :: ex_count

    CALL DATE_AND_TIME (values=time1)

    excitation_counters = 0 ! array assignment; initialization, global variable
    b_vec_c = (0.0,0.0) ! initialise vector

    PRINT*, 'In FEMFEKO_EXCITATION_VECTOR' ! DEBUGGER

    IF (CBAA_ANALYSIS) THEN
       ! Start MMB extension 24 Apr 01
       ! Coaxial feed contribution: (not fully integrated yet for commercial FEMFEKO)
       PRINT*, 'CBAA_ANALYSIS in FEMFEKO_EXCITATION_VECTOR' ! DEBUGGER
       IF (NUM_AX_cards.GT.0) THEN
          IF (COAX_WHITNEY) THEN
             CALL CBAA_MAKE_COAX_A_AND_B(1)
          ELSE     
             CALL CBAA_MAKE_COAX_BVECTOR(k0)
          END IF
          excitation_counters(1) = 1 ! temporary workaround to make femfeko calculate the directivity
       END IF

       ! A8 card contributions:
       PRINT*, 'NUM_A8_cards: ', NUM_A8_cards
       DO ex_count = 1,NUM_A8_cards
          IF (this_analysis.EQ.A8data(ex_count)%analysis_no)  THEN ! this card is part of the analysis
             CALL CBAA_MAKE_CURRENTPROBE_BVECTOR(A8data(ex_count)%prx0,A8data(ex_count)%pry0,A8data(ex_count)%prz0, &
                  A8data(ex_count)%prdir,A8data(ex_count)%prlen,A8data(ex_count)%prrad,A8data(ex_count)%prabs,          &
                  A8data(ex_count)%prphase*PI/180.0)
             excitation_counters(1) = excitation_counters(1) + 1
          END IF
       END DO

       ! A0 card contributions:
       DO ex_count = 1,NUM_A0_cards
          IF (this_analysis.EQ.A0data(ex_count)%analysis_no)  THEN ! this card is part of the analysis
             excitation_counters(2) = excitation_counters(2) + 1
             ! These variables will possibly be used for RCS calculation porposes:
             t_pw     = (A0data(ex_count)%EIR3 + REAL(theta_num)*A0data(ex_count)%DThei)*PI/180.0
             p_pw     = (A0data(ex_count)%EIR4 + REAL(phi_num)  *A0data(ex_count)%DPhii)*PI/180.0
             e_pwE    = A0data(ex_count)%EIR5*PI/180.0
             mag_pwE  = A0data(ex_count)%EIR1
             phase_pw = A0data(ex_count)%EIR2*PI/180.0
             CALL CBAA_MAKE_PLANEWAVE_BVECTOR(t_pw,p_pw,e_pwE,mag_pwE,phase_pw)
          END IF
       END DO
    END IF

    IF (GW_ANALYSIS) THEN
       IF(GW_COMPUTE_S_PARAMS) THEN
          CALL FEMFEKO_ASSIGN_PORTS(this_analysis,port_counter)
       ELSE
          CALL FEMFEKO_ASSIGN_PORTS(this_analysis)
       END IF
       CALL GW_MAKE_PORT_BVECTOR
    END IF

    IF (FD_SCAT_ANALYSIS) THEN
       ! Currently, only one incident field with one set of angles supported.
       IF(SCAT_FIELD) THEN
          CALL FD_SCAT_PRE_X_VEC ! NB! Call first to compute prescribed degrees of freedom - also sets up incident field parameters.
          IF(.NOT.TEST_INTERPOLATE_FIELD) CALL FD_SCAT_BVECTOR_SCAT(k0)
          ! CALL FD_SCAT_INC_SCAT(frequency) ! Incorrect in the context of PEC.
       ELSE 
          CALL FD_SCAT_BVECTOR_TOT(k0)
       END IF
    END IF


    ! Check that there are sources present, this must be a warning:
    IF (SUM(excitation_counters).EQ.0) THEN
       !print *,'No sources'
    END IF

    CALL DATE_AND_TIME (values=time2)
    exbuild_time = exbuild_time + TIME_DIFFERENCE(time1,time2)

  END SUBROUTINE FEMFEKO_EXCITATION_VECTOR
!*******************************************************************************


  SUBROUTINE FEMFEKO_OUTPUT_FE_RESULTS(this_analysis)
    USE feminterface, ONLY: OUTPUT_FE_RESULTS
!*******************************************************************************
! Outputs near field results for the specified analysis.
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: this_analysis      ! analysis number to be processed

    INTEGER(I4B) :: fec_count,i_mode                    ! Counters
    REAL(SP) previous_eigen_error, current_eigen_error  ! variables used for obtaining first 
                                                        ! eigenmode to be processed

    CALL DATE_AND_TIME (values=time1)

    ! Process all cards belonging to this analysis:    
    DO fec_count = 1,NUM_FE_cards
      IF (this_analysis.EQ.FEdata(fec_count)%analysis_no)  THEN ! this card is part of the analysis

        ! Confirm reading the card to be processed:
        WRITE (FILEOUT,'(/,1X,A)') 'read from memory:'
        WRITE (FILEOUT,'(1X,A,5X,5(1X,I5),6(1X,E12.5))') 'FE', &
             FEdata(fec_count)%feltyp, FEdata(fec_count)%n_x, &
             FEdata(fec_count)%n_y, FEdata(fec_count)%n_z, &
             FEdata(fec_count)%felkor, FEdata(fec_count)%x0, &
             FEdata(fec_count)%y0, FEdata(fec_count)%z0, &
             FEdata(fec_count)%delta_x, FEdata(fec_count)%delta_y, &
             FEdata(fec_count)%delta_z

! Changed DBD 02 Apr 2003
        ! The CBAA and GW cases:
        if_CBAA_or_GW_or_TD_or_FD: IF (CBAA_ANALYSIS.OR.GW_ANALYSIS.OR.TD_ANALYSIS.OR.FD_SCAT_ANALYSIS) THEN
          ! Write the additional information header (TEMP for visualization purposes):
          ! Note that the names are incorrect in terms of what is 
		  ! being written, this is due to incomplete WinFEKO support
	      ! when this was written. The references to EIGENMODE etc
		  ! should be corrected when WinFEKO incorporates full 
		  ! FEMFEKO support. 
		  WRITE (FILEOUT,'(//20X,A,5X,I6)') 'EIGENMODE NUMBER: ', 0
          WRITE (FILEOUT,'(/20X,A,5X,G14.5/)') 'EIGENVALUE =', 0.0
          WRITE (FILEOUT,'(20X,A/)') 'NORMALIZED EIGENMODE RESULTS:'
          CALL OUTPUT_FE_RESULTS(FEdata(fec_count)%feltyp,FEdata(fec_count)%n_x,       &
                                 FEdata(fec_count)%n_y, FEdata(fec_count)%n_z,         &
                                 FEdata(fec_count)%felkor, FEdata(fec_count)%x0,       &
                                 FEdata(fec_count)%y0, FEdata(fec_count)%z0,           &
                                 FEdata(fec_count)%delta_x, FEdata(fec_count)%delta_y, &
                                 FEdata(fec_count)%delta_z)
        END IF if_CBAA_or_GW_or_TD_or_FD
! End changed DBD 02 Apr 2003

        ! The EIG case:
        if_EIG: IF (REAL_EIGEN_ANALYSIS) THEN
          ! If an approximation of the eigenmode was provided, search for the closest
          ! to use as the first eigenmode:
          IF (first_eigenmode.EQ.0) THEN
            previous_eigen_error = SQRT(ABS(eigenvalues(1)-first_eigen_approx**2))
            first_eigenmode = 1 ! Will always compute at least one mode.
            DO i_mode = 2,dof
              current_eigen_error = SQRT(ABS(eigenvalues(i_mode)-first_eigen_approx**2))
              IF (current_eigen_error.LT.previous_eigen_error) THEN
                first_eigenmode = i_mode
                previous_eigen_error = current_eigen_error
              END IF
            END DO
          END IF     

          ! Write out the PREDICT_SPURIOUS_EIGENMODES results if applicable:
          IF (PREDICT_SPURIOUS_EIGENMODES) THEN
            WRITE (FILEOUT,'(10X,A,1X,I6)') 'First eigenmode computed (predicted): ',first_eigenmode
          ELSE
            WRITE (FILEOUT,'(10X,A,1X,I6)') 'First eigenmode computed (requested): ',first_eigenmode
          END IF
          WRITE (FILEOUT,'(10X,A,1X,I6)') 'Number of eigenmodes computed: ',num_eigenmodes
  
          EIGENMODE_LOOP: DO i_mode = first_eigenmode,first_eigenmode+num_eigenmodes-1

            ! Check if this mode has been calculated:
            IF (i_mode.GT.dof) THEN
              WRITE (FILEOUT,'(//20X,A,5X,I6,/1X,A,5X,I6,A)')                      & 
                'Eigenmode number ',i_mode, ' cannot be computed since ',i_mode,   & 
                ' exceeds the number of degrees of freedom.'
              WRITE (FILEOUT,'(//20X,A)') 'NO FURTHER EIGENMODES COMPUTED.'
              EXIT EIGENMODE_LOOP ! No more eigenmodes can be computed.
            END IF

            ! Write the near field data block header:
            WRITE (FILEOUT,'(//,20X,A)') 'NORMALIZED EIGENMODE RESULT:'
            WRITE (FILEOUT,'(20X,A,I14)')   'Eigenmode number = ',i_mode
            WRITE (FILEOUT,'(20X,A,G14.5)') 'Eigenvalue       = ',SQRT(eigenvalues(i_mode))

            CALL OUTPUT_FE_RESULTS(FEdata(fec_count)%feltyp,FEdata(fec_count)%n_x,       &
                                   FEdata(fec_count)%n_y, FEdata(fec_count)%n_z,         &
                                   FEdata(fec_count)%felkor, FEdata(fec_count)%x0,       &
                                   FEdata(fec_count)%y0, FEdata(fec_count)%z0,           &
                                   FEdata(fec_count)%delta_x, FEdata(fec_count)%delta_y, &
                                   FEdata(fec_count)%delta_z,i_mode)
          END DO EIGENMODE_LOOP
        END IF if_EIG

      END IF
    END DO

    CALL DATE_AND_TIME (values=time2)
    FE_time = FE_time + TIME_DIFFERENCE(time1,time2)

  END SUBROUTINE FEMFEKO_OUTPUT_FE_RESULTS
!*******************************************************************************

  SUBROUTINE FEMFEKO_OUTPUT_FF_RESULTS(this_analysis)
    USE feminterface, ONLY: OUTPUT_FF_RESULTS
!*******************************************************************************
! Writes requested calculation results to the output file.
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: this_analysis      ! analysis number to be processed

    INTEGER(I4B) :: jj      ! counter

    CALL DATE_AND_TIME (values=time1)

    ! Process FF cards belonging to this analysis:
    FF_loop: DO jj = 1,NUM_FF_cards
      IF (analysis_count.EQ.FFdata(jj)%analysis_no)  THEN ! this card is part of the analysis

         IF (.NOT.CBAA_ANALYSIS) THEN
           print *,'WARNING: FF card ignored in this analysis type'
           CYCLE FF_loop
         END IF

         ! write the card header:
         WRITE (FILEOUT,'(/,1X,A)') 'read from memory:'
         WRITE (FILEOUT,'(1X,A,5X,4(I5),4(E12.5))') 'FF', FFdata(jj)%FFreq, &
                                        FFdata(jj)%NTheta, FFdata(jj)%NPhi, &
                                        FFdata(jj)%rige, FFdata(jj)%Theta0, &
                                        FFdata(jj)%Phi0, FFdata(jj)%DTheta, &
                                        FFdata(jj)%DPhi
                                     
         SELECT CASE (FFdata(jj)%ffreq)

         CASE (0) ! No calculation
           CONTINUE ! No FF results calculated

         CASE (1) ! Use the directions specified by the FF card
           IF ((excitation_counters(1).GT.0).AND.(excitation_counters(2).EQ.0)) THEN ! only internal sources
             SELECT CASE(FFdata(jj)%rige)
             CASE (0) ! FF 1 with directivity calculation
               CALL OUTPUT_FF_RESULTS(FFdata(jj)%ntheta,FFdata(jj)%nphi,      &
                      FFdata(jj)%theta0*PI/180.0,FFdata(jj)%phi0*PI/180.0,    &
                      FFdata(jj)%dtheta*PI/180.0,FFdata(jj)%dphi*PI/180.0,2)
             CASE (1) ! FF 1 with gain calculation
               CALL OUTPUT_FF_RESULTS(FFdata(jj)%ntheta,FFdata(jj)%nphi,      &
                      FFdata(jj)%theta0*PI/180.0,FFdata(jj)%phi0*PI/180.0,    &
                      FFdata(jj)%dtheta*PI/180.0,FFdata(jj)%dphi*PI/180.0,3)
             END SELECT
          ELSE IF ((excitation_counters(1).EQ.0).AND.(excitation_counters(2).EQ.1)) THEN ! only one external source
            ! because there is only one plane wave excitation, the correct values are still available (t_pw,p_pw etc.)
            ! FF 1 with bistatic RCS calculation
            CALL OUTPUT_FF_RESULTS(FFdata(jj)%ntheta,FFdata(jj)%nphi,             &
                   FFdata(jj)%theta0*PI/180.0,FFdata(jj)%phi0*PI/180.0,           &
                   FFdata(jj)%dtheta*PI/180.0,FFdata(jj)%dphi*PI/180.0,1,mag_pwE)
          ELSE ! this combination of sources does not allow any calculation
            ! FF 1 with no calculation except the field
            CALL OUTPUT_FF_RESULTS(FFdata(jj)%ntheta,FFdata(jj)%nphi,      &
                   FFdata(jj)%theta0*PI/180.0,FFdata(jj)%phi0*PI/180.0,    &
                   FFdata(jj)%dtheta*PI/180.0,FFdata(jj)%dphi*PI/180.0,4)
          END IF

        CASE (2) ! Use the current incident plane wave direction
          IF ((excitation_counters(1).EQ.0).AND.(excitation_counters(2).EQ.1)) THEN ! only one external source
            ! FF 2 with monostatic RCS calculation
            CALL OUTPUT_FF_RESULTS(1,1,t_pw,p_pw,0.0,0.0,1,mag_pwE)
          ELSE IF (excitation_counters(2).GT.0) THEN ! this combination of sources, including some plane wave
                                                     ! source(s) does not allow any calculation under FF 2
            ! FF 2 with no calculation except the field at the last incident direction
            CALL OUTPUT_FF_RESULTS(1,1,t_pw,p_pw,0.0,0.0,4)
          ELSE
            CALL ERROR_FEMFEKO(0,4901)
          END IF

        END SELECT

      END IF
    END DO FF_loop

    CALL DATE_AND_TIME (values=time2)
    FF_time = FF_time + TIME_DIFFERENCE(time1,time2)

  END SUBROUTINE FEMFEKO_OUTPUT_FF_RESULTS
!*******************************************************************************

  SUBROUTINE FEMFEKO_POSTPROCESS(this_analysis)
! MMB added 9 Apr
    USE coax_feed
! MMB end
    USE cbaa_sys, ONLY: CBAA_COAX_ZIN_CALC, CBAA_PROBE_VOLTAGE_CALC, &
                            CBAA_COAX_POSTPRO
    USE gw_sys, ONLY: GW_S_PARAMS
!*******************************************************************************
! Do miscellaneous post processing.
!*******************************************************************************
    INTEGER(I4B), INTENT(IN) :: this_analysis          ! analysis number to be processed

    INTEGER(I4B) :: jj       ! counters
    
    CALL DATE_AND_TIME (values=time1)

    ! Do coax post processing:
	IF (CBAA_ANALYSIS) THEN
      IF (NUM_AX_cards.GT.0) THEN
        IF (COAX_WHITNEY) THEN
          CALL CBAA_COAX_ZIN_CALC
        ELSE     
          CALL CBAA_COAX_POSTPRO
        END IF
      END IF
    END IF

    ! Calculate the voltages of all current probes in this analysis:
    A8_pp_loop: DO jj = 1,NUM_A8_cards
      IF (.NOT.CBAA_ANALYSIS) THEN
        print *,'WARNING: A8 card ignored in this analysis type'
        CYCLE A8_pp_loop
      END IF
      IF (this_analysis.EQ.A8data(jj)%analysis_no)  THEN ! this card is part of the analysis
        CALL CBAA_PROBE_VOLTAGE_CALC(A8data(jj)%prx0,A8data(jj)%pry0,A8data(jj)%prz0, &
                                     A8data(jj)%prdir,A8data(jj)%prlen,               &
                                     A8data(jj)%prabs,A8data(jj)%prphase*PI/180.0,    &
                                     A8data(jj)%prvolt)
      END IF
    END DO A8_pp_loop

    ! Print S parameters, if requested. Calculate S11 and S12:

! DBD change 1 May 02
    IF (GW_ANALYSIS.AND.GW_COMPUTE_S_PARAMS) THEN
      CALL GW_OUTPUT_S_PARAMS(freqcount)
    ELSE IF (GW_ANALYSIS.AND..NOT.GW_COMPUTE_S_PARAMS) THEN
	  CALL GW_OUTPUT_FIELDS(freqcount)
	END IF
! End DBD change 1 May 02


    CALL DATE_AND_TIME (values=time2)
    postcalc_time = postcalc_time + TIME_DIFFERENCE(time1,time2)

  END SUBROUTINE FEMFEKO_POSTPROCESS
!*******************************************************************************

  SUBROUTINE INIT_TIMERS
    ! initialise timers:
    admin_time      = 0.0
    total_time      = 0.0
    reading_time    = 0.0
    geoprocess_time = 0.0
    matbuild_time   = 0.0
    exbuild_time    = 0.0
    solution_time   = 0.0
    postcalc_time   = 0.0
    FE_time         = 0.0
    FF_time         = 0.0
    ! -----------------------------------------------------------------
  END SUBROUTINE INIT_TIMERS
 
  SUBROUTINE MODE_LOOPS
    USE write_matrix
    ! Loop over the waveguide modes - or execute loop once if not 
    ! guided wave analysis. 
    MODE_M_LOOP: DO mode_m = 1,num_wg_modes_m
       MODE_N_LOOP: DO mode_n = 0,num_wg_modes_n

          ! Set up the system matrix and do a direct solve 
          ! Print the material properties:
          CALL output_materials_info
          ! Set up the system matrix:
          CALL DATE_AND_TIME (values=time1)
          IF (CBAA_ANALYSIS) CALL CBAA_SYSMAT
          IF (GW_ANALYSIS)   CALL GW_SYSMAT
          IF (FD_SCAT_ANALYSIS) CALL FD_SCAT_SYSMAT

          CALL DATE_AND_TIME (values=time2)
          matbuild_time = matbuild_time + TIME_DIFFERENCE(time1,time2)

          ! Solve the system directly if applicable:
          SELECT CASE(SOLVER_TYPE)
          CASE (0)
             IF ((CBAA_ANALYSIS.OR.GW_ANALYSIS.OR.FD_SCAT_ANALYSIS).AND.EXECUTE) THEN
                WRITE (*,'(A)',ADVANCE='NO') 'LU decomposition in progress...'
                CALL DIRECT_SOLVE(0,temp_time) ! LU decompose the matrix <A_mat_c>
                WRITE (*,'(A)') 'FINISHED3'
             END IF
             solution_time = solution_time + temp_time
          END SELECT

          ! First find the maximum NThei*NPhii A0 card, which must be used for the A0 loops:
          theta_total_points = 1
          phi_total_points   = 1 ! initialize, because if no A0's exist, the A0 loops must execute only once
          DO jj = 1,NUM_A0_cards
             ! All cards' NTheii*NPhii are either =1/=max. value, as ensured in 
             ! MESHIN.f90, thus one only have to chech for NTheii*NPhii > 1:
             IF ((A0data(jj)%NThei*A0data(jj)%NPhii.GT.1)) THEN
                theta_total_points = A0data(jj)%NThei
                phi_total_points   = A0data(jj)%NPhii
             END IF
          END DO
          ! Loop over incident excitation angles:
          CALL A0_LOOPS
       END DO MODE_N_LOOP
    END DO MODE_M_LOOP

    IF (DEBUG_WRITE_MATRIX .AND. SPARSE) THEN
       ! Write out the sparse A matrix. Note, this will only work for complex,
       ! sparese matrices. And only the last matrix built will be written.
       CALL write_sparse_c(Asparse_c, col_ind, row_ind)
    END IF
    
  END SUBROUTINE MODE_LOOPS

  SUBROUTINE PROCESS_ANALYSIS_GEO
    ! Other analysis-specific, geometry pre-processing:
    IF (GW_ANALYSIS)   CALL GW_PREPROCESSING ! must be called BEFORE NUMBER_DOF
    CALL INIT_GEOMETRY_ATTRIBUTES
    CALL PEC_OBJECT_SEARCH
    IF (CBAA_ANALYSIS) CALL CBAA_BOUNDARY_SEARCH
    IF (CBAA_ANALYSIS) CALL COAX_BOUNDARY_SEARCH(COAX_WHITNEY)
    IF (GW_ANALYSIS)   CALL PORT_BOUNDARY_SEARCH
    ! DBD addition 25 March and 21 May 2003   
    IF (TD_ANALYSIS) THEN 
       CALL ABC_BOUNDARY_SEARCH
       IF (SCAT_FIELD) CALL SCAT_TOTAL_BOUNDARY_SEARCH ! See note in routine. Must
       ! call this routine before 
       ! PML_ASSIGN_PROPERTIES
    END IF
    ! End DBD addition 25 March and 21 May 2003   
! DBD addition 16 Dec 2003   
   IF (FD_SCAT_ANALYSIS) THEN 
     CALL ABC_BOUNDARY_SEARCH
     IF (SCAT_FIELD) THEN   ! using the scattered field formulation
       IF (PEC_SPHERE_SCAT) THEN 
         CALL SCAT_TOTAL_BOUNDARY_SEARCH_SPH 
               IF (CURVILINEAR) THEN
           CALL CURVILINEAR_BOUNDARY_SEARCH_SPH
                 STOP 'UNIMPLEMENTED CODE IN ROUTINE FEMFEKO; EXECUTION STOPPED'
                 ! Note: this implementation is presently incomplete for the scattered field formulation
                 ! The prescribed BC's will also need to take curvature into account.
         END IF
           ELSE IF (PEC_GENERAL_SCAT) THEN
         CALL SCAT_TOTAL_BOUNDARY_SEARCH 
       END IF
       CALL SCAT_TOTAL_BOUNDARY_CHECK
         CALL PEC_SCATTERER_SEARCH
     ELSE                         ! using the total field formulation
       IF(CURVILINEAR.AND.PEC_SPHERE_SCAT) THEN
         CALL CURVILINEAR_BOUNDARY_SEARCH_SPH 
       END IF
         CALL SCAT_TOTAL_BOUNDARY_CHECK 
         IF (PEC_SPHERE_SCAT.OR.PEC_GENERAL_SCAT) THEN
         CALL PEC_SCATTERER_SEARCH ! Both spherical and general PEC scatterers are treated in the same fashion.
         ELSE IF (PENETRABLE_SCAT) THEN
           CONTINUE ! No action is needed, code proceeds as usual for FEM.
         ELSE
	     STOP 'IE: Unsupported analysis for total field formulation of EM scattering'
          END IF
       END IF
    END IF
    ! Do not add PEC boundaries in TD case using ABC's or in FD open region 
    ! case ditto.
    IF (.NOT.TD_ANALYSIS.AND..NOT.FD_SCAT_ANALYSIS) CALL PEC_BOUNDARY_SEARCH 

    ! End DBD addition 16 Dec 2003   
    IF (FD_SCAT_ANALYSIS) THEN
       CALL FREE_DIR_ASSIGN_SEARCH
    ELSE
       CALL FREE_ASSIGN_SEARCH
    END IF
    IF (CBAA_ANALYSIS) CALL CBAA_PREPROCESSING
    CALL QUAD_LINE_RULES_MAKE ! Added DBD 04 June 2003
    CALL QUAD_TRI_RULES_MAKE
    CALL CUBE_TET_RULES_MAKE ! Added DBD 21 Feb 02
  END SUBROUTINE PROCESS_ANALYSIS_GEO

  SUBROUTINE PROCESS_GEOMETRY
    USE output_mesh
    USE parseinput, ONLY: meshout_filename
    ! General processing of the geometry data:
    CALL DATE_AND_TIME (values=time1)
    WRITE (*,'(A)') 'Processing the geometry data...'
    CALL MESHSORT
    CALL NODEELEMENT_INDEXLIST_MAKE
    CALL FAST_EDGEMAKE
    IF (CBAA_ANALYSIS) CALL CBAA_EDGESORT
    CALL FAST_FACEMAKE 
    CALL FAST_FACE_CONNECT     
    IF (SPARSE.OR.CBAA_ANALYSIS) CALL FAST_EDGE_CONNECT ! Needed for any sparse scheme, and always for CBAA.
    CALL FAST_CONNECT_ELEMENT_TO_FACE
!!$    CALL NODEELEMENT_INDEXLIST_CLEAN
    
    IF (WRITE_MESH) THEN ! Write out mesh in new mesh format, then quit
       OPEN(unit=MESHOUT, status='NEW', file=meshout_filename)
       CALL output_mesh_nodes(MESHOUT)
       CLOSE(unit=meshout)
       PRINT*, 'Mesh written to file: ', meshout_filename
       STOP 
    END IF
  END SUBROUTINE PROCESS_GEOMETRY



  SUBROUTINE READ_FILES
    ! Reading the input and mesh file(s):
    USE parseinput
    USE input_mesh
    CALL DATE_AND_TIME (values=time1)
    WRITE (*,'(A)',ADVANCE='NO') 'Reading the input file(s)......'
    OPEN (UNIT=INFILE,STATUS='OLD',file=infilename, iostat=ios)
    CALL read_namelists(infile)
    REWIND(infile)
    SELECT CASE(mesh_format)
    CASE('oldfek')
       CALL MESHIN_FEK_OLD(infilename)
    CASE('oldfeknative')
       ! Bit of a hack to allow the .fek file to be read for control, and 
       ! the new format mesh for geometry
       CALL meshin_fek_old(infilename)
       CLOSE(meshin)
       CALL input_mesh_read(extra_meshfilename)
    CASE('newfeknative')
       ! Bit of a hack to allow the .fek file to be read for control, and 
       ! the new format mesh for geometry
       CALL meshin_fek(infilename)
       CLOSE(meshin)
       CALL input_mesh_read(extra_meshfilename)
    CASE('newfek')
       CALL MESHIN_FEK(infilename)
    CASE('native')
       STOP 'New native mesh reading not yet supported!'
    CASE default
       PRINT*, 'Unknown mesh format, ', mesh_type, ' specified'
       STOP
    END SELECT

    CALL output_mesh_info
    WRITE (*,'(A)') 'FINISHED'
    PRINT*, 'COMPUTE_EIGENMODES:          ', COMPUTE_EIGENMODES
    CALL DATE_AND_TIME (values=time2)
    reading_time = reading_time + TIME_DIFFERENCE(time1,time2)
    ! -----------------------------------------------------------------
  END SUBROUTINE READ_FILES

  SUBROUTINE S_PARAMS_CALC
!*********************************************************************
! Loops over all the ports, and sets up appropriate excitation vectors,
! solves the system and calculates the S parameters
!*********************************************************************
    USE gw_sys, ONLY: GW_S_PARAMS
    IMPLICIT NONE
    DO port_counter = 1,num_ports

       ! Set up the b-vector (excitation):
       CALL FEMFEKO_EXCITATION_VECTOR(analysis_count,cardi,cardj,&
            port_counter)
       ! Solve the matrix equation:
       SELECT CASE (SOLVER_TYPE)
       CASE (0)
          IF (CBAA_ANALYSIS.OR.GW_ANALYSIS) CALL DIRECT_SOLVE(1,temp_time) ! back substitution
       CASE (1:5)
          IF (CBAA_ANALYSIS.OR.GW_ANALYSIS) CALL ITER_SOLVE(temp_time)
       END SELECT
       solution_time = solution_time + temp_time
       ! Compute S params for this port.
       CALL GW_S_PARAMS(port_counter) 
    END DO
    ! Post processing:
    CALL FEMFEKO_POSTPROCESS(analysis_count)

    ! Note that for S-parameter analysis, near fields are not
    ! computed. 
    IF (NUM_FE_cards.NE.0) CALL ERROR_FEMFEKO(0,4902) ! Print warning.
  END SUBROUTINE S_PARAMS_CALC
END MODULE FEMFEKO



