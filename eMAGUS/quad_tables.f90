! Last changed 21 Feb 2002, DBD. Some typos corrected. 
! Tetrahedral rules added. 

MODULE quad_tables
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This module contains datastractures that store quadrature (2D numerical 
! integration) and cubature (3D numerical integration) rules for various 
! geometries.
! Created: Oct 2001, MMB.
! Extended to include 3D rules: Feb 2002, DBD.
! Extended to include 1D rules: 04 Jun 2003, DBD.
!*******************************************************************************

  ! Datastructures for line rules:
  TYPE QUAD_LINE_RULE
    REAL(SP), DIMENSION(:,:), POINTER :: rule
  END TYPE QUAD_LINE_RULE
  TYPE(QUAD_LINE_RULE), DIMENSION(:), ALLOCATABLE :: quad_line_rules

  INTEGER(I4B) gauss_line_points      ! Number of quad. points to be used for 
                                      ! quadrature on lines. 


  ! Datastructures for triangle rules:
  TYPE QUAD_TRI_RULE
    REAL(SP), DIMENSION(:,:), POINTER :: rule
  END TYPE QUAD_TRI_RULE
  TYPE(QUAD_TRI_RULE), DIMENSION(:), ALLOCATABLE :: quad_tri_rules

  INTEGER(I4B) gauss_points      ! Number of quad. points to be used for 
                                 ! quadrature on triangles. 

  ! Datastructures for tetrahedral rules:
  TYPE CUBE_TET_RULE
    REAL(SP), DIMENSION(:,:), POINTER :: rule
  END TYPE CUBE_TET_RULE
  TYPE(CUBE_TET_RULE), DIMENSION(:), ALLOCATABLE :: cube_tet_rules

  INTEGER(I4B) tet_gauss_points  ! Number of cubature. points to be used for 
                                 ! cubature on tetrahedrons.
								 
!*******************************************************************************
! Internal routines:
! (A routine is either listed as an interface, or in the private statement.)

INTERFACE QUAD_LINE_RULES_CLEAN
  MODULE PROCEDURE QUAD_LINE_RULES_CLEAN
END INTERFACE

INTERFACE QUAD_LINE_RULES_MAKE
  MODULE PROCEDURE QUAD_LINE_RULES_MAKE
END INTERFACE

INTERFACE QUAD_TRI_RULES_CLEAN
  MODULE PROCEDURE QUAD_TRI_RULES_CLEAN
END INTERFACE

INTERFACE QUAD_TRI_RULES_MAKE
  MODULE PROCEDURE QUAD_TRI_RULES_MAKE
END INTERFACE

INTERFACE CUBE_TET_RULES_CLEAN
  MODULE PROCEDURE CUBE_TET_RULES_CLEAN
END INTERFACE

INTERFACE CUBE_TET_RULES_MAKE
  MODULE PROCEDURE CUBE_TET_RULES_MAKE
END INTERFACE


!*******************************************************************************

CONTAINS
!*******************************************************************************


SUBROUTINE QUAD_LINE_RULES_CLEAN
  IMPLICIT NONE
!*******************************************************************************
! Deallocates the storage of the triagular quadrature rules.
!*******************************************************************************

  DEALLOCATE(quad_line_rules(1)%rule)
  DEALLOCATE(quad_line_rules(2)%rule)
  DEALLOCATE(quad_line_rules(3)%rule)
  DEALLOCATE(quad_line_rules(4)%rule)
  
  DEALLOCATE(quad_line_rules)

END SUBROUTINE QUAD_LINE_RULES_CLEAN
!*******************************************************************************

SUBROUTINE QUAD_LINE_RULES_MAKE
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine assignes the integration point data for a line (1D).
! Structure of one row: [simplex1 simplex2 weight]
! Adapted from the Gauss-Legendre rules given in J Jin, 
! "The Finite Element Method in Electromagnetics", 2nd edn, 2002. 
! The rules are for integration from -1 to +1. The simplex integration is 
! from 0 to 1. The linear mapping which accomplishes this change of variables is:
!
! simplex1 = xi/2+1/2. 
!
! From the properties of the simplex coordinates, 
! simplex2 = 1-simplex1.
!
! Note also the weights must be divided by 2.  
!*******************************************************************************

  ! First allocate the space for the rules in general (1:max number of points):
  ALLOCATE(quad_line_rules(1:4))

  ! Assign the 1 point rule data, degree of precision 1. 
  ALLOCATE(quad_line_rules(1)%rule(1,3))
  quad_line_rules(1)%rule(1,1:3) = (/ 0.5, 0.5, 2.0/2.0_SP /)

  ! Assign the 2 point rule data: degree of precision 3
  ALLOCATE(quad_line_rules(2)%rule(2,3))
  quad_line_rules(2)%rule(1,1:3) = (/ 0.2113248654,  0.7886751346, 1.0/2.0_SP /)
  quad_line_rules(2)%rule(2,1:3) = (/ 0.7886751346,  0.2113248654, 1.0/2.0_SP /)

  ! Assign the 3 point rule data: degree of precision 5
  ALLOCATE(quad_line_rules(3)%rule(3,3))
  quad_line_rules(3)%rule(1,1:3) = (/ 0.1127016654,   0.8872983346, 0.5555555556/2.0_SP /)
  quad_line_rules(3)%rule(2,1:3) = (/ 0.5,            0.5,          0.8888888889/2.0_SP /)
  quad_line_rules(3)%rule(3,1:3) = (/ 0.8872983346,   0.1127016654, 0.5555555556/2.0_SP /)
  
  ! Assign the 4 point rule data: degree of precision 7
  ALLOCATE(quad_line_rules(4)%rule(4,3))
  quad_line_rules(4)%rule(1,1:3) = (/ 0.0694318442,   0.9305681558, 0.3478548451/2.0_SP /) 
  quad_line_rules(4)%rule(2,1:3) = (/ 0.3300094782,   0.6699905218, 0.6521451549/2.0_SP /) 
  quad_line_rules(4)%rule(3,1:3) = (/ 0.6699905218,   0.3300094782, 0.6521451549/2.0_SP /) 
  quad_line_rules(4)%rule(4,1:3) = (/ 0.9305681558,   0.0694318442, 0.3478548451/2.0_SP /) 
  
  
END SUBROUTINE QUAD_LINE_RULES_MAKE
!*******************************************************************************

SUBROUTINE QUAD_TRI_RULES_CLEAN
  IMPLICIT NONE
!*******************************************************************************
! Deallocates the storage of the triangular quadrature rules.
!*******************************************************************************

  DEALLOCATE(quad_tri_rules(1)%rule)
  DEALLOCATE(quad_tri_rules(3)%rule)
  DEALLOCATE(quad_tri_rules(4)%rule)
  DEALLOCATE(quad_tri_rules(6)%rule)
  DEALLOCATE(quad_tri_rules(7)%rule)
  DEALLOCATE(quad_tri_rules(13)%rule)
  DEALLOCATE(quad_tri_rules(25)%rule)

  DEALLOCATE(quad_tri_rules)

END SUBROUTINE QUAD_TRI_RULES_CLEAN
!*******************************************************************************


SUBROUTINE QUAD_TRI_RULES_MAKE
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine assignes the integration point data for the triangular domain.
! Structure of one row: [simplex1 simplex2 simplex3 weight]
! Simplex quadrature points and weights for triangles.
! Given with 15 significant digits to accommodate higher precision in future.
! See: Cowper, "Gaussian quadrature formulas for triangles", 
! Intl. Journal Numer. Meth. Eng., 1973, pp. 95-100.
! See: Dunavant, "High degree efficient symmetrical gaussian quadrature
! rules for the triangle", Intl. Journal Numer. Meth. Eng., 1985, vol. 21, 
! pp. 1129-1148.
!*******************************************************************************

  ! First allocate the space for the rules in general (1:max number of points):
  ALLOCATE(quad_tri_rules(1:25))

  ! Assign the 1 point rule data:
  ALLOCATE(quad_tri_rules(1)%rule(1,4))
  quad_tri_rules(1)%rule(1,1:4) = (/ 0.3333333333, 0.3333333333, 0.3333333333, 1.0 /)

  ! Assign the 3 point rule data: degree of precision 2
  ALLOCATE(quad_tri_rules(3)%rule(3,4))
  quad_tri_rules(3)%rule(1,1:4) = (/ 0.6666666666, 0.1666666666, 0.1666666666, 0.33333333333 /)
  quad_tri_rules(3)%rule(2,1:4) = (/ 0.1666666666, 0.6666666666, 0.1666666666, 0.33333333333 /)
  quad_tri_rules(3)%rule(3,1:4) = (/ 0.1666666666, 0.1666666666, 0.6666666666, 0.33333333333 /)

  ! Assign the 4 point rule data: degree of precision 3
  ALLOCATE(quad_tri_rules(4)%rule(4,4))
  quad_tri_rules(4)%rule(1,1:4) = (/ 0.3333333333, 0.3333333333, 0.3333333333,      -0.5625 /) 
  quad_tri_rules(4)%rule(2,1:4) = (/          0.6,          0.2,          0.2, 0.5208333333 /) 
  quad_tri_rules(4)%rule(3,1:4) = (/          0.2,          0.6,          0.2, 0.5208333333 /) 
  quad_tri_rules(4)%rule(4,1:4) = (/          0.2,          0.2,          0.6, 0.5208333333 /) 

  ! Assign the 6 point rule data: degree of precision 4
  ALLOCATE(quad_tri_rules(6)%rule(6,4))
  quad_tri_rules(6)%rule(1,1:4) = (/ 0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322 /)
  quad_tri_rules(6)%rule(2,1:4) = (/ 0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322 /)
  quad_tri_rules(6)%rule(3,1:4) = (/ 0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322 /)
  quad_tri_rules(6)%rule(4,1:4) = (/ 0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011 /)
  quad_tri_rules(6)%rule(5,1:4) = (/ 0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011 /)
  quad_tri_rules(6)%rule(6,1:4) = (/ 0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011 /)

  ! Assign the 7 point rule data: [Dunavant] degree of precision 5
  ALLOCATE(quad_tri_rules(7)%rule(7,4))
  quad_tri_rules(7)%rule(1,1:4) = (/ 0.333333333, 0.333333333, 0.333333333,        0.225 /) 
  quad_tri_rules(7)%rule(2,1:4) = (/ 0.059715872, 0.470142064, 0.470142064,  0.132394153 /) 
  quad_tri_rules(7)%rule(3,1:4) = (/ 0.470142064, 0.059715872, 0.470142064,  0.132394153 /) 
  quad_tri_rules(7)%rule(4,1:4) = (/ 0.470142064, 0.470142064, 0.059715872,  0.132394153 /) 
  quad_tri_rules(7)%rule(5,1:4) = (/ 0.797426985,  0.10128651,  0.10128651, 0.1259391805 /) 
  quad_tri_rules(7)%rule(6,1:4) = (/  0.10128651, 0.797426985,  0.10128651, 0.1259391805 /) 
  quad_tri_rules(7)%rule(7,1:4) = (/  0.10128651,  0.10128651, 0.797426985, 0.1259391805 /) 

  ! Assign the 13 point rule data: [Dunavant] degree of precision 7
  ALLOCATE(quad_tri_rules(13)%rule(13,4))
  quad_tri_rules(13)%rule(1,1:4)  = (/    0.3333333333,    0.3333333333,    0.3333333333, -0.1495700444677 /) 
  quad_tri_rules(13)%rule(2,1:4)  = (/    0.4793080678,   0.26034596608,   0.26034596608,   0.175615257433 /) 
  quad_tri_rules(13)%rule(3,1:4)  = (/   0.26034596608,    0.4793080678,   0.26034596608,   0.175615257433 /) 
  quad_tri_rules(13)%rule(4,1:4)  = (/   0.26034596608,   0.26034596608,    0.4793080678,   0.175615257433 /) 
  quad_tri_rules(13)%rule(5,1:4)  = (/  0.869739794196,  0.065130102902,  0.065130102902,   0.053347235609 /) 
  quad_tri_rules(13)%rule(6,1:4)  = (/  0.065130102902,  0.869739794196,  0.065130102902,   0.053347235609 /) 
  quad_tri_rules(13)%rule(7,1:4)  = (/  0.065130102902,  0.065130102902,  0.869739794196,   0.053347235609 /) 
  quad_tri_rules(13)%rule(8,1:4)  = (/  0.048690315425,  0.312865496005, 0.6384441885698,  0.0771137608903 /)
  quad_tri_rules(13)%rule(9,1:4)  = (/  0.048690315425, 0.6384441885698,  0.312865496005,  0.0771137608903 /) 
  quad_tri_rules(13)%rule(10,1:4) = (/ 0.6384441885698,  0.048690315425,  0.312865496005,  0.0771137608903 /) 
  quad_tri_rules(13)%rule(11,1:4) = (/  0.312865496005,  0.048690315425, 0.6384441885698,  0.0771137608903 /) 
  quad_tri_rules(13)%rule(12,1:4) = (/  0.312865496005, 0.6384441885698,  0.048690315425,  0.0771137608903 /) 
  quad_tri_rules(13)%rule(13,1:4) = (/ 0.6384441885698,  0.312865496005,  0.048690315425,  0.0771137608903 /) 

  ! Assign the 25 point rule data: [Dunavant] degree of precision 10
  ALLOCATE(quad_tri_rules(25)%rule(25,4))
  quad_tri_rules(25)%rule(1,1:4)  = (/   0.33333333333,   0.33333333333,   0.33333333333,     0.0908179904 /) 
  quad_tri_rules(25)%rule(2,1:4)  = (/  0.028844733233,  0.485577633384,  0.485577633384,  0.0367259577565 /) 
  quad_tri_rules(25)%rule(3,1:4)  = (/  0.485577633384,  0.028844733233,  0.485577633384,  0.0367259577565 /) 
  quad_tri_rules(25)%rule(4,1:4)  = (/  0.485577633384,  0.485577633384,  0.028844733233,  0.0367259577565 /) 
  quad_tri_rules(25)%rule(5,1:4)  = (/   0.78103684903,  0.109481575485,  0.109481575485, 0.04532105943553 /) 
  quad_tri_rules(25)%rule(6,1:4)  = (/  0.109481575485,   0.78103684903,  0.109481575485, 0.04532105943553 /) 
  quad_tri_rules(25)%rule(7,1:4)  = (/  0.109481575485,  0.109481575485,   0.78103684903, 0.04532105943553 /) 
  quad_tri_rules(25)%rule(8,1:4)  = (/ 0.1417072194149, 0.3079398387641,  0.550352941821,  0.0727579168454 /)
  quad_tri_rules(25)%rule(9,1:4)  = (/ 0.1417072194149,  0.550352941821, 0.3079398387641,  0.0727579168454 /) 
  quad_tri_rules(25)%rule(10,1:4) = (/ 0.3079398387641, 0.1417072194149,  0.550352941821,  0.0727579168454 /) 
  quad_tri_rules(25)%rule(11,1:4) = (/  0.550352941821, 0.1417072194149, 0.3079398387641,  0.0727579168454 /) 
  quad_tri_rules(25)%rule(12,1:4) = (/ 0.3079398387641,  0.550352941821, 0.1417072194149,  0.0727579168454 /) 
  quad_tri_rules(25)%rule(13,1:4) = (/  0.550352941821, 0.3079398387641, 0.1417072194149,  0.0727579168454 /) 
  quad_tri_rules(25)%rule(14,1:4) = (/ 0.0250035347627, 0.2466725606399, 0.7283239045974, 0.02832724253106 /)
  quad_tri_rules(25)%rule(15,1:4) = (/ 0.0250035347627, 0.7283239045974, 0.2466725606399, 0.02832724253106 /)
  quad_tri_rules(25)%rule(16,1:4) = (/ 0.2466725606399, 0.0250035347627, 0.7283239045974, 0.02832724253106 /)
  quad_tri_rules(25)%rule(17,1:4) = (/ 0.7283239045974, 0.0250035347627, 0.2466725606399, 0.02832724253106 /)
  quad_tri_rules(25)%rule(18,1:4) = (/ 0.2466725606399, 0.7283239045974, 0.0250035347627, 0.02832724253106 /)
  quad_tri_rules(25)%rule(19,1:4) = (/ 0.7283239045974, 0.2466725606399, 0.0250035347627, 0.02832724253106 /)
  quad_tri_rules(25)%rule(20,1:4) = (/ 0.0095408154003, 0.0668032510122, 0.9236559335875, 0.00942166696373 /)
  quad_tri_rules(25)%rule(21,1:4) = (/ 0.0095408154003, 0.9236559335875, 0.0668032510122, 0.00942166696373 /)
  quad_tri_rules(25)%rule(22,1:4) = (/ 0.0668032510122, 0.0095408154003, 0.9236559335875, 0.00942166696373 /)
  quad_tri_rules(25)%rule(23,1:4) = (/ 0.9236559335875, 0.0095408154003, 0.0668032510122, 0.00942166696373 /)
  quad_tri_rules(25)%rule(24,1:4) = (/ 0.0668032510122, 0.9236559335875, 0.0095408154003, 0.00942166696373 /)
  quad_tri_rules(25)%rule(25,1:4) = (/ 0.9236559335875, 0.0668032510122, 0.0095408154003, 0.00942166696373 /)

END SUBROUTINE QUAD_TRI_RULES_MAKE
!*******************************************************************************

SUBROUTINE CUBE_TET_RULES_CLEAN
  IMPLICIT NONE
!*******************************************************************************
! Deallocates the storage of the triagular quadrature rules.
!*******************************************************************************

  DEALLOCATE(cube_tet_rules(4)%rule)
  DEALLOCATE(cube_tet_rules(11)%rule)
  DEALLOCATE(cube_tet_rules(24)%rule)
  DEALLOCATE(cube_tet_rules(45)%rule)
  DEALLOCATE(cube_tet_rules)

END SUBROUTINE CUBE_TET_RULES_CLEAN
!*******************************************************************************

SUBROUTINE CUBE_TET_RULES_MAKE
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This routine assignes the integration point data for the tetrahedral domain.
! Structure of one row: [simplex1 simplex2 simplex3 simplex 4 weight]
! Symmetric simplex quadrature points and weights for tetrahedrons.
! Given with as many significant digits as available in the references, 
! to accommodate higher precision in future.
! See: Y. Jinyun, "Symmetric gaussian quadrature formulae for tetrahedronal
! regions", Computer methods in  Applied Mechanics and Engineering",
! 43, 1984, pp. 349--353. Note that eta, xi, zeta and nu here correspond to 
! alpha, beta, nu and delta in this reference respectively and 
! and P. Keast "Moderate-degree tetrahedral quadrature formulas", 
! Computer methods in  Applied Mechanics and Engineering, 55, 1986, pp.
! 339-348.
! 
! Author: DBD, 21 Feb 2002. Similar structure and use to 2D rules. 
! Extended DBD 25 Jul 2005; 6th and 8th order rules - as entered & tested by MMB - added.
!*******************************************************************************

  ! Following used for 11 point rule.
  REAL(SP), PARAMETER :: zz = 1.0_SP/4.0_SP
  REAL(SP), PARAMETER :: aa = 0.714285714285714285E-01
  REAL(SP), PARAMETER :: bb = 0.785714285714285714E+00
  REAL(SP), PARAMETER :: cc = 0.399403576166799219E+00
  REAL(SP), PARAMETER :: dd = 0.100596423833200785E+00
  REAL(SP), DIMENSION(5) :: qpoint1,qpoint2,qpoint3,qpoint4,qpoint5,qpoint6,qpoint7

  ! First allocate the space for the rules in general (1:max number of points):
  ALLOCATE(cube_tet_rules(1:45))

  ALLOCATE(cube_tet_rules(4)%rule(4,5)) 
  ! 4 point rule, degree of precision 2 [Jinyun]
  cube_tet_rules(4)%rule(1,1:5) = (/0.5854101966249685, 0.138196601125015, 0.138196601125015, 0.138196601125015, 1.0_SP/4.0_SP /)
  cube_tet_rules(4)%rule(2,1:5) = (/0.138196601125015, 0.5854101966249685, 0.138196601125015, 0.138196601125015, 1.0_SP/4.0_SP /)
  cube_tet_rules(4)%rule(3,1:5) = (/0.138196601125015, 0.138196601125015, 0.5854101966249685, 0.138196601125015, 1.0_SP/4.0_SP /)
  cube_tet_rules(4)%rule(4,1:5) = (/0.138196601125015, 0.138196601125015, 0.138196601125015, 0.5854101966249685, 1.0_SP/4.0_SP /)

  ALLOCATE(cube_tet_rules(11)%rule(11,5))   
  ! Eleven point rule; degree of precision 4 [Keast; 1st N=4 rule]
  ! Note that Keast's rule omit a factor of 6 (volume of a "unitary" tet). 
  cube_tet_rules(11)%rule(1,1:5) =  (/ zz, zz, zz, zz, 6.0_SP*(-0.1315555555555555550E-01) /)   
  cube_tet_rules(11)%rule(2,1:5) =  (/ aa, aa, aa, bb, 6.0_SP*0.7622222222222222222E-02 /)
  cube_tet_rules(11)%rule(3,1:5) =  (/ aa, aa, bb, aa, 6.0_SP*0.7622222222222222222E-02 /)
  cube_tet_rules(11)%rule(4,1:5) =  (/ aa, bb, aa, aa, 6.0_SP*0.7622222222222222222E-02 /)
  cube_tet_rules(11)%rule(5,1:5) =  (/ bb, aa, aa, aa, 6.0_SP*0.7622222222222222222E-02 /)
  cube_tet_rules(11)%rule(6,1:5) =  (/ cc, cc, dd, dd, 6.0_SP*0.2488888888888888880E-01 /)
  cube_tet_rules(11)%rule(7,1:5) =  (/ cc, dd, cc, dd, 6.0_SP*0.2488888888888888880E-01 /)
  cube_tet_rules(11)%rule(8,1:5) =  (/ cc, dd, dd, cc, 6.0_SP*0.2488888888888888880E-01 /)
  cube_tet_rules(11)%rule(9,1:5) =  (/ dd, dd, cc, cc, 6.0_SP*0.2488888888888888880E-01 /)
  cube_tet_rules(11)%rule(10,1:5) = (/ dd, cc, dd, cc, 6.0_SP*0.2488888888888888880E-01 /)
  cube_tet_rules(11)%rule(11,1:5) = (/ dd, cc, cc, dd, 6.0_SP*0.2488888888888888880E-01 /)
  

  ALLOCATE(cube_tet_rules(24)%rule(24,5))
  ! Rule added by MM Botha, 2002-2003.
  ! Assign the order 6 rule:
  ! 24 point rule; degree of precision 6 [Keast; 1986, CMAME]
  ! Note that Keast's rule omit a factor of 6 (volume of a "unitary" tet). 
  qpoint1 = (/ 0.214602871259151684_SP     , 0.214602871259151684_SP     , 0.214602871259151684_SP     , &
               0.356191386222544953_SP     , 6.0_SP*0.665379170969464506E-02_SP /)
  qpoint2 = (/ 0.406739585346113397E-01_SP , 0.406739585346113397E-01_SP , 0.406739585346113397E-01_SP , &
               0.877978124396165982_SP     , 6.0_SP*0.167953517588677620E-02_SP /)
  qpoint3 = (/ 0.322337890142275646_SP     , 0.322337890142275646_SP     , 0.322337890142275646_SP     , &
               0.329863295731730594E-01_SP , 6.0_SP*0.922619692394239843E-02_SP /)
  qpoint4 = (/ 0.636610018750175299E-01_SP , 0.636610018750175299E-01_SP , 0.269672331458315867_SP     , &
               0.603005664791649076_SP     , 6.0_SP*0.803571428571428248E-02_SP /)

  cube_tet_rules(24)%rule(1,1:5)   = (/ qpoint1(1) , qpoint1(2) , qpoint1(3) , qpoint1(4) , qpoint1(5) /)   
  cube_tet_rules(24)%rule(2,1:5)   = (/ qpoint1(2) , qpoint1(3) , qpoint1(4) , qpoint1(1) , qpoint1(5) /)   
  cube_tet_rules(24)%rule(3,1:5)   = (/ qpoint1(3) , qpoint1(4) , qpoint1(1) , qpoint1(2) , qpoint1(5) /)   
  cube_tet_rules(24)%rule(4,1:5)   = (/ qpoint1(4) , qpoint1(1) , qpoint1(2) , qpoint1(3) , qpoint1(5) /)   
  cube_tet_rules(24)%rule(5,1:5)   = (/ qpoint2(1) , qpoint2(2) , qpoint2(3) , qpoint2(4) , qpoint2(5) /)   
  cube_tet_rules(24)%rule(6,1:5)   = (/ qpoint2(2) , qpoint2(3) , qpoint2(4) , qpoint2(1) , qpoint2(5) /)   
  cube_tet_rules(24)%rule(7,1:5)   = (/ qpoint2(3) , qpoint2(4) , qpoint2(1) , qpoint2(2) , qpoint2(5) /)   
  cube_tet_rules(24)%rule(8,1:5)   = (/ qpoint2(4) , qpoint2(1) , qpoint2(2) , qpoint2(3) , qpoint2(5) /)   

  cube_tet_rules(24)%rule(9,1:5)   = (/ qpoint3(1) , qpoint3(2) , qpoint3(3) , qpoint3(4) , qpoint3(5) /)   
  cube_tet_rules(24)%rule(10,1:5)  = (/ qpoint3(2) , qpoint3(3) , qpoint3(4) , qpoint3(1) , qpoint3(5) /)   
  cube_tet_rules(24)%rule(11,1:5)  = (/ qpoint3(3) , qpoint3(4) , qpoint3(1) , qpoint3(2) , qpoint3(5) /)   
  cube_tet_rules(24)%rule(12,1:5)  = (/ qpoint3(4) , qpoint3(1) , qpoint3(2) , qpoint3(3) , qpoint3(5) /)   

  cube_tet_rules(24)%rule(13,1:5)  = (/ qpoint4(1) , qpoint4(2) , qpoint4(3) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(14,1:5)  = (/ qpoint4(1) , qpoint4(2) , qpoint4(4) , qpoint4(3) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(15,1:5)  = (/ qpoint4(1) , qpoint4(3) , qpoint4(2) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(16,1:5)  = (/ qpoint4(1) , qpoint4(4) , qpoint4(2) , qpoint4(3) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(17,1:5)  = (/ qpoint4(1) , qpoint4(3) , qpoint4(4) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(18,1:5)  = (/ qpoint4(1) , qpoint4(4) , qpoint4(3) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(19,1:5)  = (/ qpoint4(3) , qpoint4(1) , qpoint4(2) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(20,1:5)  = (/ qpoint4(4) , qpoint4(1) , qpoint4(2) , qpoint4(3) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(21,1:5)  = (/ qpoint4(3) , qpoint4(1) , qpoint4(4) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(22,1:5)  = (/ qpoint4(4) , qpoint4(1) , qpoint4(3) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(23,1:5)  = (/ qpoint4(3) , qpoint4(4) , qpoint4(1) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(24)%rule(24,1:5)  = (/ qpoint4(4) , qpoint4(3) , qpoint4(1) , qpoint4(2) , qpoint4(5) /)   

  ALLOCATE(cube_tet_rules(45)%rule(45,5))
  ! Rule added by MM Botha, 2002-2003.
  ! Assign the order 8 rule:
  ! 45 point rule; degree of precision 8 [Keast; 1986, CMAME]
  ! Note that Keast's rule omit a factor of 6 (volume of a "unitary" tet). 
  ! !!!!*N*O*T*E*!!!! The negative weight on the first point is indicated positive in the reference,
  ! but in the text it is indicated that the rule does contain negative weights. Bu setting that one negative,
  ! the rule sums to unity. SO: hopefully this is correct....!
  qpoint1 = (/ 0.25_SP , 0.25_SP , 0.25_SP , 0.25_SP , -6.0_SP*(0.393270066412926145E-01_SP) /) 
  qpoint2 = (/ 0.127470936566639015_SP , 0.127470936566639015_SP , 0.127470936566639015_SP , &
               0.617587190300082967_SP , 6.0_SP*(0.408131605934270525E-02_SP) /) 
  qpoint3 = (/ 0.320788303926322960E-01_SP , 0.320788303926322960E-01_SP , 0.320788303926322960E-01_SP , &
               0.903763508822103123_SP     , 6.0_SP*(0.658086773304341943E-03_SP) /) 
  qpoint4 = (/ 0.497770956432810185E-01_SP , 0.497770956432810185E-01_SP , 0.450222904356718978_SP , &
               0.450222904356718978_SP     , 6.0_SP*(0.438425882512284693E-02_SP) /) 
  qpoint5 = (/ 0.183730447398549945_SP , 0.183730447398549945_SP , 0.316269552601450060_SP , &
               0.316269552601450060_SP , 6.0_SP*(0.138300638425098166E-01_SP) /) 
  qpoint6 = (/ 0.231901089397150906_SP , 0.231901089397150906_SP , 0.229177878448171174E-01_SP , &
               0.513280033360881072_SP , 6.0_SP*(0.424043742468372453E-02_SP) /) 
  qpoint7 = (/ 0.379700484718286102E-01_SP , 0.379700484718286102E-01_SP , 0.730313427807538396_SP , &
               0.193746475248804382_SP     , 6.0_SP*(0.223873973961420164E-02_SP) /) 

  cube_tet_rules(45)%rule(1,1:5)   = (/ qpoint1(1) , qpoint1(2) , qpoint1(3) , qpoint1(4) , qpoint1(5) /)   

  cube_tet_rules(45)%rule(2,1:5)   = (/ qpoint2(1) , qpoint2(2) , qpoint2(3) , qpoint2(4) , qpoint2(5) /)   
  cube_tet_rules(45)%rule(3,1:5)   = (/ qpoint2(2) , qpoint2(3) , qpoint2(4) , qpoint2(1) , qpoint2(5) /)   
  cube_tet_rules(45)%rule(4,1:5)   = (/ qpoint2(3) , qpoint2(4) , qpoint2(1) , qpoint2(2) , qpoint2(5) /)   
  cube_tet_rules(45)%rule(5,1:5)   = (/ qpoint2(4) , qpoint2(1) , qpoint2(2) , qpoint2(3) , qpoint2(5) /)   

  cube_tet_rules(45)%rule(6,1:5)   = (/ qpoint3(1) , qpoint3(2) , qpoint3(3) , qpoint3(4) , qpoint3(5) /)   
  cube_tet_rules(45)%rule(7,1:5)   = (/ qpoint3(2) , qpoint3(3) , qpoint3(4) , qpoint3(1) , qpoint3(5) /)   
  cube_tet_rules(45)%rule(8,1:5)   = (/ qpoint3(3) , qpoint3(4) , qpoint3(1) , qpoint3(2) , qpoint3(5) /)   
  cube_tet_rules(45)%rule(9,1:5)   = (/ qpoint3(4) , qpoint3(1) , qpoint3(2) , qpoint3(3) , qpoint3(5) /)   

  cube_tet_rules(45)%rule(10,1:5)  = (/ qpoint4(1) , qpoint4(2) , qpoint4(3) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(45)%rule(11,1:5)  = (/ qpoint4(1) , qpoint4(3) , qpoint4(2) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(45)%rule(12,1:5)  = (/ qpoint4(1) , qpoint4(3) , qpoint4(4) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(45)%rule(13,1:5)  = (/ qpoint4(3) , qpoint4(1) , qpoint4(2) , qpoint4(4) , qpoint4(5) /)   
  cube_tet_rules(45)%rule(14,1:5)  = (/ qpoint4(3) , qpoint4(1) , qpoint4(4) , qpoint4(2) , qpoint4(5) /)   
  cube_tet_rules(45)%rule(15,1:5)  = (/ qpoint4(3) , qpoint4(4) , qpoint4(1) , qpoint4(2) , qpoint4(5) /)   

  cube_tet_rules(45)%rule(16,1:5)  = (/ qpoint5(1) , qpoint5(2) , qpoint5(3) , qpoint5(4) , qpoint5(5) /)   
  cube_tet_rules(45)%rule(17,1:5)  = (/ qpoint5(1) , qpoint5(3) , qpoint5(2) , qpoint5(4) , qpoint5(5) /)   
  cube_tet_rules(45)%rule(18,1:5)  = (/ qpoint5(1) , qpoint5(3) , qpoint5(4) , qpoint5(2) , qpoint5(5) /)   
  cube_tet_rules(45)%rule(19,1:5)  = (/ qpoint5(3) , qpoint5(1) , qpoint5(2) , qpoint5(4) , qpoint5(5) /)   
  cube_tet_rules(45)%rule(20,1:5)  = (/ qpoint5(3) , qpoint5(1) , qpoint5(4) , qpoint5(2) , qpoint5(5) /)   
  cube_tet_rules(45)%rule(21,1:5)  = (/ qpoint5(3) , qpoint5(4) , qpoint5(1) , qpoint5(2) , qpoint5(5) /)   

  cube_tet_rules(45)%rule(22,1:5)  = (/ qpoint6(1) , qpoint6(2) , qpoint6(3) , qpoint6(4) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(23,1:5)  = (/ qpoint6(1) , qpoint6(2) , qpoint6(4) , qpoint6(3) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(24,1:5)  = (/ qpoint6(1) , qpoint6(3) , qpoint6(2) , qpoint6(4) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(25,1:5)  = (/ qpoint6(1) , qpoint6(4) , qpoint6(2) , qpoint6(3) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(26,1:5)  = (/ qpoint6(1) , qpoint6(3) , qpoint6(4) , qpoint6(2) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(27,1:5)  = (/ qpoint6(1) , qpoint6(4) , qpoint6(3) , qpoint6(2) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(28,1:5)  = (/ qpoint6(3) , qpoint6(1) , qpoint6(2) , qpoint6(4) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(29,1:5)  = (/ qpoint6(4) , qpoint6(1) , qpoint6(2) , qpoint6(3) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(30,1:5)  = (/ qpoint6(3) , qpoint6(1) , qpoint6(4) , qpoint6(2) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(31,1:5)  = (/ qpoint6(4) , qpoint6(1) , qpoint6(3) , qpoint6(2) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(32,1:5)  = (/ qpoint6(3) , qpoint6(4) , qpoint6(1) , qpoint6(2) , qpoint6(5) /)   
  cube_tet_rules(45)%rule(33,1:5)  = (/ qpoint6(4) , qpoint6(3) , qpoint6(1) , qpoint6(2) , qpoint6(5) /)   

  cube_tet_rules(45)%rule(34,1:5)  = (/ qpoint7(1) , qpoint7(2) , qpoint7(3) , qpoint7(4) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(35,1:5)  = (/ qpoint7(1) , qpoint7(2) , qpoint7(4) , qpoint7(3) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(36,1:5)  = (/ qpoint7(1) , qpoint7(3) , qpoint7(2) , qpoint7(4) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(37,1:5)  = (/ qpoint7(1) , qpoint7(4) , qpoint7(2) , qpoint7(3) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(38,1:5)  = (/ qpoint7(1) , qpoint7(3) , qpoint7(4) , qpoint7(2) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(39,1:5)  = (/ qpoint7(1) , qpoint7(4) , qpoint7(3) , qpoint7(2) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(40,1:5)  = (/ qpoint7(3) , qpoint7(1) , qpoint7(2) , qpoint7(4) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(41,1:5)  = (/ qpoint7(4) , qpoint7(1) , qpoint7(2) , qpoint7(3) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(42,1:5)  = (/ qpoint7(3) , qpoint7(1) , qpoint7(4) , qpoint7(2) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(43,1:5)  = (/ qpoint7(4) , qpoint7(1) , qpoint7(3) , qpoint7(2) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(44,1:5)  = (/ qpoint7(3) , qpoint7(4) , qpoint7(1) , qpoint7(2) , qpoint7(5) /)   
  cube_tet_rules(45)%rule(45,1:5)  = (/ qpoint7(4) , qpoint7(3) , qpoint7(1) , qpoint7(2) , qpoint7(5) /)   


END SUBROUTINE CUBE_TET_RULES_MAKE

END MODULE quad_tables





