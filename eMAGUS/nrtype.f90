MODULE nrtype
  IMPLICIT NONE
!*******************************************************************************
! This module contains definitions for various kind types and named constants
! It is an abbrevaited version of nrtype.f90 in "Numerical Recipes in Fortran 90", 2nd edn
! WH Press et al. Appendix C. Note that Appendix C is in the public domain.
!*******************************************************************************

  ! Symbolic names for kind types of 4-, 2- and 1-byte integers:
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
!!!
!!! This should be defined as the pointer size for a given architecture. eg
!!! for IA32 it should be a 32 bit int (I4B), for AMD64/SGI, it should be 
!!! a 64 bit int ( SELECTED_INT_KIND(RANGE(1)*2) seems to work)
!!!
  INTEGER, PARAMETER :: PTR = I4B
!!!
  ! Symbolic names for kind types of single- and double-precision reals:
  INTEGER, PARAMETER :: SP = 4
  INTEGER, PARAMETER :: DP = 8

  ! Symbolic names for kind types of single- and double-precision complex:
  INTEGER, PARAMETER :: SPC = 4
  INTEGER, PARAMETER :: DPC = 8

  ! Symbolic names for kind types of default logical:
  INTEGER, PARAMETER :: LGT = 4
 
  ! Frequently used mathematical constants:
  REAL(SP),     PARAMETER :: PI     = 3.14159265359
  REAL(SP),     PARAMETER :: EPS    = 1.0E-5           ! Small number for tests. 
  COMPLEX(SPC), PARAMETER :: j      = (0.0,1.0)        ! sqrt(-1)
  COMPLEX(SPC), PARAMETER :: ZERO_C = (0.0,0.0)        ! Zero.
  REAL(SP),     PARAMETER :: eps_0  = 8.8541878176E-12 ! Free space constants
  REAL(SP),     PARAMETER :: mu_0   = PI*4E-7          ! ditto.
  REAL(SP),     PARAMETER :: c_0    = 2.997925E+8      ! ditto.
  REAL(SP),     PARAMETER :: Z_zero = 376.73        ! ditto.
  REAL(SP),     PARAMETER :: D2R    = PI/180.0      ! Conversion factor degrees to radians.
  REAL(SP),     PARAMETER :: R2D    = 180.0/PI      ! Conversion factor degrees to radians.
 
!*******************************************************************************

END MODULE nrtype
