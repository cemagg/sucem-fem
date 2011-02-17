MODULE math_tools
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Internal routines:
! (A routine is either listed as an interface, or in the private statement.)

INTERFACE CROSS_PRODUCT
  MODULE PROCEDURE CROSS_PRODUCT_REAL, CROSS_PRODUCT_COMPLEX
END INTERFACE

INTERFACE DET_DIM2
  MODULE PROCEDURE DET_DIM2_REAL, DET_DIM2_COMPLEX
END INTERFACE

INTERFACE SORT_VECTOR
  MODULE PROCEDURE SORT_VECTOR_INT, SORT_VECTOR_SINGLE
END INTERFACE


!PRIVATE
!*******************************************************************************

CONTAINS
!*******************************************************************************


FUNCTION CROSS_PRODUCT_COMPLEX(A,B) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! Functions to compute the cross product of two real-valued or 
! complex-valued vectors. 
! DBD 9 May 1999.
!******************************************************************************
   COMPLEX(SPC), INTENT(IN), DIMENSION (3) :: A,B
   COMPLEX(SPC), DIMENSION (3) :: CROSS_PRODUCT_COMPLEX
   CROSS_PRODUCT_COMPLEX(1) = A(2)*B(3) - A(3)*B(2)
   CROSS_PRODUCT_COMPLEX(2) = A(3)*B(1) - A(1)*B(3)
   CROSS_PRODUCT_COMPLEX(3) = A(1)*B(2) - A(2)*B(1)
END FUNCTION CROSS_PRODUCT_COMPLEX
!*******************************************************************************


FUNCTION CROSS_PRODUCT_REAL(A,B) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! Functions to compute the cross product of two real-valued or 
! complex-valued vectors. 
! DBD 9 May 1999.
!******************************************************************************
   REAL(SP), INTENT(IN), DIMENSION (3) :: A,B
   REAL(SP), DIMENSION (3) :: CROSS_PRODUCT_REAL
   CROSS_PRODUCT_REAL(1) = A(2)*B(3) - A(3)*B(2)
   CROSS_PRODUCT_REAL(2) = A(3)*B(1) - A(1)*B(3)
   CROSS_PRODUCT_REAL(3) = A(1)*B(2) - A(2)*B(1)
END FUNCTION CROSS_PRODUCT_REAL
!*******************************************************************************


FUNCTION DET_DIM2_REAL(A) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! This function computes the determinant of a matrix of dimension 2. DBD 1997
!*******************************************************************************
   REAL(SP), INTENT(IN), DIMENSION (2,2) :: A
   REAL(SP) DET_DIM2_REAL
   DET_DIM2_REAL = A(1,1) * A(2,2) - A(1,2) * A(2,1)  

END FUNCTION DET_DIM2_REAL
!*******************************************************************************

FUNCTION DET_DIM2_COMPLEX(A) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! This function computes the determinant of a matrix of dimension 2.
! with complex entries.  DBD 25 Aug 2004.
!*******************************************************************************
   COMPLEX(SPC), INTENT(IN), DIMENSION (2,2) :: A
   COMPLEX(SPC) DET_DIM2_COMPLEX
   DET_DIM2_COMPLEX = A(1,1) * A(2,2) - A(1,2) * A(2,1)  

END FUNCTION DET_DIM2_COMPLEX
!*******************************************************************************


FUNCTION DET_DIM3(A) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! This function computes the determinant of a matrix of dimension 3. DBD 1997
!*******************************************************************************
   REAL(SP), INTENT(IN), DIMENSION (3,3) :: A
   REAL(SP), DIMENSION(2,2):: MINOR_11, MINOR_12, MINOR_13
   REAL(SP) DET_DIM3

   MINOR_11 = A(2:3,2:3)
   MINOR_12(:,1) = A(2:3,1)  
   MINOR_12(:,2) = A(2:3,3)
   MINOR_13 = A(2:3,1:2)

   DET_DIM3 =   A(1,1)*DET_DIM2(MINOR_11) - A(1,2)*DET_DIM2(MINOR_12) & 
              + A(1,3)*DET_DIM2(MINOR_13)

END FUNCTION DET_DIM3
!*******************************************************************************


FUNCTION DET_DIM4(A) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! This function computes the determinant of a matrix of dimension 4. DBD 1997
!*******************************************************************************
   REAL(SP), INTENT(IN), DIMENSION (4,4) :: A
   REAL(SP), DIMENSION(3,3):: MINOR_11, MINOR_12, MINOR_13, MINOR_14
   REAL(SP) DET_DIM4 

   MINOR_11 = A(2:4,2:4)
   MINOR_12(:,1) = A(2:4,1)  
   MINOR_12(:,2:3) = A(2:4,3:4)
   MINOR_13(:,1:2) = A(2:4,1:2)  
   MINOR_13(:,3) = A(2:4,4)
   MINOR_14 = A(2:4,1:3)

   DET_DIM4 =   A(1,1)*DET_DIM3(MINOR_11) - A(1,2)*DET_DIM3(MINOR_12) & 
              + A(1,3)*DET_DIM3(MINOR_13) - A(1,4)*DET_DIM3(MINOR_14)

END FUNCTION DET_DIM4
!*******************************************************************************


SUBROUTINE FIND_IN_LIST(list,listsize,element,listpos)
  USE nrtype
  IMPLICIT NONE
!*************************************************************************
! This routine recieves integer <list> of length <listsize>. It searches
! for <element>, and returns its position in <listpos>. <listpos=0> 
! indicates that the element is not in the list.
! 9 May 99 00:25 RHG
! 17 Nov 14h45 DBD. USE statement changed to conform to MIPS 7.2.1
!*************************************************************************
  INTEGER(I4B),INTENT(IN) :: listsize
  INTEGER(I4B), DIMENSION(listsize),INTENT(IN) :: list
  INTEGER(I4B), INTENT(IN)  :: element
  INTEGER(I4B), INTENT(OUT) :: listpos
  INTEGER(I4B) :: search 

  listpos = 0  ! This is the default value
               ! and indicates that element was not found in list
                  
  DO search = 1,listsize
    IF (list(search).EQ.element) THEN
      listpos = search
      EXIT
    END IF
  END DO
 
END SUBROUTINE FIND_IN_LIST
!*******************************************************************************


FUNCTION PHASE(field)
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! A function to return the phase of a complex number 
! The value lies in the interval (-180,+180] degrees.
! DBD June 99, revised 17 Feb 2000. 
!*******************************************************************************
   COMPLEX(SPC), INTENT(IN):: field
   REAL(SP) PHASE
   IF (ABS(REAL(field)).LT.EPS.AND.ABS(AIMAG(field)).LT.EPS) THEN
     PHASE = 0.0_SP
   ELSE
     PHASE = 180.0_SP/PI * ATAN2(AIMAG(field),REAL(field))
   END IF
END FUNCTION PHASE
!*******************************************************************************


SUBROUTINE SORT_VECTOR_INT(VECTOR,N)
  USE nrtype
  IMPLICIT NONE
!******************************************************************************
! Sorting subroutines: sorts the vector VECTOR of length N 
! into ascending order. Presently implemented for integer and real SP
! vectors. Operator overloading done in INTERFACE block.
! The sort algorithm used is Insertion Sort; see "Introduction to 
! Algorithms", Cormen et al., MIT Press, 1996, pp2-4.
! Author: D B Davidson
! 3 June 1999.
!******************************************************************************
  INTEGER(I4B) , INTENT (IN) :: N
  INTEGER(I4B), INTENT (INOUT), DIMENSION(N)  :: VECTOR
  INTEGER(I4B) i, jj, key
  INTEGER(I4B), DIMENSION(N) :: A  ! temporary variable
  A = VECTOR
  ! Now perform insertion sort on temporary variable
  DO jj=2,SIZE(A)
     key = A(jj)  
     ! Insert A(jj) into the sorted sequence A(1..jj-1)
     i = jj-1
     DO WHILE (i.GT.0.AND.A(i).GT.key)
       A(i+1) = A(i)
       i = i-1 
     END DO 
     A(i+1) = key
  END DO ! Sorting finished.
  VECTOR = A ! Return sorted vector
  RETURN
END SUBROUTINE SORT_VECTOR_INT
!******************************************************************************


SUBROUTINE SORT_VECTOR_SINGLE(VECTOR,N)
  USE nrtype
  IMPLICIT NONE
!******************************************************************************
! Sorting subroutines: sorts the vector VECTOR of length N 
! into ascending order. Presently implemented for integer and real SP
! vectors. Operator overloading done in INTERFACE block.
! The sort algorithm used is Insertion Sort; see "Introduction to 
! Algorithms", Cormen et al., MIT Press, 1996, pp2-4.
! Author: D B Davidson
! 3 June 1999.
!******************************************************************************
  INTEGER(I4B) , INTENT (IN) :: N
  REAL(SP), INTENT (INOUT), DIMENSION(N)  :: VECTOR
  INTEGER(I4B) i, jj
  REAL(SP) key
  REAL(SP), DIMENSION(N) :: A  ! temporary variable
  A = VECTOR
  ! Now perform insertion sort on temporary variable
  DO jj=2,SIZE(A)
     key = A(jj)  
     ! Insert A(jj) into the sorted sequence A(1..jj-1)
     i = jj-1
     DO WHILE (i.GT.0.AND.A(i).GT.key)
       A(i+1) = A(i)
       i = i-1 
     END DO 
     A(i+1) = key
  END DO ! Sorting finished.
  VECTOR = A ! Return sorted vector
  RETURN
END SUBROUTINE SORT_VECTOR_SINGLE
!******************************************************************************


FUNCTION TIME_DIFFERENCE(time1,time2)
  USE nrtype
  IMPLICIT NONE
!*******************************************************************************
! Calculates the difference, in seconds (real), between two times, <time1-time2>,
! as returned by DATE_AND_TIME. This routine does not take leap years into account.
! MMB 25 Jan 2001
!*******************************************************************************
  INTEGER(I4B), DIMENSION(8), INTENT(IN) :: time1,time2
  REAL(SP) :: TIME_DIFFERENCE

  INTEGER(I4B) :: year1,year2,month1,month2,day1,day2,hour1,hour2,min1,min2, &
                  sec1,sec2,milsec1,milsec2
  INTEGER(I4B), DIMENSION(12) :: month_data
  INTEGER(I4B) :: count   ! counter

  month_data = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! First create integer variables containing the two years, months, days, 
  ! hours, minutes, seconds and milliseconds:
  year1   = time1(1)
  month1  = time1(2)
  day1    = time1(3)
  hour1   = time1(5)
  min1    = time1(6)
  sec1    = time1(7)
  milsec1 = time1(8)
  year2   = time2(1)
  month2  = time2(2)
  day2    = time2(3)
  hour2   = time2(5)
  min2    = time2(6)
  sec2    = time2(7)
  milsec2 = time2(8)

  ! Now subtract year1 from both year counts and the calculate the number of seconds
  ! for both times (integer, ignoring milsec):
  year2 = year2 - year1
  year1 = year1 - year1

  ! Converts years and months to days:  
  day2 = day2 + 365*year2
  DO count = 1,12
    IF (month1.GE.count) day1 = day1 + month_data(count)
    IF (month2.GE.count) day2 = day2 + month_data(count)
  END DO

  ! Converts days, hours, and minutes to seconds:
  sec1 = 86400*day1 + 3600*hour1 + 60*min1 + sec1
  sec2 = 86400*day2 + 3600*hour2 + 60*min2 + sec2

  ! Calculate the difference:
  IF (milsec2.GE.milsec1) THEN
    TIME_DIFFERENCE = REAL(sec2-sec1) + 0.001*REAL(milsec2-milsec1)
  ELSE
    TIME_DIFFERENCE = REAL(sec2-sec1-1) + 0.001*REAL(1000+milsec2-milsec1)
  END IF

END FUNCTION TIME_DIFFERENCE
!*******************************************************************************


FUNCTION VECTOR_LENGTH(A) 
   USE nrtype
   IMPLICIT NONE
!*******************************************************************************
! Function to compute the Euclidean length of a 3D geometrical vector 
! DBD 20 March 2001.
!******************************************************************************
   REAL(SP), INTENT(IN), DIMENSION (3) :: A
   REAL(SP) VECTOR_LENGTH
   REAL(SP) TEMP
   TEMP = DOT_PRODUCT(A,A)
   IF (TEMP.LT.0.0_SP) THEN
     STOP 'IE: In function VECTOR_LENGTH.'
   END IF
   VECTOR_LENGTH  = SQRT(TEMP)
END FUNCTION VECTOR_LENGTH
!******************************************************************************


END MODULE math_tools
