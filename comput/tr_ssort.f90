!*DECK TR_SSORT
subroutine TR_SSORT (X, Y, N, KFLAG)
  use iso_c_binding, only: fp => c_double
!***BEGIN PROLOGUE  TR_SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (TR_SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   TR_SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           if-THEN-else-end if.  (RWC, WRB)
!***end PROLOGUE  TR_SSORT
!     .. Scalar Arguments ..
  integer :: KFLAG, N
!     .. Array Arguments ..
  REAL X(*), Y(*)
!     .. Local Scalars ..
  REAL R, T, TT, TTY, TY
  integer :: I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
  integer :: IL(21), IU(21)
!     .. External subroutines ..
!     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  TR_SSORT
  NN = N
  if (NN .LT. 1) THEN
     write(6,*) ' ?tr_ssort:  n.gt.0 expected, n = ',n
     return
  end if
!
  KK = ABS(KFLAG)
  if (KK.NE.1 .AND. KK.NE.2) THEN
     write(6,*) ' ?tr_ssort: kflag = ',kflag,' not in {-2,-1,1,2}'
     return
  end if
!
!     Alter array X to get decreasing order if needed
!
  if (KFLAG .LE. -1) THEN
     do  I=1,NN
        X(I) = -X(I)
     enddo
  end if
  !
  if (KK .EQ. 2) goto 100
!
!     Sort X only
!
  M = 1
  I = 1
  J = NN
  R = 0.375E0_fp
  !
20 if (I .EQ. J) goto 60
  if (R .LE. 0.5898437E0_fp) THEN
     R = R+3.90625E-2_fp
  else
     R = R-0.21875E0_fp
  end if
  !
30 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = X(IJ)
  !
  !     If first element of array is greater than T, interchange with T
  !
  if (X(I) .GT. T) THEN
     X(IJ) = X(I)
     X(I) = T
     T = X(IJ)
  end if
  L = J
  !
  !     If last element of array is less than than T, interchange with T
  !
  if (X(J) .LT. T) THEN
     X(IJ) = X(J)
     X(J) = T
     T = X(IJ)
     !
     !        If first element of array is greater than T, interchange with T
     !
     if (X(I) .GT. T) THEN
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
     end if
  end if
  !
  !     Find an element in the second half of the array which is smaller
  !     than T
  !
40 L = L-1
  if (X(L) .GT. T) goto 40
  !
  !     Find an element in the first half of the array which is greater
  !     than T
  !
50 K = K+1
  if (X(K) .LT. T) goto 50
  !
  !     Interchange these elements
  !
  if (K .LE. L) THEN
     TT = X(L)
     X(L) = X(K)
     X(K) = TT
     goto 40
  end if
  !
  !     Save upper and lower subscripts of the array yet to be sorted
  !
  if (L-I .GT. J-K) THEN
     IL(M) = I
     IU(M) = L
     I = K
     M = M+1
  else
     IL(M) = K
     IU(M) = J
     J = L
     M = M+1
  end if
  goto 70
  !
  !     Begin again on another portion of the unsorted array
  !
60 M = M-1
  if (M .EQ. 0) goto 190
  I = IL(M)
  J = IU(M)
  !
70 if (J-I .GE. 1) goto 30
  if (I .EQ. 1) goto 20
  I = I-1
  !
80 I = I+1
  if (I .EQ. J) goto 60
  T = X(I+1)
  if (X(I) .LE. T) goto 80
  K = I
  !
90 X(K+1) = X(K)
  K = K-1
  if (T .LT. X(K)) goto 90
  X(K+1) = T
  goto 80
  !
  !     Sort X and carry Y along
  !
  100 M = 1
  I = 1
  J = NN
  R = 0.375E0
  !
110 if (I .EQ. J) goto 150
  if (R .LE. 0.5898437E0_fp) THEN
     R = R+3.90625E-2_fp
  else
     R = R-0.21875E0_fp
  end if
  !
120 K = I
  !
  !     Select a central element of the array and save it in location T
  !
  IJ = I + INT((J-I)*R)
  T = X(IJ)
  TY = Y(IJ)
  !
  !     If first element of array is greater than T, interchange with T
  !
  if (X(I) .GT. T) THEN
     X(IJ) = X(I)
     X(I) = T
     T = X(IJ)
     Y(IJ) = Y(I)
     Y(I) = TY
     TY = Y(IJ)
  end if
  L = J
  !
  !     If last element of array is less than T, interchange with T
  !
  if (X(J) .LT. T) THEN
     X(IJ) = X(J)
     X(J) = T
     T = X(IJ)
     Y(IJ) = Y(J)
     Y(J) = TY
     TY = Y(IJ)
     !
     !        If first element of array is greater than T, interchange with T
     !
     if (X(I) .GT. T) THEN
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        Y(IJ) = Y(I)
        Y(I) = TY
        TY = Y(IJ)
     end if
  end if
  !
  !     Find an element in the second half of the array which is smaller
  !     than T
  !
130 L = L-1
  if (X(L) .GT. T) goto 130
  !
  !     Find an element in the first half of the array which is greater
  !     than T
  !
140 K = K+1
  if (X(K) .LT. T) goto 140
  !
  !     Interchange these elements
  !
  if (K .LE. L) THEN
     TT = X(L)
     X(L) = X(K)
     X(K) = TT
     TTY = Y(L)
     Y(L) = Y(K)
     Y(K) = TTY
     goto 130
  end if
  !
  !     Save upper and lower subscripts of the array yet to be sorted
  !
  if (L-I .GT. J-K) THEN
     IL(M) = I
     IU(M) = L
     I = K
     M = M+1
  else
     IL(M) = K
     IU(M) = J
     J = L
     M = M+1
  end if
  goto 160
  !
  !     Begin again on another portion of the unsorted array
  !
150 M = M-1
  if (M .EQ. 0) goto 190
  I = IL(M)
  J = IU(M)
  !
160 if (J-I .GE. 1) goto 120
  if (I .EQ. 1) goto 110
  I = I-1
  !
170 I = I+1
  if (I .EQ. J) goto 150
  T = X(I+1)
  TY = Y(I+1)
  if (X(I) .LE. T) goto 170
  K = I
  !
180 X(K+1) = X(K)
  Y(K+1) = Y(K)
  K = K-1
  if (T .LT. X(K)) goto 180
  X(K+1) = T
  Y(K+1) = TY
  goto 170
  !
  !     Clean up
  !
190 if (KFLAG .LE. -1) THEN
     do I=1,NN
        X(I) = -X(I)
     enddo
  endif
  return
end subroutine TR_SSORT
