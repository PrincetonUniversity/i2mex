!*DECK TR_DSORT
subroutine TR_DSORT (DX, DY, N, KFLAG)
  use iso_c_binding, only: fp => c_double
!***BEGIN PROLOGUE  TR_DSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      doUBLE PRECISION (SSORT-S, TR_DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   TR_DSORT sorts array DX and optionally makes the same interchanges in
!   array DY.  The array DX may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      DX - array of values to be sorted   (usually abscissas)
!      DY - array to be (optionally) carried along
!      N  - number of values in array DX to be sorted
!      KFLAG - control parameter
!            =  2  means sort DX in increasing order and carry DY along.
!            =  1  means sort DX in increasing order (ignoring DY)
!            = -1  means sort DX in decreasing order (ignoring DY)
!            = -2  means sort DX in decreasing order and carry DY along.
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
!   901012  Declared all variables; changed X,Y to DX,DY; changed
!           code to parallel SSORT. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           if-THEN-else-end if.  (RWC, WRB)
!***end PROLOGUE  TR_DSORT
!     .. Scalar Arguments ..
  integer :: KFLAG, N
!     .. Array Arguments ..
  doUBLE PRECISION DX(*), DY(*)
  !     .. Local Scalars ..
  doUBLE PRECISION R, T, TT, TTY, TY
  integer :: I, IJ, J, K, KK, L, M, NN
  !     .. Local Arrays ..
  integer :: IL(21), IU(21)
  !     .. External subroutines ..
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
  !***FIRST EXECUTABLE STATEMENT  TR_DSORT
  NN = N
  if (NN .LT. 1) THEN
     write(6,*) ' ?tr_dsort:  n.gt.0 expected, n = ',n
     return
  end if
!
  KK = ABS(KFLAG)
  if (KK.NE.1 .AND. KK.NE.2) THEN
     write(6,*) ' ?tr_dsort: kflag = ',kflag,' not in {-2,-1,1,2}'
     return
  end if
  !
  !     Alter array DX to get decreasing order if needed
  !
  if (KFLAG .LE. -1) THEN
     DX(1:NN) = -DX(1:NN)
  end if
!
  if (KK .EQ. 2) goto 100
!
!     Sort DX only
!
  M = 1
  I = 1
  J = NN
  R = 0.375_fp
!
20 if (I .EQ. J) goto 60
  if (R .LE. 0.5898437_fp) THEN
     R = R+3.90625e-2_fp
  else
     R = R-0.21875_fp
  end if
!
30 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = DX(IJ)
!
!     If first element of array is greater than T, interchange with T
!
  if (DX(I) .GT. T) THEN
     DX(IJ) = DX(I)
     DX(I) = T
     T = DX(IJ)
  end if
  L = J
!
!     If last element of array is less than than T, interchange with T
!
  if (DX(J) .LT. T) THEN
     DX(IJ) = DX(J)
     DX(J) = T
     T = DX(IJ)
     !
     !        If first element of array is greater than T, interchange with T
     !
     if (DX(I) .GT. T) THEN
        DX(IJ) = DX(I)
        DX(I) = T
        T = DX(IJ)
     end if
  end if
!
!     Find an element in the second half of the array which is smaller
!     than T
!
40 L = L-1
  if (DX(L) .GT. T) goto 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
50 K = K+1
  if (DX(K) .LT. T) goto 50
  !
!     Interchange these elements
!
  if (K .LE. L) THEN
     TT = DX(L)
     DX(L) = DX(K)
     DX(K) = TT
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
  T = DX(I+1)
  if (DX(I) .LE. T) goto 80
  K = I
  !
90 DX(K+1) = DX(K)
  K = K-1
  if (T .LT. DX(K)) goto 90
  DX(K+1) = T
  goto 80
  !
  !     Sort DX and carry DY along
  !
100 M = 1
  I = 1
  J = NN
  R = 0.375_fp
  !
110 if (I .EQ. J) goto 150
  if (R .LE. 0.5898437_fp) THEN
     R = R+3.90625e-2_fp
  else
     R = R-0.21875_fp
  end if
  !
120 K = I
  !
  !     Select a central element of the array and save it in location T
  !
  IJ = I + INT((J-I)*R)
  T = DX(IJ)
  TY = DY(IJ)
  !
  !     If first element of array is greater than T, interchange with T
  !
  if (DX(I) .GT. T) THEN
     DX(IJ) = DX(I)
     DX(I) = T
     T = DX(IJ)
     DY(IJ) = DY(I)
     DY(I) = TY
     TY = DY(IJ)
  end if
  L = J
  !
  !     If last element of array is less than T, interchange with T
  !
  if (DX(J) .LT. T) THEN
     DX(IJ) = DX(J)
     DX(J) = T
     T = DX(IJ)
     DY(IJ) = DY(J)
     DY(J) = TY
     TY = DY(IJ)
     !
     !        If first element of array is greater than T, interchange with T
     !
     if (DX(I) .GT. T) THEN
        DX(IJ) = DX(I)
        DX(I) = T
        T = DX(IJ)
        DY(IJ) = DY(I)
        DY(I) = TY
        TY = DY(IJ)
     end if
  end if
  !
  !     Find an element in the second half of the array which is smaller
  !     than T
  !
130 L = L-1
  if (DX(L) .GT. T) goto 130
  !
  !     Find an element in the first half of the array which is greater
  !     than T
  !
140 K = K+1
  if (DX(K) .LT. T) goto 140
  !
  !     Interchange these elements
  !
  if (K .LE. L) THEN
     TT = DX(L)
     DX(L) = DX(K)
     DX(K) = TT
     TTY = DY(L)
     DY(L) = DY(K)
     DY(K) = TTY
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
  T = DX(I+1)
  TY = DY(I+1)
  if (DX(I) .LE. T) goto 170
  K = I
  !
180 DX(K+1) = DX(K)
  DY(K+1) = DY(K)
  K = K-1
  if (T .LT. DX(K)) goto 180
  DX(K+1) = T
  DY(K+1) = TY
  goto 170
  !
  !     Clean up
  !
190 if (KFLAG .LE. -1) THEN
     DX(1:NN) = -DX(1:NN)
  end if
  return
end subroutine TR_DSORT
