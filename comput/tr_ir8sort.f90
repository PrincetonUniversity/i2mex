!DECK TR_IR8SORT
subroutine TR_IR8SORT (IDX, DY, N, KFLAG)
  !***BEGIN PROLOGUE  TR_IR8SORT
  !***PURPOSE  Sort an array and optionally make the same interchanges in
  !            an auxiliary array.  The array may be sorted in increasing
  !            or decreasing order.  A slightly modified QUICKSORT
  !            algorithm is used.
  !
  !  >> modified version (dmc) -- sort by integer key; carry real(fp) DY array
  !
  !***LIBRARY   SLATEC
  !***CATEGORY  N6A2B
  !***TYPE      doUBLE PRECISION (SSORT-S, TR_DSORT-D, ISORT-I)
  !***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
  !***AUTHOR  Jones, R. E., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***DESCRIPTION
  !
  !   TR_IR8SORT sorts array IDX and optionally makes the same interchanges in
  !   array DY.  The array IDX may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      IDX - array of values to be sorted   (usually abscissas)
  !      DY - array to be (optionally) carried along
  !      N  - number of values in array IDX to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort IDX in increasing order and carry DY along.
  !            =  1  means sort IDX in increasing order (ignoring DY)
  !            = -1  means sort IDX in decreasing order (ignoring DY)
  !            = -2  means sort IDX in decreasing order and carry DY along.
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
  !***END PROLOGUE  TR_IR8SORT
  !     .. Scalar Arguments ..
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  !
  integer KFLAG, N
  !     .. Array Arguments ..
  integer IDX(*)
  real(fp) DY(*)
  !     .. Local Scalars ..
  real(fp) R, TTY, TY
  integer T, TT
  integer I, IJ, J, K, KK, L, M, NN
  !     .. Local Arrays ..
  integer IL(21), IU(21)
  !     .. External Subroutines ..
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
  !***FIRST EXECUTABLE STATEMENT  TR_IR8SORT
  NN = N
  if (NN .LT. 1) THEN
    write(6,*) ' ?tr_ir8sort:  n.gt.0 expected, n = ',n
    RETURN
  end if
  !
  KK = ABS(KFLAG)
  if (KK.NE.1 .AND. KK.NE.2) THEN
    write(6,*) ' ?tr_ir8sort: kflag = ',kflag,' not in {-2,-1,1,2}'
    RETURN
  end if
  !
  !     Alter array IDX to get decreasing order if needed
  !
  if (KFLAG .LE. -1) THEN
    do I=1,NN
      IDX(I) = -IDX(I)
    end do
  end if
  !
  if (KK .EQ. 2) GO TO 100
  !
  !     Sort IDX only
  !
  M = 1
  I = 1
  J = NN
  R = 0.375D0
  !
20 if (I .EQ. J) GO TO 60
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
  T = IDX(IJ)
  !
  !     If first element of array is greater than T, interchange with T
  !
  if (IDX(I) .GT. T) THEN
    IDX(IJ) = IDX(I)
    IDX(I) = T
    T = IDX(IJ)
  end if
  L = J
  !
  !     If last element of array is less than than T, interchange with T
  !
  if (IDX(J) .LT. T) THEN
    IDX(IJ) = IDX(J)
    IDX(J) = T
    T = IDX(IJ)
    !
    !        If first element of array is greater than T, interchange with T
    !
    if (IDX(I) .GT. T) THEN
      IDX(IJ) = IDX(I)
      IDX(I) = T
      T = IDX(IJ)
    end if
  end if
  !
  !     Find an element in the second half of the array which is smaller
  !     than T
  !
40 L = L-1
  if (IDX(L) .GT. T) GO TO 40
  !
  !     Find an element in the first half of the array which is greater
  !     than T
  !
50 K = K+1
  if (IDX(K) .LT. T) GO TO 50
  !
  !     Interchange these elements
  !
  if (K .LE. L) THEN
    TT = IDX(L)
    IDX(L) = IDX(K)
    IDX(K) = TT
    GO TO 40
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
  GO TO 70
  !
  !     Begin again on another portion of the unsorted array
  !
60 M = M-1
  if (M .EQ. 0) GO TO 190
  I = IL(M)
  J = IU(M)
  !
70 if (J-I .GE. 1) GO TO 30
  if (I .EQ. 1) GO TO 20
  I = I-1
  !
80 I = I+1
  if (I .EQ. J) GO TO 60
  T = IDX(I+1)
  if (IDX(I) .LE. T) GO TO 80
  K = I
  !
90 IDX(K+1) = IDX(K)
  K = K-1
  if (T .LT. IDX(K)) GO TO 90
  IDX(K+1) = T
  GO TO 80
  !
  !     Sort IDX and carry DY along
  !
100 M = 1
  I = 1
  J = NN
  R = 0.375D0
  !
110 if (I .EQ. J) GO TO 150
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
  T = IDX(IJ)
  TY = DY(IJ)
  !
  !     If first element of array is greater than T, interchange with T
  !
  if (IDX(I) .GT. T) THEN
    IDX(IJ) = IDX(I)
    IDX(I) = T
    T = IDX(IJ)
    DY(IJ) = DY(I)
    DY(I) = TY
    TY = DY(IJ)
  end if
  L = J
  !
  !     If last element of array is less than T, interchange with T
  !
  if (IDX(J) .LT. T) THEN
    IDX(IJ) = IDX(J)
    IDX(J) = T
    T = IDX(IJ)
    DY(IJ) = DY(J)
    DY(J) = TY
    TY = DY(IJ)
    !
    !        If first element of array is greater than T, interchange with T
    !
    if (IDX(I) .GT. T) THEN
      IDX(IJ) = IDX(I)
      IDX(I) = T
      T = IDX(IJ)
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
  if (IDX(L) .GT. T) GO TO 130
  !
  !     Find an element in the first half of the array which is greater
  !     than T
  !
140 K = K+1
  if (IDX(K) .LT. T) GO TO 140
  !
  !     Interchange these elements
  !
  if (K .LE. L) THEN
    TT = IDX(L)
    IDX(L) = IDX(K)
    IDX(K) = TT
    TTY = DY(L)
    DY(L) = DY(K)
    DY(K) = TTY
    GO TO 130
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
  GO TO 160
  !
  !     Begin again on another portion of the unsorted array
  !
150 M = M-1
  if (M .EQ. 0) GO TO 190
  I = IL(M)
  J = IU(M)
  !
160 if (J-I .GE. 1) GO TO 120
  if (I .EQ. 1) GO TO 110
  I = I-1
  !
170 I = I+1
  if (I .EQ. J) GO TO 150
  T = IDX(I+1)
  TY = DY(I+1)
  if (IDX(I) .LE. T) GO TO 170
  K = I
  !
180 IDX(K+1) = IDX(K)
  DY(K+1) = DY(K)
  K = K-1
  if (T .LT. IDX(K)) GO TO 180
  IDX(K+1) = T
  DY(K+1) = TY
  GO TO 170
  !
  !     Clean up
  !
190 if (KFLAG .LE. -1) THEN
    do I=1,NN
      IDX(I) = -IDX(I)
    end do
  end if
  RETURN
end subroutine TR_IR8SORT
