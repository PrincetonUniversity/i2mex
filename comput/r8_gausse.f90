!******************** START FILE GAUSSE.FOR ; GROUP FRSUBS ******************
!-------------------------------------
!  THIS FILE WILL HOLD subroutineS OF FRANTIC NEUTRALS CODE
!  WHICH doN'T USE COMMON BLOCKS AND WHICH THEREFORE MAY BE
!  EXCLUDED FROM THE "OLYMPUSIZATION" PROCESS FOR NOW, ANYWAY.
!     D. MC CUNE  15 MAY 1981
!-------------------------------------
subroutine r8_gausse(A,IROW,N,M,X,EPS,IW)
  !	11/7/76
  !	THIS PROGRAM WAS WRITTEN BY HARRY H. TOWNER OF THE PRINCETON
  !	PLASMA PHYSICS LAB.  THIS ROUTINE WILL SOLVE A SYSTEM OF LINEAR
  !	EQUATIONS BY USING GAUSS ELIMATION AND BACKWORD SUBSITUTION.
  !	PARAMETER LIST:
  !	A - THE N BY M AUGMENTED MATRIX WITH THE RIGHT HAND VECTOR
  !		IN COLUMN M.
  !	IROW - THE ROW DIMENSION OF A AS DEFINED BY THE CALLING PROGRAM.
  !	N - THE # EQUATIONS.
  !	M - N+1.
  !	X - SOLUTION VECTOR WITH N ELEMENTS.
  !	EPS - EVERY PIVOT ELEMENT MUST BE GREATER THAN EPS.  if
  !		A PIVOT ELEMENT IS SMALLER THE EQUATIONS ARE
  !		DECLARED SINGULAR.
  !	IW - THE LOGICAL UNIT USED FOR PRINTING.
  !	NOTE THAT THE MATRIX A IS DESTROYED.
  !============
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, pi, twopi 
  implicit none
  integer :: irow,n,m,iw,nn,j,maxrow,i,jj,k,l1
  real(fp) :: eps,amax,absa,store,anorm,asub,sum
  real(fp), dimension(irow,m) :: a
  real(fp), dimension(n) :: x
  !============
  NN=N-1
  !	I FIRST REDUCE MATRIX A.
  do J=1,NN
    !	MAXROW WILL GIVE THE ROW IN WHICH THE MAXIMUM ELEMENT WAS FOUND.
    MAXROW=J
    AMAX=zero
    !	I NOW SEARCH FOR THE MAXIMUM ELEMENT IN ROW J AND GREATER.
    do I=J,N
      ABSA=ABS(A(I,J))
      if(ABSA-AMAX.gt.0) then 
         AMAX=ABSA
         MAXROW=I
      endif
    end do
    !	NOW THE MAXIMUM ELEMENT IN COL. J HAS THE VALUE OF AMAX AND
    !	IS IN ROW=MAXROW.
    !	I NOW CHECK TO SEE if EQUATIONS ARE SINGULAR.
    if(AMAX.GT.EPS) GO TO 275
    WRITE(IW,50) MAXROW,J
50  FORMAT(1X,'A PIVOT ELEMENT WAS FOUND TO BE LESS THAN EPS', &
         ' THE ELEMENT IS (',I3,',',I3,')')
    !  INSERT CALL TO TRANSP abort ROUTINE
    write(iw,*) ' [... error in GAUSSE ...]'
    call bad_exit
    !	if MAXROW doES NOT EQUAL J SWITCH ROWS.
275 if(MAXROW.EQ.J) GO TO 290
    do I=J,M
      STORE=A(J,I)
      A(J,I)=A(MAXROW,I)
      A(MAXROW,I)=STORE
    end do
    !	ROWS ARE NOW SWITCHED.
    !	I NOW NORMALIZE ROW J SO THAT A(J,J)=1.
290 ANORM=A(J,J)
    do I=J,M
      A(J,I)=A(J,I)/ANORM
    end do
    !	I NOW SUBTRACT ROW J FROM ROWS GREATER THEN J SUCH THAT
    !	THERE ARE ZEROS IN COLUMN J.
    JJ=J
300 JJ=JJ+1
    ASUB=A(JJ,J)
    do I=J,M
      A(JJ,I)=A(JJ,I)-A(J,I)*ASUB
    end do
    if(JJ.LT.N) GO TO 300
  end do
  !	I NOW PREFORM BACKWARD SUBSITION.
  X(N)=A(N,M)/A(N,N)
  do K=1,NN
    SUM=zero
    L1=N-K
    do I=1,K
      SUM=SUM+A(L1,L1+I)*X(L1+I)
    end do
    X(L1)=A(L1,M)-SUM
  end do
  RETURN
END subroutine r8_gausse
