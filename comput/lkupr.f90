       FUNCTION LKUPR(X,TABLE,N)
!
! Look up a real number X in a TABLE of length N.
! returns the value lkupr such that:
!
!	table(lkupr) <= x < table(lkupr+1)
!
! assumes that table is in increasing order and does a binary search.
!
      integer :: LKUPR,n    
      REAL X,TABLE(N)
!
!========================================UPPER CASE & PRETTY PRINT RTM MAY 1988
!
      if(X .LT. TABLE(1))THEN
          LKUPR=0
          goto 70000		!doNE
      end if
!
      if(X .GE. TABLE(N))THEN
          LKUPR=N
          goto 70000		!doNE
      end if
!
      N1=1
      N2=N
!
100   continue	!LOOK AGAIN LOOP
!
      N3=(N2+N1)/2
!
      if(X .GE. TABLE(N3))THEN
          N1=N3
      else
          N2=N3
      end if
!
      if(N2 .EQ. N1+1)THEN
          LKUPR=N1
          goto 70000	    !doNE
      else
          goto 100	    !LOOK AGAIN
      end if

70000 continue	!All returns from here
!
      return
      end
