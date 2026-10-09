!******************** START FILE IRNS.FOR ; GROUP RNG ******************
      FUNCTION IRNS(IDUMMY)
!  11/29/78
!  THIS FUNCTION WAS WRITTEN BY HARRY H. TOWNER OF THE PRINCETON PLASMA
!  PHYSICS LAB.  IRNS WILL return A UNIQUE integer :: SUITABLE FOR SETTING
!  THE RANdoM NUMBER SEED.
!  5/01 RGA -- replace with f90 standard call
      integer :: iv(8)
      call date_and_time(values=iv)
      secnds = 60*(60*iv(5)+iv(6))+iv(7)+iv(8)/1000.

      IRNS=1.E4*SECNDS
      IRNS = 2*IRNS + 1	!MAKE SURE IT'S ODD, RTM 26 FEB 86
      return
      end
!******************** end FILE IRNS.FOR ; GROUP RNG ******************
