!******************** START FILE XINTZ0.FOR ; GROUP IXCALC ******************
!-----------------------------------------------------------------
!  XINTZ0
!
!  INTERPOLATE ZONE BDY FCN TO VALUES AT ZONE CENTER.
!   Assume value is 0.0 at very center FZB(0) (if it existed)
!     tbt 6/1


      subroutine XINTZ0(FZB,FZC,N)
use iso_c_binding, only: fp => c_double


!	Updates:
!	tbt 06/01/95  At GA. Copied from XintzC

!	-----------------
      REAL FZC(N),FZB(N)
      Real Zfactor
!	-----------------

!       Interpolation AT CENTER - assume = 0
        FZB0 = 0.0

      FZC(1)= 0.5*(FZB0+FZB(1))
!
      do 20 I=2,N
      IM1=I-1
      FZC(I)=0.5*(FZB(IM1)+FZB(I))
 20   continue
      return
      end
!******************** end FILE XINTZ0.FOR ; GROUP IXCALC ******************
