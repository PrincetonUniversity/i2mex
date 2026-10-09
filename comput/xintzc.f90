!******************** START FILE XINTZC.FOR ; GROUP IXCALC ******************
!-----------------------------------------------------------------
!  XINTZC
!
!  INTERPOLATE ZONE BDY FCN TO VALUES AT ZONE CENTER.
!   ONE SLIGHT EXTRAPOLATION AT CENTER-- ASSUME APPROACHING A ZERO
!   RATE OF CHANGE THERE


      subroutine XINTZC(FZB,FZC,N)
use iso_c_binding, only: fp => c_double


!	Updates:
!	tbt 03/08/94 - Changed Zfactor to .125 from .3333333

!	-----------------
      REAL FZC(N),FZB(N)
      Real Zfactor
!	-----------------

!       EXTRAPOLATION AT CENTER
      Zfactor = .125          ! Was .33333  - for parabolic fit.

      FZC(1)=FZB(1)- Zfactor*(FZB(2)-FZB(1))
!
      do 20 I=2,N
      IM1=I-1
      FZC(I)=0.5*(FZB(IM1)+FZB(I))
 20   continue
      return
      end
!******************** end FILE XINTZC.FOR ; GROUP IXCALC ******************
