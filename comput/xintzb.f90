!******************** START FILE XINTZB.FOR ; GROUP IXCALC ******************
!-------------------------------------------------------------
!  XINTZB
!
!  INTERPOLATE ZONE-CENTERED VARIABLE TO ZONE BOUNDARY VALUES
!
!  ONE SLIGHT EXTRAPOLATION AT OUTER EDGE
!
      subroutine XINTZB(FZC,FZB,N)
use iso_c_binding, only: fp => c_double
!
      REAL FZC(N),FZB(N)
!
      INM1=N-1
      do 10 I=1,INM1
         IP1=I+1
         FZB(I)=0.5*(FZC(I)+FZC(IP1))
 10   continue

!       Linear EXTRAPOLATION AT EDGE
      FZB(N)=FZC(N)+0.5*(FZC(N)-FZC(INM1))
!
      return
      end
!******************** end FILE XINTZB.FOR ; GROUP IXCALC ******************
