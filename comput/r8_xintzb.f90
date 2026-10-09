!
!  INTERPOLATE ZONE-CENTERED VARIABLE TO ZONE BOUNDARY VALUES
!
!  ONE SLIGHT EXTRAPOLATION AT OUTER EDGE
!
subroutine R8_XINTZB(FZC,FZB,N)
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: n,inm1,i,ip1
  real(fp), dimension(N) :: FZC, FZB

  INM1=N-1
  do I=1,INM1
    IP1=I+1
    FZB(I)=0.5_fp*(FZC(I)+FZC(IP1))
  end do

  ! Linear EXTRAPOLATION AT EDGE
  FZB(N) = FZC(N) + 0.5_fp*(FZC(N)-FZC(INM1))

  return
end subroutine R8_XINTZB
