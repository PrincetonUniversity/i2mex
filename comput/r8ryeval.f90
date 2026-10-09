subroutine r8ryeval(r0spl,rmspl,ymspl,zsnthtk,zcsthtk,nmom,zr99,zy99)
  !
  ! evaluate moments position:  symmetric eq. moments formula
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: nmom, imul
  real(fp) :: zr99, zy99
  real(fp) :: r0spl
  real(fp), dimension(nmom) :: rmspl, ymspl
  real(fp), dimension(nmom) :: zcsthtk, zsnthtk
  !
  !--------------------------------
  !
  ZR99=R0SPL
  ZY99 = 0.0_fp
  do IMUL=1,NMOM
    ZR99=ZR99+RMSPL(IMUL)*ZCSTHTK(IMUL)
    ZY99=ZY99+YMSPL(IMUL)*ZSNTHTK(IMUL)
  end do
  return
end subroutine r8ryeval
