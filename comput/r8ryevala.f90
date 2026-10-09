subroutine r8ryevala(rmsplc,ymsplc,rmspls,ymspls, &
     zsnthtk,zcsthtk,nmom,zr99,zy99)
  !
  !  evaluate moments position:  asymmetric eq. moments formula
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: nmom, imul
  real(fp) :: zr99, zy99
  real(fp), dimension(0:nmom) :: rmsplc, ymsplc, rmspls, ymspls
  real(fp), dimension(nmom) :: zcsthtk, zsnthtk
  !
  !--------------------------------
  !
  ZR99 = RMSPLC(0)
  ZY99 = YMSPLC(0)
  do IMUL=1,NMOM
    ZR99 = ZR99 + RMSPLC(IMUL)*ZCSTHTK(IMUL) &
         + RMSPLS(IMUL)*ZSNTHTK(IMUL)
    ZY99 = ZY99 + YMSPLC(IMUL)*ZCSTHTK(IMUL) &
         + YMSPLS(IMUL)*ZSNTHTK(IMUL)
  end do
  !
  return
end subroutine r8ryevala

