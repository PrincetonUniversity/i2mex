subroutine r8jacevala(DRMDROC,DYMDROC,DRMDROS,DYMDROS, &
     RMSPLC,YMSPLC,RMSPLS,YMSPLS,ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
  !
  !  evaluate 2x2 jacobian using moments and derivatives
  !    ** asymmetric formula **
  !
  !  evaluate moments position:  asymmetric eq. moments formula
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer nmom,imul
  real(fp) :: zc,zs
  real(fp) :: rmsplc(0:nmom),ymsplc(0:nmom),rmspls(0:nmom),ymspls(0:nmom)
  real(fp) :: drmdroc(0:nmom),dymdroc(0:nmom)
  real(fp) :: drmdros(0:nmom),dymdros(0:nmom)
  real(fp) :: zcsthtk(nmom),zsnthtk(nmom)
  real(fp) :: zjac(2,2)
  !
  ZJAC(1,1)= DRMDROC(0)
  ZJAC(1,2)= 0.0_fp
  !
  ZJAC(2,1)= DYMDROC(0)
  ZJAC(2,2)= 0.0_fp
  !
  do IMUL=1,NMOM
    ZC = ZCSTHTK(IMUL)
    ZS = ZSNTHTK(IMUL)
    !
    ZJAC(1,1) = ZJAC(1,1) + DRMDROC(IMUL)*ZC + DRMDROS(IMUL)*ZS
    ZJAC(1,2) = ZJAC(1,2) - RMSPLC(IMUL)*IMUL*ZS + RMSPLS(IMUL)*IMUL*ZC
    !
    ZJAC(2,1) = ZJAC(2,1) + DYMDROC(IMUL)*ZC + DYMDROS(IMUL)*ZS
    ZJAC(2,2) = ZJAC(2,2) - YMSPLC(IMUL)*IMUL*ZS + YMSPLS(IMUL)*IMUL*ZC
  end do
  !
  return
end subroutine r8jacevala
