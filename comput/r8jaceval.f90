subroutine r8jaceval(DR0DRO,DRMDRO,DYMDRO,RMSPL,YMSPL, &
     ZSNTHTK,ZCSTHTK,NMOM,ZJAC)
  !
  !  evaluate 2x2 jacobian using moments and derivatives
  !    ** symmetric formula **
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: nmom,imul
  real(fp) :: rmspl(nmom),ymspl(nmom)
  real(fp) :: dr0dro,drmdro(nmom),dymdro(nmom)
  real(fp) :: zcsthtk(nmom),zsnthtk(nmom)
  !
  real(fp) :: zjac(2,2)
  !
  !--------------------------------
  !
  ZJAC(1,1)=DR0DRO
  ZJAC(1,2)=0.0_fp
  !
  ZJAC(2,1)=0.0_fp
  ZJAC(2,2)=0.0_fp
  !
  do IMUL=1,NMOM
    ZJAC(1,1)=ZJAC(1,1)+DRMDRO(IMUL)*ZCSTHTK(IMUL)
    ZJAC(1,2)=ZJAC(1,2)-RMSPL(IMUL)*IMUL*ZSNTHTK(IMUL)
    !
    ZJAC(2,1)=ZJAC(2,1)+DYMDRO(IMUL)*ZSNTHTK(IMUL)
    ZJAC(2,2)=ZJAC(2,2)+YMSPL(IMUL)*IMUL*ZCSTHTK(IMUL)
  end do
  !
  return
end subroutine r8jaceval
