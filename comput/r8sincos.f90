subroutine r8sincos(ZTHETA,JANG,SNTHTK,CSTHTK)
  !
  !	THIS subroutine CALCULATES:
  !		SIN(ZTHETA)
  !		COS(ZTHETA)
  !		SIN(2*ZTHETA)
  !		COS(2*ZTHETA)
  !	ETC., UP TO JANG*ZTHETA
  !
  !
  !		ZTHETA IS INPUT ANGLE
  !		JANG IS HIGHEST N*ZTHETA TO CALCULATE
  !		SNTHTK IS THE ARRAY CONTAINING SIN(N*ZTHETA)
  !		CSTHTK IS THE ARRAY CONTAINING COS(N*ZTHETA)
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: jang, i
  real(fp) :: ztheta, zthetap, zsinp, zcosp
  real(fp), dimension(jang) :: snthtk, csthtk
  !
  !  LOCAL MEMORY (DMC 6 JUL 1994)
  !
  DATA ZTHETAP/0.0E0_fp/
  DATA ZSINP/0.0E0_fp/
  DATA ZCOSP/1.0E0_fp/
  !
  SAVE ZTHETAP,ZSINP,ZCOSP
  !
  !--------------------------------------------------------------------
  !
  !  DMC -- USE LOCAL MEMORY FOR SPEED
  !
  if(ZTHETA.NE.ZTHETAP) then
    !
    !  EVALUATE SIN,COS
    !
    SNTHTK(1)=SIN(ZTHETA)
    CSTHTK(1)=COS(ZTHETA)
    ZTHETAP=ZTHETA
    ZSINP=SNTHTK(1)
    ZCOSP=CSTHTK(1)
  else
    !
    !  REUSE PREVIOUS RESULTS
    !
    SNTHTK(1)=ZSINP
    CSTHTK(1)=ZCOSP
  end if
  !
  do I=2,JANG
    SNTHTK(I)=SNTHTK(I-1)*CSTHTK(1)+CSTHTK(I-1)*SNTHTK(1)
    CSTHTK(I)=CSTHTK(I-1)*CSTHTK(1)-SNTHTK(I-1)*SNTHTK(1)
  end do

  return
end subroutine r8sincos
