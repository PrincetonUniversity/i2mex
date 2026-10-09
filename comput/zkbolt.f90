!******************** START FILE ZKBOLT.FOR ; GROUP KAPAI ******************
!
!
!@@@
!  BOLTONS KI(NU-HAT,EPS) FUNCTION (THESIS PP 74--78)
!
FUNCTION ZKBOLT(ZNUHAT,ZEPS)
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: ZKBOLT
  real(fp) :: zeps,znuhat
  real(fp) :: zse,zeps32,ztanc,ztbnc,ztcnc,ztdnc,znh,zkinc,ztaps
  real(fp) :: ztbps,ztcps,ztdps,zkips
  real(fp), dimension(3) :: zanc,zbnc,zcnc,zdnc ! N-C COEFFICIENTS FROM TABLE VI
  real(fp), dimension(3) :: zaps,zbps,zcps,zdps ! P-S COEFFICIENTS FROM TABLE VI
  !
  DATA ZANC/2.441_fp,-3.87_fp,2.19_fp/
  DATA ZBNC/.9362_fp,-3.109_fp,4.087_fp/
  DATA ZCNC/.241_fp,3.40_fp,-2.54_fp/
  DATA ZDNC/.2664_fp,-.352_fp,.44_fp/
  !
  DATA ZAPS/.364_fp,-2.76_fp,2.21_fp/
  DATA ZBPS/.553_fp,2.41_fp,-3.42_fp/
  DATA ZCPS/1.18_fp,.292_fp,1.07_fp/
  DATA ZDPS/.0188_fp,.180_fp,-.127_fp/
  !
  ZSE=SQRT(ZEPS)
  ZEPS32=ZSE*ZSE*ZSE
  ZTANC=.66_fp+((ZANC(3)*ZSE+ZANC(2))*ZSE+ZANC(1))*ZSE
  ZTBNC=((ZBNC(3)*ZEPS+ZBNC(2))*ZEPS+ZBNC(1))/(ZEPS**.75_fp)
  ZTCNC=((ZCNC(3)*ZEPS+ZCNC(2))*ZEPS+ZCNC(1))/(ZEPS32)
  ZTDNC=((ZDNC(3)*ZEPS+ZDNC(2))*ZEPS+ZDNC(1))/(ZEPS32)
  !
  !  N-C PART
  !
  ZNH=SQRT(ZNUHAT)
  ZKINC=ZTANC/(1._fp+ZTBNC*ZNH+ZTCNC*ZNUHAT+ZTDNC*ZNUHAT*ZNUHAT)
  !
  ZTAPS=((ZAPS(3)*ZSE+ZAPS(2))*ZSE+ZAPS(1))*ZEPS32
  ZTBPS=((ZBPS(3)*ZSE+ZBPS(2))*ZSE+ZBPS(1))*ZEPS32
  ZTCPS=((ZCPS(3)*ZSE+ZCPS(2))*ZSE+ZCPS(1))
  ZTDPS=((ZDPS(3)*ZSE+ZDPS(2))*ZSE+ZDPS(1))
  !
  !  P-S PART
  !
  ZKIPS=1.57_fp*ZEPS32 + (ZTAPS+ZTBPS*ZNH)/ &
       (1._fp+ZTCPS*ZNUHAT**1.5_fp+ZTDPS*ZNUHAT**2.5_fp)
  !
  ZKBOLT=ZKINC+ZKIPS
  !
  RETURN
END FUNCTION ZKBOLT
