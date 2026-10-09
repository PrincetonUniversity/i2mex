!******************** START FILE WXFINT.FOR ; GROUP FILTR6 ******************
!
!   WXFINT
!
!  COMPUTES INTEGRAL OF W(T)*F(T)*DT
!   WHERE W(T)=1-ABS(XCEN-T)/DELTA IS THE WEIGHTING FUNCTION
!   AND F(T) IS THE PIECEWISE LINEAR INTERPOLATION FUCNTION
!   OF THE DATA SERIES.
!
!   INTEGRATE FROM X1 TO X2
!   (XF1,YF1),(XF2,YF2) ARE THE ENDPOINTS OF A LINEAR
!   SEGMENT OF F(T);  XCEN SHOULD NOT LIE BETWEEN XF1 AND
!   XF2, X1 AND X2 SHOULD EQUAL OR LIE
!   BETWEEN XF1 AND XF2
!
FUNCTION r8wxfint(X1,X2,XF1,YF1,XF2,YF2,XCEN,DELTA)
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: ZERO, HALF, ONE
  implicit none
  real(fp) :: r8wxfint
  ! Arguments
  real(fp), intent(in) :: x1, x2
  real(fp), intent(in) :: xf1, yf1, xf2, yf2
  real(fp), intent(in) :: xcen
  real(fp), intent(in) :: delta
  ! Local variables
  real(fp) :: zs, xa, xb, xfa, zc0, zdi, zc1, zc2
  integer :: isign
  !
  ZS=(YF2-YF1)/(XF2-XF1)
  XA=max((X1-XCEN),-DELTA)
  XB=min((X2-XCEN),DELTA)
  XFA=XF1-XCEN
  !
  ZC0=YF1-XFA*ZS
  ZDI=ONE/DELTA
  ISIGN=-1
  IF(XA.LT.ZERO) ISIGN=1
  ZC1=HALF*(ZS+ISIGN*ZC0*ZDI)
  ZC2=.3333333333333333333333333_fp*ISIGN*ZS*ZDI
  !
  r8wxfint=ZDI*(XB-XA)*(ZC0+(XB+XA)*ZC1+(XB*XB+XA*XB+XA*XA)*ZC2)
  !
  RETURN
END FUNCTION r8wxfint
