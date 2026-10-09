function xintrp(ZT,TIM,DAT,NUM)

  use iso_c_binding, only: fp => c_double, c_int
  implicit NONE

  real(fp) :: xintrp
  integer, intent(in) :: num
  real(fp), intent(in), dimension(num) :: tim
  real(fp), intent(in), dimension(num) :: dat
  real(fp), intent(in) :: zt

  integer :: itm1,i,ip1,ip0
  integer, save :: iistart
 
  IF(ZT.LE.TIM(1)) THEN
     xintrp=DAT(1)
  ELSE IF(ZT.GE.TIM(NUM)) THEN
     xintrp=DAT(NUM)
  ELSE
     IF(IISTART.LE.0) IISTART=1
     IF(IISTART.GT.(num-1)) IISTART=1
     do I=1,num-1
        IP1=I+IISTART
        IF(IP1.GT.NUM) IP1=IP1-num+1
        IP0=IP1-1
        IF((ZT.GE.TIM(IP0)).AND.(ZT.LE.TIM(IP1))) exit 
     enddo
     xintrp=DAT(IP0)+(ZT-TIM(IP0))*(DAT(IP1)-DAT(IP0))/(TIM(IP1)-TIM(IP0))
     IISTART=IP0
  ENDIF

end function xintrp
