!& STEPYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  STEP vsn adapted from TFTRG1:TFTRYR.FOR
!  I don't have all the calender year boundary shot numbers!
!
subroutine STEPYR(NSHOT,ZYEAR)
  implicit none
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen,ilen1,ilen3
  !
  !  A CAREFUL EXAMINATION OF THE RECORD YIELDS...
  !
  zyear=' '
  ilen=len(zyear)
  !
  ilen1 = max(1,ilen-1) ! RGA, quick fix for out of range error
  ilen3 = max(1,ilen-3)
  if(ilen.gt.3) then
    ZYEAR(ilen3:ilen)='2018'
  else
    ZYEAR(ilen1:ilen)='18'
  end if
  !
  return
end subroutine STEPYR
 
