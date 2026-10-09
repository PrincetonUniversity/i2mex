!------------------------------------------------------------------
!  ST80YR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
subroutine STARYR(NSHOT,ZYEAR)
  implicit none
  !
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen,ilen1,ilen3
  !
  zyear=' '
  ilen=len(zyear)
  ilen1 = max(1,ilen-1) ! RGA, quick fix for out of range error
  ilen3 = max(1,ilen-3)
  if(ilen.gt.3) then
     ZYEAR(ilen3:ilen)='2023'
  else
     ZYEAR(ilen1:ilen)='23'
  end if
  !
  return
end subroutine STARYR
