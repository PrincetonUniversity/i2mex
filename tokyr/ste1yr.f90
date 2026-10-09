!------------------------------------------------------------------
!  STE1YR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
subroutine STE1YR(NSHOT,ZYEAR)
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
     ZYEAR(ilen3:ilen)='2025'
     if(NSHOT.GE.100000) ZYEAR(ilen3:ilen)='2026'
  else
     ZYEAR(ilen1:ilen)='25'
     if(NSHOT.GE.100000) ZYEAR(ilen1:ilen)='26'
  end if
  !
  return
end subroutine STE1YR
