!------------------------------------------------------------------
!  HL2MYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  HL2M vsn adapted from TFTRYR.FOR
!  04/15/13 C. Ludescher-Furth
!
subroutine HL2MYR(NSHOT,ZYEAR)
  implicit none
  !
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    if(NSHOT.GE.1) then
      ZYEAR(ilen-3:ilen)='2015'
    else
      ZYEAR(ilen-3:ilen)='2014'
    end if
  else
    if(NSHOT.GE.1) then
      ZYEAR(ilen-1:ilen)='15'
    else
      ZYEAR(ilen-1:ilen)='14'
    end if
  end if
  !
  return
end subroutine HL2MYR
