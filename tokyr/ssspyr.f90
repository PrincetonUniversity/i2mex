!  SSSP vsn adapted from TFTRYR.FOR
subroutine SSSPYR(NSHOT,ZYEAR)
  implicit none
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  integer :: ilen
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    if(NSHOT.GE. 10000 .and. NSHOT.LE.19999) then
      ZYEAR(ilen-3:ilen)='2014'
    else
      ZYEAR(ilen-3:ilen)='2013'
    end if
  else
    if(NSHOT.GE.10000 .and. NSHOT.LE.19999) then
      ZYEAR(ilen-1:ilen)='14'
    else
      ZYEAR(ilen-1:ilen)='13'
    end if
  end if
  !
  return
end subroutine SSSPYR
