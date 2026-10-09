!------------------------------------------------------------------
!
!  ARIS vsn adapted from TFTRYR.FOR
!  06/08/12 C. Ludescher-Furth
!
subroutine ARISYR(NSHOT,ZYEAR)
  implicit none
  !
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  integer :: ilen
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    if(NSHOT.GE.10000) then
      ZYEAR(ilen-3:ilen)='2012'
    else
      ZYEAR(ilen-3:ilen)='2011'
    end if
  else
    if(NSHOT.GE.10000) then
      ZYEAR(ilen-1:ilen)='12'
    else
      ZYEAR(ilen-1:ilen)='11'
    end if
  end if
  !
  return
end subroutine ARISYR

 
