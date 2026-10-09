!------------------------------------------------------------------
!
! RGA: LTX year from Craig Jacobson email 14Nov2013  
!
subroutine LTXYR(NSHOT,ZYEAR)
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
    if (nshot .lt. 108000) then
      ZYEAR(ilen-3:ilen)='2012'
    else if (nshot .le. 108105) then
      ZYEAR(ilen-3:ilen)='2013'
    else if (nshot .le. 110632) then
      ZYEAR(ilen-3:ilen)='2014'               
    else
      ZYEAR(ilen-3:ilen)='2015'
    end if
  else
    if (nshot .lt. 108000) then
      ZYEAR(ilen-1:ilen)='12'
    else if (nshot .le. 108105) then
      ZYEAR(ilen-1:ilen)='13'
    else if (nshot .le. 110632) then
      ZYEAR(ilen-1:ilen)='14'               
    else
      ZYEAR(ilen-1:ilen)='15'
    end if
  end if
  !
  return
end subroutine LTXYR
