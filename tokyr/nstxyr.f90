!& NSTXYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  NSTX vsn adapted from TFTRG1:TFTRYR.FOR
!  I don't have all the calender year boundary shot numbers!
!
subroutine NSTXYR(NSHOT,ZYEAR)
  implicit none
  !
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
    ZYEAR(ilen3:ilen)='2000'
    if(NSHOT.GE.104515) ZYEAR(ilen3:ilen)='2001'
    if(NSHOT.GE.106641) ZYEAR(ilen3:ilen)='2002'
    if(NSHOT.GE.109080) ZYEAR(ilen3:ilen)='2003'
    if(NSHOT.GE.110945) ZYEAR(ilen3:ilen)='2004'
    if(NSHOT.GE.115000) ZYEAR(ilen3:ilen)='2005'
    if(NSHOT.GE.119000) ZYEAR(ilen3:ilen)='2006'
    if(NSHOT.GE.122001) ZYEAR(ilen3:ilen)='2007'
    if(NSHOT.GE.125858) ZYEAR(ilen3:ilen)='2008'
    if(NSHOT.GE.130890) ZYEAR(ilen3:ilen)='2009'
    if(NSHOT.GE.136452) ZYEAR(ilen3:ilen)='2010'
    if(NSHOT.GE.142633) ZYEAR(ilen3:ilen)='2013'
  else
    ZYEAR(ilen1:ilen)='00'
    if(NSHOT.GE.104515) ZYEAR(ilen1:ilen)='01'
    if(NSHOT.GE.106641) ZYEAR(ilen1:ilen)='02'
    if(NSHOT.GE.109080) ZYEAR(ilen1:ilen)='03'
    if(NSHOT.GE.110945) ZYEAR(ilen1:ilen)='04'
    if(NSHOT.GE.115000) ZYEAR(ilen1:ilen)='05'
    if(nshot.ge.119000) ZYEAR(ilen1:ilen)='06'
    if(NSHOT.GE.122001) ZYEAR(ilen1:ilen)='07'
    if(NSHOT.GE.125858) ZYEAR(ilen1:ilen)='08'
    if(NSHOT.GE.130890) ZYEAR(ilen1:ilen)='09'
    if(NSHOT.GE.136452) ZYEAR(ilen1:ilen)='10'
    if(NSHOT.GE.142633) ZYEAR(ilen1:ilen)='13'
  end if
  !
  return
end subroutine NSTXYR
