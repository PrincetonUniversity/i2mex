!------------------------------------------------------------------
!  KSTRYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  KSTR vsn adapted from TFTRYR.FOR
!  04/15/13 C. Ludescher-Furth
!           shot info from Robert Budny
!
subroutine KSTRYR(NSHOT,ZYEAR)
  implicit none
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    if(NSHOT.GE.23181) then
      ZYEAR(ilen-3:ilen)='2020'
    else if(NSHOT.GE.21779) then
      ZYEAR(ilen-3:ilen)='2019'
    else if(NSHOT.GE.19455) then
      ZYEAR(ilen-3:ilen)='2018'
    else if(NSHOT.GE.17638) then
      ZYEAR(ilen-3:ilen)='2017'
    else if(NSHOT.GE.14576) then
      ZYEAR(ilen-3:ilen)='2016'
    else if(NSHOT.GE.11740) then
      ZYEAR(ilen-3:ilen)='2015'
    else if(NSHOT.GE.9470) then
      ZYEAR(ilen-3:ilen)='2014'
    else if(NSHOT.GE.8356) then
      ZYEAR(ilen-3:ilen)='2013'
    else if(NSHOT.GE.6471) then
      ZYEAR(ilen-3:ilen)='2012'
    else if(NSHOT.GE.4469) then
      ZYEAR(ilen-3:ilen)='2011'
    else if(NSHOT.GE.2372) then
      ZYEAR(ilen-3:ilen)='2010'
    else
      ZYEAR(ilen-3:ilen)='2009'
    end if
  else
    if(NSHOT.GE.23181) then
      ZYEAR(ilen-1:ilen)='20'
    else if(NSHOT.GE.21779) then
      ZYEAR(ilen-1:ilen)='19'
    else if(NSHOT.GE.19455) then
      ZYEAR(ilen-1:ilen)='18'
    else if(NSHOT.GE.17638) then
      ZYEAR(ilen-1:ilen)='17'
    else if(NSHOT.GE.14576) then
      ZYEAR(ilen-1:ilen)='16'
    else if(NSHOT.GE.11740) then
      ZYEAR(ilen-1:ilen)='15'
    else if(NSHOT.GE.9470) then
      ZYEAR(ilen-1:ilen)='14'
    else if(NSHOT.GE.8356) then
      ZYEAR(ilen-1:ilen)='13'
    else if(NSHOT.GE.6471) then
      ZYEAR(ilen-1:ilen)='12'
    else if(NSHOT.GE.4469) then
      ZYEAR(ilen-1:ilen)='11'
    else if(NSHOT.GE.6471) then
      ZYEAR(ilen-1:ilen)='10'
    else
      ZYEAR(ilen-1:ilen)='09'
    end if
  end if
  return
end subroutine KSTRYR
