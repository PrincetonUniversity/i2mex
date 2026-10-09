!------------------------------------------------------------------
!  EASTYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  EAST vsn adapted from TFTRYR.FOR
!  10/12/04 C. Ludescher-Furth
!           shot info from Siye Ding
!
subroutine EASTYR(NSHOT,ZYEAR)
  implicit none
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  integer :: ilen
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    if(NSHOT.GE.130818) then
      ZYEAR(ilen-3:ilen)='2024'
    else if(NSHOT.GE.119713) then
      ZYEAR(ilen-3:ilen)='2023'
    else if(NSHOT.GE.106985) then
      ZYEAR(ilen-3:ilen)='2022'
    else if(NSHOT.GE.96942) then
      ZYEAR(ilen-3:ilen)='2021'
    else if(NSHOT.GE.94342) then
      ZYEAR(ilen-3:ilen)='2020'
    else if(NSHOT.GE.82191) then
      ZYEAR(ilen-3:ilen)='2019'
    else if(NSHOT.GE.75809) then
      ZYEAR(ilen-3:ilen)='2018'
    else if(NSHOT.GE.71842) then
      ZYEAR(ilen-3:ilen)='2017'
    else if(NSHOT.GE.58226) then
      ZYEAR(ilen-3:ilen)='2016'
    else if(NSHOT.GE.53047) then
      ZYEAR(ilen-3:ilen)='2015'
    else if(NSHOT.GE.44327) then
      ZYEAR(ilen-3:ilen)='2014'
    else if(NSHOT.GE.36813) then
      ZYEAR(ilen-3:ilen)='2012'
    else if(NSHOT.GE.14705) then
      ZYEAR(ilen-3:ilen)='2010'
    else if(NSHOT.GE.9989) then
      ZYEAR(ilen-3:ilen)='2009'
    else if(NSHOT.GE.5672) then
      ZYEAR(ilen-3:ilen)='2008'
    else if(NSHOT.GE.3080) then
      ZYEAR(ilen-3:ilen)='2007'
    else if(NSHOT.GE.1138) then
      ZYEAR(ilen-3:ilen)='2006'
    else
      ZYEAR(ilen-3:ilen)='2005'
    end if
  else
    if(NSHOT.GE.130818) then
       ZYEAR(ilen-3:ilen)='24'
    else if(NSHOT.GE.119713) then
       ZYEAR(ilen-3:ilen)='23'
    else if(NSHOT.GE.106985) then
       ZYEAR(ilen-3:ilen)='22'
    else if(NSHOT.GE.96942) then
      ZYEAR(ilen-3:ilen)='21'
    else if(NSHOT.GE.94342) then
      ZYEAR(ilen-1:ilen)='20'
    else if(NSHOT.GE.82191) then
      ZYEAR(ilen-1:ilen)='19'
    else if(NSHOT.GE.75809) then
      ZYEAR(ilen-1:ilen)='18'
    else if(NSHOT.GE.71842) then
      ZYEAR(ilen-1:ilen)='17'
    else if(NSHOT.GE.58226) then
      ZYEAR(ilen-1:ilen)='16'
    else if(NSHOT.GE.53047) then
      ZYEAR(ilen-1:ilen)='15'
    else if(NSHOT.GE.44327) then
      ZYEAR(ilen-1:ilen)='14'
    else if(NSHOT.GE.36813) then
      ZYEAR(ilen-1:ilen)='12'
    else if(NSHOT.GE.14705) then
      ZYEAR(ilen-1:ilen)='10'
    else if(NSHOT.GE.9989) then
      ZYEAR(ilen-1:ilen)='09'
    else if(NSHOT.GE.5672) then
      ZYEAR(ilen-1:ilen)='08'
    else if(NSHOT.GE.3080) then
      ZYEAR(ilen-1:ilen)='07'
    else if(NSHOT.GE.1138) then
      ZYEAR(ilen-1:ilen)='06'
    else
      ZYEAR(ilen-1:ilen)='05'
    end if
  end if
  !
  return
end subroutine EASTYR
