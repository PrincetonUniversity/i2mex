!------------------------------------------------------------------
!  AUGDYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  AUGD vsn adapted from TFTRYR.FOR
!  10/12/04 C. Ludescher-Furth
!           shot info from G. Tardini (IPP)
!
subroutine AUGDYR(NSHOT,ZYEAR)
  implicit none
  !
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen
  !
  !  A CAREFUL EXAMINATION OF THE RECORD YIELDS...
  !
  zyear=' '
  ilen=len(zyear)
  !
  if(ilen.gt.3) then
    ZYEAR(ilen-3:ilen)='1999'
    if(NSHOT.GE.12980) ZYEAR(ilen-3:ilen)='2000'
    if(NSHOT.GE.14049) ZYEAR(ilen-3:ilen)='2001'
    if(NSHOT.GE.15011) ZYEAR(ilen-3:ilen)='2002'
    if(NSHOT.GE.16624) ZYEAR(ilen-3:ilen)='2003'
    if(NSHOT.GE.18244) ZYEAR(ilen-3:ilen)='2004'
    if(NSHOT.GE.19630) ZYEAR(ilen-3:ilen)='2005'
    if(NSHOT.GE.20713) ZYEAR(ilen-3:ilen)='2006'
    if(NSHOT.GE.21482) ZYEAR(ilen-3:ilen)='2007'
    if(NSHOT.GE.22586) ZYEAR(ilen-3:ilen)='2008'
    if(NSHOT.GE.24190) ZYEAR(ilen-3:ilen)='2009'
    if(NSHOT.GE.25890) ZYEAR(ilen-3:ilen)='2010'
    if(NSHOT.GE.26085) ZYEAR(ilen-3:ilen)='2011'
    if(NSHOT.GE.27405) ZYEAR(ilen-3:ilen)='2012'
    if(NSHOT.GE.29145) ZYEAR(ilen-3:ilen)='2013'
    if(NSHOT.GE.30151) ZYEAR(ilen-3:ilen)='2014'
    if(NSHOT.GE.31778) ZYEAR(ilen-3:ilen)='2015'
    if(NSHOT.GE.32805) ZYEAR(ilen-3:ilen)='2016'
    if(NSHOT.GE.33725) ZYEAR(ilen-3:ilen)='2017'
    if(NSHOT.GE.34997) ZYEAR(ilen-3:ilen)='2018'
    if(NSHOT.GE.35267) ZYEAR(ilen-3:ilen)='2019'
    if(NSHOT.GE.37001) ZYEAR(ilen-3:ilen)='2020'
    if(NSHOT.GE.38340) ZYEAR(ilen-3:ilen)='2021'
    if(NSHOT.GE.39965) ZYEAR(ilen-3:ilen)='2022'
  else
    ZYEAR(ilen-1:ilen)='99'
    if(NSHOT.GE.12980) ZYEAR(ilen-1:ilen)='00'
    if(NSHOT.GE.14049) ZYEAR(ilen-1:ilen)='01'
    if(NSHOT.GE.15011) ZYEAR(ilen-1:ilen)='02'
    if(NSHOT.GE.16624) ZYEAR(ilen-1:ilen)='03'
    if(NSHOT.GE.18244) ZYEAR(ilen-1:ilen)='04'
    if(NSHOT.GE.19630) ZYEAR(ilen-1:ilen)='05'
    if(NSHOT.GE.20713) ZYEAR(ilen-1:ilen)='06'
    if(NSHOT.GE.21482) ZYEAR(ilen-1:ilen)='07'
    if(NSHOT.GE.22586) ZYEAR(ilen-1:ilen)='08'
    if(NSHOT.GE.24190) ZYEAR(ilen-1:ilen)='09'
    if(NSHOT.GE.25890) ZYEAR(ilen-1:ilen)='10'
    if(NSHOT.GE.26085) ZYEAR(ilen-1:ilen)='11'
    if(NSHOT.GE.27405) ZYEAR(ilen-1:ilen)='12'
    if(NSHOT.GE.29145) ZYEAR(ilen-1:ilen)='13'
    if(NSHOT.GE.30151) ZYEAR(ilen-1:ilen)='14'
    if(NSHOT.GE.31778) ZYEAR(ilen-1:ilen)='15'
    if(NSHOT.GE.32805) ZYEAR(ilen-1:ilen)='16'
    if(NSHOT.GE.33725) ZYEAR(ilen-1:ilen)='17'
    if(NSHOT.GE.34997) ZYEAR(ilen-1:ilen)='18'
    if(NSHOT.GE.35267) ZYEAR(ilen-1:ilen)='19'
    if(NSHOT.GE.37001) ZYEAR(ilen-1:ilen)='20'
    if(NSHOT.GE.38340) ZYEAR(ilen-1:ilen)='21'
    if(NSHOT.GE.39965) ZYEAR(ilen-1:ilen)='22'
  end if
  !
  return
end subroutine AUGDYR
