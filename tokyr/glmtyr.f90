!------------------------------------------------------------------
!  GLMTYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  GLMT vsn adapted from TFTRYR.FOR
!
subroutine GLMTYR(NSHOT,ZYEAR)
  implicit none
  integer, intent(in) :: NSHOT
  character(len=*), intent(inout) :: ZYEAR
  !
  integer :: ilen

  zyear=' '
  ilen=len(zyear)
  ! from Gleb Kurskiev
  !  2074   2000      #1
  !  1180   2001      #2075
  !  1539   2002      #3256
  !  2498   2003      #4796
  !  4258   2004      #7295
  !  4274   2005      #11554
  !  2848   2006      #15829
  !  1845   2007      #18678
  !  2704   2008      #20524
  !  1907   2009      #23229
  !  2444   2010      #25137
  !  1705   2011      #27582
  !  1321   2012      #29288
  !  2135   2013      #30610
  !  2389   2014      #32746
  !  1036   2015      #35136
  !
  if(ilen.gt.3) then
    if(NSHOT.GE.23229) then
      ZYEAR(ilen-3:ilen)='2009'
    else if(NSHOT.GE.20524) then
      ZYEAR(ilen-3:ilen)='2008'
    else if(NSHOT.GE.18678) then
      ZYEAR(ilen-3:ilen)='2007'
    else if(NSHOT.GE.15829) then
      ZYEAR(ilen-3:ilen)='2006'
    else if(NSHOT.GE.11554) then
      ZYEAR(ilen-3:ilen)='2005'
    else if(NSHOT.GE.7295) then
      ZYEAR(ilen-3:ilen)='2004'
    else if(NSHOT.GE.4796) then
      ZYEAR(ilen-3:ilen)='2003'
    else if(NSHOT.GE.3256) then
      ZYEAR(ilen-3:ilen)='2002'
    else if(NSHOT.GE.2075) then
      ZYEAR(ilen-3:ilen)='2001'
    end if
  else
    if(NSHOT.GE.35136) then
      ZYEAR(ilen-1:ilen)='15'
    else if(NSHOT.GE.32746) then
      ZYEAR(ilen-1:ilen)='14'
    else if(NSHOT.GE.30610) then
      ZYEAR(ilen-1:ilen)='13'
    else if(NSHOT.GE.29288) then
      ZYEAR(ilen-1:ilen)='12'
    else if(NSHOT.GE.27582) then
      ZYEAR(ilen-1:ilen)='11'
    else if(NSHOT.GE.25137) then
      ZYEAR(ilen-1:ilen)='10'
    else
      ZYEAR(ilen-1:ilen)='09'
    end if
  end if
  !
  return
end subroutine GLMTYR

 
