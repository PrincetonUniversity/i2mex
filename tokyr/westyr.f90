!& WESTYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  WEST vsn adapted from TFTRG1:TFTRYR.FOR
!  I don't have all the calender year boundary shot numbers!
!
subroutine WESTYR(NSHOT,ZYEAR)
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
  ilen1 = max(1,ilen-1) ! RGA, quick fix for out of range error
  ilen3 = max(1,ilen-3)
  if(ilen.gt.3) then
     ZYEAR(ilen-3:ilen)='2019'
     if(NSHOT.GE.56052) ZYEAR(ilen3:ilen)='2020'
     if(NSHOT.GE.56544) ZYEAR(ilen3:ilen)='2021'
     if(NSHOT.GE.56921) ZYEAR(ilen3:ilen)='2022'
     if(NSHOT.GE.58001) ZYEAR(ilen3:ilen)='2023'
  else
     ZYEAR(ilen-3:ilen)='19'
     if(NSHOT.GE.56052) ZYEAR(ilen3:ilen)='20'
     if(NSHOT.GE.56544) ZYEAR(ilen3:ilen)='21'
     if(NSHOT.GE.56921) ZYEAR(ilen3:ilen)='22'
     if(NSHOT.GE.58001) ZYEAR(ilen3:ilen)='23'
  end if
  !
  return
end subroutine WESTYR
