!& TFTRYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
subroutine TFTRYR(NSHOT,ZYEAR)
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
    ZYEAR(ilen-3:ilen)='1983'
    if(NSHOT.GE.8839) ZYEAR(ilen-3:ilen)='1984'
    if(NSHOT.GE.12000) ZYEAR(ilen-3:ilen)='1985'
    if(NSHOT.GE.17220) ZYEAR(ilen-3:ilen)='1986'
    if(NSHOT.GE.26699) ZYEAR(ilen-3:ilen)='1987'
    if(NSHOT.GE.31951) ZYEAR(ilen-3:ilen)='1988'
    if(NSHOT.GE.39117) ZYEAR(ilen-3:ilen)='1989'
    if(NSHOT.GE.45121) ZYEAR(ilen-3:ilen)='1990'
    if(NSHOT.GE.56486) ZYEAR(ilen-3:ilen)='1991'
    if(NSHOT.GE.60108) ZYEAR(ilen-3:ilen)='1992'
    if(NSHOT.GE.69346) ZYEAR(ilen-3:ilen)='1993'
    if(NSHOT.GE.73489) ZYEAR(ilen-3:ilen)='1994'
    if(NSHOT.GE.81729) ZYEAR(ilen-3:ilen)='1995'
    if(NSHOT.GE.90000) ZYEAR(ilen-3:ilen)='1996'
    if(nshot.ge.101951) zyear(ilen-3:ilen)='1997'
  else
    ZYEAR(ilen-1:ilen)='83'
    if(NSHOT.GE.8839) ZYEAR(ilen-1:ilen)='84'
    if(NSHOT.GE.12000) ZYEAR(ilen-1:ilen)='85'
    if(NSHOT.GE.17220) ZYEAR(ilen-1:ilen)='86'
    if(NSHOT.GE.26699) ZYEAR(ilen-1:ilen)='87'
    if(NSHOT.GE.31951) ZYEAR(ilen-1:ilen)='88'
    if(NSHOT.GE.39117) ZYEAR(ilen-1:ilen)='89'
    if(NSHOT.GE.45121) ZYEAR(ilen-1:ilen)='90'
    if(NSHOT.GE.56486) ZYEAR(ilen-1:ilen)='91'
    if(NSHOT.GE.60108) ZYEAR(ilen-1:ilen)='92'
    if(NSHOT.GE.69346) ZYEAR(ilen-1:ilen)='93'
    if(NSHOT.GE.73489) ZYEAR(ilen-1:ilen)='94'
    if(NSHOT.GE.81729) ZYEAR(ilen-1:ilen)='95'
    if(NSHOT.GE.90000) ZYEAR(ilen-1:ilen)='96'
    if(nshot.ge.101951) zyear(ilen-1:ilen)='97'
  end if
  !
  return
end subroutine TFTRYR
