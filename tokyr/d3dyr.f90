!& D3DYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  D3D vsn adapted from TFTRG1:TFTRYR.FOR
!
subroutine D3DYR(NSHOT,ZYEAR)
  implicit none
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
    ZYEAR(ilen3:ilen)='1983'
    if(NSHOT.GE.40560) ZYEAR(ilen3:ilen)='1984'
    if(NSHOT.GE.50531) ZYEAR(ilen3:ilen)='1985'
    if(NSHOT.GE.50690) ZYEAR(ilen3:ilen)='1986'
    if(NSHOT.GE.53583) ZYEAR(ilen3:ilen)='1987'
    if(NSHOT.GE.57504) ZYEAR(ilen3:ilen)='1988'
    if(NSHOT.GE.62639) ZYEAR(ilen3:ilen)='1989'
    if(NSHOT.GE.68163) ZYEAR(ilen3:ilen)='1990'
    if(NSHOT.GE.70801) ZYEAR(ilen3:ilen)='1991'
    if(NSHOT.GE.74182) ZYEAR(ilen3:ilen)='1992'
    if(NSHOT.GE.76295) ZYEAR(ilen3:ilen)='1993'
    if(NSHOT.GE.80293) ZYEAR(ilen3:ilen)='1994'
    if(NSHOT.GE.84468) ZYEAR(ilen3:ilen)='1995'
    if(NSHOT.GE.88061) ZYEAR(ilen3:ilen)='1996'
    if(NSHOT.GE.90895) ZYEAR(ilen3:ilen)='1997'
    if(NSHOT.GE.94458) ZYEAR(ilen3:ilen)='1998'
    if(NSHOT.GE.97674) ZYEAR(ilen3:ilen)='1999'
    if(NSHOT.GE.100513) ZYEAR(ilen3:ilen)='2000'
    if(NSHOT.GE.104793) ZYEAR(ilen3:ilen)='2001'
    if(NSHOT.GE.108779) ZYEAR(ilen3:ilen)='2002'
    if(NSHOT.GE.111982) ZYEAR(ilen3:ilen)='2003'
    if(NSHOT.GE.116638) ZYEAR(ilen3:ilen)='2004'
    if(NSHOT.GE.121238) ZYEAR(ilen3:ilen)='2005'
    if(NSHOT.GE.124000) ZYEAR(ilen3:ilen)='2006'
    if(NSHOT.GE.127146) ZYEAR(ilen3:ilen)='2007'
    if(NSHOT.GE.131027) ZYEAR(ilen3:ilen)='2008'
    if(NSHOT.GE.134884) ZYEAR(ilen3:ilen)='2009'
    if(NSHOT.GE.140817) ZYEAR(ilen3:ilen)='2010'
    if(NSHOT.GE.143160) ZYEAR(ilen3:ilen)='2011'
    if(NSHOT.GE.147780) ZYEAR(ilen3:ilen)='2012'
    if(NSHOT.GE.151236) ZYEAR(ilen3:ilen)='2013'
    if(NSHOT.GE.155933) ZYEAR(ilen3:ilen)='2014'
    if(NSHOT.GE.160923) ZYEAR(ilen3:ilen)='2015'
    if(NSHOT.GE.164750) ZYEAR(ilen3:ilen)='2016'
    if(NSHOT.GE.168439) ZYEAR(ilen3:ilen)='2017'
    if(NSHOT.GE.174574) ZYEAR(ilen3:ilen)='2018'
    if(NSHOT.GE.177828) ZYEAR(ilen3:ilen)='2019'
    if(NSHOT.GE.181663) ZYEAR(ilen3:ilen)='2020'
  else
    ZYEAR(ilen1:ilen)='83'
    if(NSHOT.GE.40560) ZYEAR(ilen1:ilen)='84'
    if(NSHOT.GE.50531) ZYEAR(ilen1:ilen)='85'
    if(NSHOT.GE.50690) ZYEAR(ilen1:ilen)='86'
    if(NSHOT.GE.53583) ZYEAR(ilen1:ilen)='87'
    if(NSHOT.GE.57504) ZYEAR(ilen1:ilen)='88'
    if(NSHOT.GE.62639) ZYEAR(ilen1:ilen)='89'
    if(NSHOT.GE.68163) ZYEAR(ilen1:ilen)='90'
    if(NSHOT.GE.70801) ZYEAR(ilen1:ilen)='91'
    if(NSHOT.GE.74182) ZYEAR(ilen1:ilen)='92'
    if(NSHOT.GE.76295) ZYEAR(ilen1:ilen)='93'
    if(NSHOT.GE.80293) ZYEAR(ilen1:ilen)='94'
    if(NSHOT.GE.84468) ZYEAR(ilen1:ilen)='95'
    if(NSHOT.GE.88061) ZYEAR(ilen1:ilen)='96'
    if(NSHOT.GE.90895) ZYEAR(ilen1:ilen)='97'
    if(NSHOT.GE.94458) ZYEAR(ilen1:ilen)='98'
    if(NSHOT.GE.97674) ZYEAR(ilen1:ilen)='99'
    if(NSHOT.GE.100513) ZYEAR(ilen1:ilen)='00'
    if(NSHOT.GE.104793) ZYEAR(ilen1:ilen)='01'
    if(NSHOT.GE.108779) ZYEAR(ilen1:ilen)='02'
    if(NSHOT.GE.111982) ZYEAR(ilen1:ilen)='03'
    if(NSHOT.GE.116638) ZYEAR(ilen1:ilen)='04'
    if(NSHOT.GE.121238) ZYEAR(ilen1:ilen)='05'
    if(NSHOT.GE.124000) ZYEAR(ilen1:ilen)='06'
    if(NSHOT.GE.127146) ZYEAR(ilen1:ilen)='07'
    if(NSHOT.GE.131027) ZYEAR(ilen1:ilen)='08'
    if(NSHOT.GE.134884) ZYEAR(ilen1:ilen)='09'
    if(NSHOT.GE.140817) ZYEAR(ilen1:ilen)='10'
    if(NSHOT.GE.143160) ZYEAR(ilen1:ilen)='11'
    if(NSHOT.GE.147780) ZYEAR(ilen1:ilen)='12'
    if(NSHOT.GE.151236) ZYEAR(ilen1:ilen)='13'
    if(NSHOT.GE.155933) ZYEAR(ilen1:ilen)='14'
    if(NSHOT.GE.160923) ZYEAR(ilen1:ilen)='15'
    if(NSHOT.GE.164750) ZYEAR(ilen1:ilen)='16'
    if(NSHOT.GE.168439) ZYEAR(ilen1:ilen)='17'
    if(NSHOT.GE.174574) ZYEAR(ilen1:ilen)='18'
    if(NSHOT.GE.177828) ZYEAR(ilen1:ilen)='19'
    if(NSHOT.GE.181663) ZYEAR(ilen1:ilen)='20'
  end if
  !
  return
end subroutine D3DYR
