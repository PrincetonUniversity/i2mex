!------------------------------------------------------------------
!  ITERYR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
!  July 27, 2009 CLF: ITER vsn adapted from D3DYR.FOR
!
subroutine iteryr(nshot, zyear)
  !
  implicit none
  !
  integer, intent(in) :: nshot
  character(len=*), intent(out) :: zyear
  !
  integer :: ilen, ilen1, ilen3
  !
  zyear=' '
  ilen=len(zyear)
  !
  ilen1 = max(1,ilen-1) 
  ilen3 = max(1,ilen-3)
  if(ilen.gt.3) then
    zyear(ilen3:ilen)='2013'
    if (nshot .ge. 31000 .and. nshot .le. 39999) then
      zyear(ilen3:ilen)='2011'
    else if (nshot .ge. 50000 .and. nshot .lt. 51000) then
      zyear(ilen3:ilen)='2013'
    else if (nshot .eq. 11000 .or. nshot .eq. 31000) then
      zyear(ilen3:ilen)='2009'
    else if (nshot .le. 2000) then
      zyear(ilen3:ilen)='1995'
    else if (nshot .eq. 3000) then
      zyear(ilen3:ilen)='2007'
    else if (nshot .eq. 97048) then
      zyear(ilen3:ilen)='1995'
    else if (nshot .eq. 10000) then
      zyear(ilen3:ilen)='1991' 
    else if (nshot .eq. 10001) then 
      zyear(ilen3:ilen)='1995'
    else if (nshot .ge. 20000 .and. nshot .lt. 20150) then
      zyear(ilen3:ilen)='2006'
    else if (nshot .ge. 20150 .and. nshot .lt. 20200) then
      zyear(ilen3:ilen)='2009'
    else if (nshot .ge. 20200 .and. nshot .lt. 20300) then
      zyear(ilen3:ilen)='2006'
    else if (nshot .ge. 20300 .and. nshot .lt. 20400) then
      zyear(ilen3:ilen)='2007'
    else if (nshot .ge. 20400 .and. nshot .le. 20800) then
      zyear(ilen3:ilen)='2008'
    else if (nshot .eq. 20900 ) then
      zyear(ilen3:ilen)='2006'
    else if (nshot .eq. 28000 ) then
      zyear(ilen3:ilen)='2002'
    else if (nshot .eq. 30000 ) then
      zyear(ilen3:ilen)='2004'
    else if (nshot .ge. 40100 .and. nshot .le. 40600) then
      zyear(ilen3:ilen)='2005'
    else if (nshot .eq. 40000 ) then
      zyear(ilen3:ilen)='2004'
    else if (nshot .ge. 40700 .and. nshot .le. 40900) then
      zyear(ilen3:ilen)='2006'
    else if (nshot .eq. 60000 ) then
      zyear(ilen3:ilen)='2006'
    else if (nshot .eq. 70000 ) then
      zyear(ilen3:ilen)='2009'
    else if (nshot .ge. 80000 ) then
      zyear(ilen3:ilen)='2008'
    end if
  else
    zyear(ilen1:ilen)='13'
    if (nshot .ge. 31000 .and. nshot .le. 39999) then
      zyear(ilen1:ilen)='11'
    else if (nshot .ge. 50000 .and. nshot .lt. 51000) then
      zyear(ilen1:ilen)='13'
    else if (nshot .eq. 11000 .or. nshot .eq. 31000) then
      zyear(ilen1:ilen)='09'
    else if (nshot .le. 2000) then
      zyear(ilen1:ilen)='95'
    else if (nshot .eq. 3000) then
      zyear(ilen1:ilen)='07'
    else if (nshot .eq. 97048) then
      zyear(ilen1:ilen)='95'
    else if (nshot .eq. 10000) then
      zyear(ilen1:ilen)='91' 
    else if (nshot .eq. 10001) then 
      zyear(ilen1:ilen)='95'
    else if (nshot .ge. 20000 .and. nshot .lt. 20150) then
      zyear(ilen3:ilen)='06'
    else if (nshot .ge. 20150 .and. nshot .lt. 20200) then
      zyear(ilen3:ilen)='09'
    else if (nshot .ge. 20200 .and. nshot .lt. 20300) then
      zyear(ilen3:ilen)='06'
    else if (nshot .ge. 20300 .and. nshot .lt. 20400) then
      zyear(ilen3:ilen)='07'
    else if (nshot .ge. 20400 .and. nshot .le. 20800) then
      zyear(ilen3:ilen)='08'
    else if (nshot .eq. 20900) then
      zyear(ilen1:ilen)='06'
    else if (nshot .eq. 28000) then
      zyear(ilen1:ilen)='02'
    else if (nshot .eq. 30000) then
      zyear(ilen1:ilen)='04'
    else if (nshot .ge. 40100 .and. nshot .le. 40600) then
      zyear(ilen1:ilen)='05'
    else if (nshot .eq. 40000) then
      zyear(ilen1:ilen)='04'
    else if (nshot .ge. 40700 .and. nshot .le. 40900) then
      zyear(ilen1:ilen)='06'
    else if (nshot .eq. 60000) then
      zyear(ilen1:ilen)='06'
    else if (nshot .eq. 70000) then
      zyear(ilen1:ilen)='09'
    else if (nshot .ge. 80000) then
      zyear(ilen1:ilen)='08'
    end if
  end if

  return
end subroutine iteryr
