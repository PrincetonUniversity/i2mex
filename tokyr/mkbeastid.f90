subroutine mkbeastid(runid,ipulse,cpulse)
  ! make BEAST pulse number from runid
  ! Try-Character of runid is replaced by integer
  ! A=1,..., Z=26
  implicit none
  character(len=*),  intent(in)    :: runid
  ! output
  integer,           intent(inout) :: ipulse  ! mds pulse number as integer
  character(len=10), intent(inout) :: cpulse  ! mds pulse number as charcater
  ! local
  integer         lr, ls, ios
  integer         runi
  character(len=2) :: runc

  call checkid(runid,lr,ls)
  if (lr==0) then
    cpulse = " " 
    ipulse = 0
    return
  end if

  if (lr .eq. ls ) then      !  old style runid
    cpulse = runid(1:lr)
    if (lr .eq. 5) then
      read(cpulse,'(i5)',iostat=ios) ipulse
    else
      read(cpulse,'(i4)',iostat=ios) ipulse
    endif
  else
    runi=ichar(runid(ls+1:ls+1))-64
    runc=char(runi/10+48) // char(mod(runi,10)+48)
    ! drop leading zero
    if (runid(1:1) .eq. '0') then
      cpulse = runid(2:ls)//runc//runid(ls+2:lr)
    else
      cpulse = runid(1:ls)//runc//runid(ls+2:lr)
    endif
    read(cpulse,'(i10)',iostat=ios) ipulse
  endif
 
  if (ios .ne. 0) ipulse=0
  return
end subroutine mkbeastid
