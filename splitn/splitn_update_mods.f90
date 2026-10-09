subroutine splitn_update_mod(itime,ztime,ierr)

  ! modify the time of update #itime
  ! return an error code if (a) block #itime does not exist or (b) the 
  !   modification would cause the update time sequence to be non-monotonic.

  use splitn_module
  implicit NONE

  integer, intent(in) :: itime   ! block index
  real*8, intent(in) :: ztime    ! new time desired for block
  integer, intent(out) :: ierr   ! return code, 0=OK

  real*8 :: zlim_low,zlim_high
  integer :: ii,jj,kk
  character*20 dstr
  character*32 znam32
  !-------------------------
  ierr = 0

  if((itime.le.0).or.(itime.gt.nupdate)) then
     write(6,*) ' ?splitn_update_mod: update block #',itime,' does not exist.'
     ierr = 1
     return
  endif

  ! check monotonicity

  if(itime.eq.1) then
     zlim_low = tinit
  else
     zlim_low = tup(itime-1)
  endif

  if(tinit.ge.ztime) then
     write(6,*) ' ?splitn_update_mod: update time: ',ztime,' seconds is'
     write(6,*) '  at/before namelist TINIT = ',tinit,' and so cannot be used.'
     ierr = 1
     return
  endif

  if(zlim_low.ge.ztime) then
     write(6,*) ' ?splitn_update_mod: update time: ',ztime,' seconds is'
     write(6,*) '  at/before the preceding update block time: ',zlim_low
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  zlim_high=tlarge
  if(itime.lt.nupdate) zlim_high=tup(itime+1)

  if(ztime.ge.zlim_high) then
     write(6,*) ' ?splitn_update_mod: update time: ',ztime,' seconds is'
     write(6,*) '  at/after the following update block time: ',zlim_high
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  ! OK...................

  znam32='~UPDATE_TIME'

  call splitn_put_chk('splitn_update_mod',znam32,jj,ierr)
  if(ierr.ne.0) return

  ! the update time value "belongs" to the preceding block...

  if(itime.le.1) then
     ii=varlist(jj)%addr   ! 1st time -> main body
     kk=varlist(jj)%nlinadr
  else
     ! nth time -> (n-1)'th block
     ii=nr8 + (itime-2)*nr8_st + varlist(jj)%addr_st
     kk=nreal+nint+nlog+nr8+nchv + &
          (itime-2)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  dbuf(ii)=ztime   ! reset here
  tup(itime)=ztime

  ! now reset the corresponding text

  dstr=' '
  write(dstr,'(1pd19.12)') ztime
  call splitn_fput_clean(dstr,'d')
  call splitn_put_sc_str(kk,znam32,dstr)

end subroutine splitn_update_mod

subroutine splitn_update_append_N(ztimes,ntimes,ierr)

  ! append N update blocks at the indicated times.
  ! requirements:
  !    ztimes(1).gt.[last existing update_time in file]
  !    ztimes in strict ascending order

  use splitn_module
  implicit NONE

  integer, intent(in) :: ntimes ! # of update times
  real*8, intent(in) :: ztimes(ntimes) ! update times in strict ascending order

  integer, intent(out) :: ierr   ! return code, 0=OK

  !-------------------------
  integer :: itime, ier2
  !-------------------------

  ierr = 0

  do itime=2,ntimes
     if(ztimes(itime).le.ztimes(itime-1)) then
        write(6,*) ' ?splitn_update_append_N: update times not ascending:'
        write(6,*) '    ztimes(',itime-1,') = ',ztimes(itime-1)
        write(6,*) '    ztimes(',itime,') = ',ztimes(itime)
        ierr = 1
        exit
     endif
  enddo

  if(ierr.ne.0) return

  do itime=1,ntimes

     call splitn_update_append(ztimes(itime),ier2)

  enddo

end subroutine splitn_update_append_N

subroutine splitn_update_append(ztime,ierr)

  ! append an update block at the indicated time.
  ! return an error if the addition would cause a non-monotonic time sequence.

  use splitn_module
  implicit NONE

  real*8, intent(in) :: ztime    ! time desired for new block
  integer, intent(out) :: ierr   ! return code, 0=OK

  character*20 dstr
  character*32 znam32
  integer :: iadr,jj,iadl
  !------------------------

  ierr = 0

  if(nupdate.gt.0) then
     if(ztime.le.tup(nupdate)) then
        write(6,*) ' ?splitn_update_append: the specified time: ',ztime
        write(6,*) '  is less than or equal to the last update time already'
        write(6,*) '  in the namelist: ',tup(nupdate)
        write(6,*) '  The times must be maintained in strict increasing order.'
        ierr=1
        return
     endif
  endif

  ! OK....................

  znam32 = '~UPDATE_TIME'

  call splitn_put_chk('splitn_update_mod',znam32,jj,ierr)
  if(ierr.ne.0) return

  if(nupdate.eq.0) then
     iadr = varlist(jj)%addr
     iadl = varlist(jj)%nlinadr
  else
     iadr = nr8 + (nupdate-1)*nr8_st + varlist(jj)%addr_st
     iadl=nreal+nint+nlog+nr8+nchv + &
          (nupdate-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + &
          varlist(jj)%nlinadr_st
  endif

  dbuf(iadr)=ztime

  ! now append the corresponding text

  kupdate=nupdate    ! internally ~update_time belongs in preceding block
  nupdate=nupdate+1
  nxblock=nxblock+1
  if(nxblock.gt.nupdate_max) call spxb_expand

  tup(nupdate)=ztime
 
  ! now append the corresponding text

  linup(nupdate)=nlines+2

  dstr=' '
  write(dstr,'(1pd19.12)') ztime
  call splitn_fput_clean(dstr,'d')

  call clear_fields
  newline = nlines + 1
  call add_nl_line(' ')

  newline = nlines + 1
  aline = ' ~update_time = '//trim(dstr)

  call parse_line(aline,.TRUE.,ierr)
  if(ierr.ne.0) return

  call add_nl_line(aline)

  call clear_fields
  newline = nlines + 1
  call add_nl_line(' ')

  kupdate = nupdate

  call splitn_put_edit_mark(0)

end subroutine splitn_update_append

subroutine splitn_update_remove(ztime,ierr)

  ! remove update blocks at/after the indicated time.
  ! report an error if no such blocks exist.

  use splitn_module
  implicit NONE

  real*8, intent(in) :: ztime    ! time desired for new block
  integer, intent(out) :: ierr   ! return code, 0=OK

  integer :: ii,jj,iupdate,istart

  !-------------------------
  ierr = 0

  if(nupdate.eq.0) then
     write(6,*) ' %splitn_update_remove: there are no updates to remove.'
     return
  endif

  if(ztime.gt.tup(nupdate)) then
     write(6,*) ' %splitn_update_remove: the removal time t=',ztime,' seconds'
     write(6,*) '  is after the last existing update time: ',tup(nupdate)
     write(6,*) '  so, no action is taken.'
     return
  endif

  do ii=nupdate,1,-1
     if(tup(ii).ge.ztime) iupdate=ii
  enddo

  ! OK....................
  ! mark corresponding text lines as removed...

  istart=linup(iupdate)
  do
     if(istart.eq.1) exit
     if(isblank(istart-1)) then
        istart = istart-1
     else
        exit
     endif
  enddo

  do ii=istart,nlines
     jj=ordl(ii)
     lenl(jj)=0  ! mark text line as deleted
  enddo

  ! make sure text references are zeroed out also

  istart = nreal+nint+nlog+nr8+nchv + &
       (iupdate-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + 1

  do ii=istart,maxlines(0)
     ilines(ii)=0
     krepeat(ii)=0
     ivrange(:,ii)=0
     irrange(:,ii)=0
  enddo

  tup(iupdate:)=tlarge
  linup(iupdate:)=ilarge

  nxblock=iupdate
  nupdate=iupdate - 1
  kupdate=nupdate

  ! mark

  call splitn_put_edit_mark(1)

end subroutine splitn_update_remove

subroutine splitn_update_insert(itime,ztime,ierr)

  ! before current update block (itime) insert a new update block at
  ! time ztime.
  ! return an error code if (a) block #itime does not exist or (b) the 
  !   modification would cause the update time sequence to be non-monotonic.

  ! method: use a temporary file...

  use splitn_module
  implicit NONE

  integer, intent(in) :: itime   ! block index, before which to insert.
  real*8, intent(in) :: ztime    ! time desired for new block to be inserted.
  integer, intent(out) :: ierr   ! return code, 0=OK

  real*8 :: zlim_low,zlim_high
  integer :: ilun                ! I/O channel for temporary file
  character*150 ztmpfil
  integer :: ilz,ii,jj,idum
  character*20 dstr

  character*50 zprog
  logical :: istarted

  !-------------------------
  ierr = 0

  if((itime.le.0).or.(itime.gt.nupdate)) then
     write(6,*) ' ?splitn_update_insert: block #',itime,' does not exist.'
     ierr = 1
     return
  endif

  ! check monotonicity

  if(itime.eq.1) then
     zlim_low = tinit
  else
     zlim_low = tup(itime-1)
  endif

  if(tinit.ge.ztime) then
     write(6,*) ' ?splitn_update_insert: insert time: ',ztime,' seconds is'
     write(6,*) '  at/before namelist TINIT = ',tinit,' and so cannot be used.'
     ierr = 1
     return
  endif

  if(zlim_low.ge.ztime) then
     write(6,*) ' ?splitn_update_insert: insert time: ',ztime,' seconds is'
     write(6,*) '  at/before the preceding update block time: ',zlim_low
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  zlim_high=tup(itime)

  if(ztime.ge.zlim_high) then
     write(6,*) ' ?splitn_update_insert: update time: ',ztime,' seconds is'
     write(6,*) '  at/after the insert target block time: ',zlim_high
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  ! OK...................

  dstr=' '
  write(dstr,'(1pd19.12)') ztime
  call splitn_fput_clean(dstr,'d')

  ! find temporary file

  call find_io_unit(ilun)
  call tmpfile_d('usplitn',ztmpfil,ilz)

  open(unit=ilun,file=ztmpfil(1:ilz),status='unknown',iostat=ierr)
  if(ierr.ne.0) then
     write(6,*) ' ?splitn_update_insert: file operation failed.'
     return
  endif

  do ii=1,nlines

     if(ii.eq.linup(itime)) then
        write(ilun,*) ' ' 
        write(ilun,'(" ~update_time = ",a,"  ! inserted by ",a)') &
             trim(dstr),edit_program
        write(ilun,*) ' ' 
     endif

     jj=ordl(ii)
     if(lenl(jj).ne.0) then
        write(ilun,'(A)') textnl(jj)(1:abs(lenl(jj)))
     endif
  enddo

  close(unit=ilun)

  istarted = edit_started
  zprog = edit_program

  call splitn_read(ztmpfil,ierr)
  if(ierr.ne.0) then
     write(6,*) ' ?? '//trim(edit_program)// &
          ' -- file readback error; cannot continue!'
     call bad_exit
  else
     call fdelete(ztmpfil,idum)
  endif

  edit_enabled=.TRUE.
  edit_program=zprog
  edit_started=istarted

  call splitn_put_edit_mark(0)  ! put in program tag...

end subroutine splitn_update_insert

subroutine splitn_update_insert_N(itime,ztimes,ntimes,ierr)

  ! before current update block (itime) insert a new set of update blocks
  ! times ztimes(1:ntimes).
  ! return an error code if (a) block #itime does not exist or (b) the 
  !   modification would cause the update time sequence to be non-monotonic.

  ! method: use a temporary file...

  use splitn_module
  implicit NONE

  integer, intent(in) :: itime   ! block index, before which to insert.
  integer, intent(in) :: ntimes  ! #of times
  real*8, intent(in) :: ztimes(ntimes)   ! times desired for new blocks 
  integer, intent(out) :: ierr   ! return code, 0=OK

  real*8 :: zlim_low,zlim_high
  integer :: ilun                ! I/O channel for temporary file
  character*150 ztmpfil
  integer :: ilz,ii,jj,idum,jtime
  character*20, dimension(:), allocatable :: dstr

  character*50 zprog
  logical :: istarted

  !-------------------------
  ierr = 0

  if((itime.le.0).or.(itime.gt.nupdate)) then
     write(6,*) ' ?splitn_update_insert: block #',itime,' does not exist.'
     ierr = 1
  endif

  ! check monotonicity of passed list of times
  do jtime=2,ntimes
     if(ztimes(jtime).le.ztimes(jtime-1)) then
        write(6,*) ' ?splitn_update_insert_N: update times not ascending:'
        write(6,*) '    ztimes(',jtime-1,') = ',ztimes(jtime-1)
        write(6,*) '    ztimes(',jtime,') = ',ztimes(jtime)
        ierr = 1
        exit
     endif
  enddo
  if(ierr.ne.0) return

  ! check monotonicity expected after insertion

  if(itime.eq.1) then
     zlim_low = tinit
  else
     zlim_low = tup(itime-1)
  endif

  if(tinit.ge.ztimes(1)) then
     write(6,*) ' ?splitn_update_insert: insert time: ',ztimes(1),' seconds is'
     write(6,*) '  at/before namelist TINIT = ',tinit,' and so cannot be used.'
     ierr = 1
     return
  endif

  if(zlim_low.ge.ztimes(1)) then
     write(6,*) ' ?splitn_update_insert: insert time: ',ztimes(1),' seconds is'
     write(6,*) '  at/before the preceding update block time: ',zlim_low
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  zlim_high=tup(itime)

  if(ztimes(ntimes).ge.zlim_high) then
     write(6,*) ' ?splitn_update_insert: update time: ',ztimes(ntimes), &
          ' seconds is'
     write(6,*) '  at/after the insert target block time: ',zlim_high
     write(6,*) '  and so cannot be used.'
     ierr = 1
     return
  endif

  ! OK...................

  allocate(dstr(ntimes))

  dstr=' '
  do jtime=1,ntimes
     write(dstr(jtime),'(1pd19.12)') ztimes(jtime)
     call splitn_fput_clean(dstr(jtime),'d')
  enddo

  ! find temporary file

  call find_io_unit(ilun)
  call tmpfile_d('usplitn',ztmpfil,ilz)

  open(unit=ilun,file=ztmpfil(1:ilz),status='unknown',iostat=ierr)
  if(ierr.ne.0) then
     write(6,*) ' ?splitn_update_insert: file operation failed.'
     return
  endif

  do ii=1,nlines

     if(ii.eq.linup(itime)) then
        write(ilun,*) ' ' 
        do jtime=1,ntimes
           write(ilun,'(" ~update_time = ",a,"  ! inserted by ",a)') &
                trim(dstr(jtime)),edit_program
           write(ilun,*) ' ' 
        enddo
     endif

     jj=ordl(ii)
     if(lenl(jj).ne.0) then
        write(ilun,'(A)') textnl(jj)(1:abs(lenl(jj)))
     endif
  enddo

  close(unit=ilun)

  istarted = edit_started
  zprog = edit_program

  call splitn_read(ztmpfil,ierr)
  if(ierr.ne.0) then
     write(6,*) ' ?? '//trim(edit_program)// &
          ' -- file readback error; cannot continue!'
     call bad_exit
  else
     call fdelete(ztmpfil,idum)
  endif

  edit_enabled=.TRUE.
  edit_program=zprog
  edit_started=istarted

  call splitn_put_edit_mark(0)  ! put in program tag...

end subroutine splitn_update_insert_N
