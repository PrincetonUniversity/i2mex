subroutine splitn_seek_update(tmin,filename,inum_found,inum_skipped,ierr)

  ! Look for a file with update blocks for an existing namelist; apply the
  ! updates if found.

  ! Read the set of update blocks -- screen out those before time tmin (sec).
  ! If there is a single block with no time label, assign it the time tmin.
  ! If there are multiple blocks, all must have associated time labels
  !   " ~update_time = <time-value> ";
  ! a leading block with no time label will be taken to be at time tmin, if
  ! there are no subsequent blocks with t < tmin that are screened out.

  ! The ~update_time recorded times in the file must be in ascending order.

  ! NOTE: if no file is found this is considered normal: no updates found.

  use splitn_module
  implicit NONE

  !----------------------------------------
  !  arguments:

  real*8, intent(in) :: tmin   ! minimum time -- screen times before tmin.
  ! if there are data records unassociated with an ~update_time, tmin is
  ! assigned to these records.

  character*(*), intent(in) :: filename   ! file containing updates

  integer, intent(out) :: inum_found      ! number of update blocks found

  integer, intent(out) :: inum_skipped    ! number of update blocks skipped
  !                                       ! (t < tmin); will also be set to
  !                                       ! one (1) if file exists but
  !                                       ! inum_found=0.

  integer, intent(out) :: ierr            ! exit status (0=OK)

  !  an error status could occur (a) if a temporary file I/O operation fails
  !  or (b) if the file contains updates to non-updatable namelist quantities,
  !  or (c) if there is a syntax error or times are in non-ascending order in
  !  the input file.

  !----------------------------------------
  !  local:
  integer :: ilun1,ilun2,io_list(2),ios
  character*20 z20
  character*150 tmp_filename
  integer :: irec,irec0,irec_t0,ic,ilen,ilenv,iexcl,itilda,ieqs
  integer :: icmt,irec_d0,iblock
  logical :: tmin_leader
  integer :: in1,in2,iv1,ilz
  real*8 :: ztime,ztime0,ztimep
  real*8 :: zdtmin = 1.0d-5  ! seconds
  !----------------------------------------

  inum_found = 0
  inum_skipped = 0
  ierr = 0

  if(.not.have_database) then
     write(6,*) ' ?splitn_seek_update: a base namelist must be read prior'
     write(6,*) '  to calling this routine.'
     ierr = 1
     return
  endif

  call find_io_unit_list(2,io_list)
  ilun1 = io_list(1)
  ilun2 = io_list(2)

  open(unit=ilun1,file=filename,status='old',iostat=ios)
  if(ios.ne.0) then
#ifdef __DEBUG
     write(6,*) ' %splitn_seek_update: update file not found: ',trim(filename)
     write(6,*) '  it is assumed that no updates are intended.'
#endif
     return     ! normally, this is done silently: no file found, no update.
  endif

  !  file found: scan for update times
  irec = 0      ! record count
  irec_t0 = 0   ! record count of first qualified ~time_update record.
  irec_d0 = 0   ! record count of first qualified data record, if (a) no prior
                ! ~time_update record is found and (b) no subsequent 
                ! ~time_update record is screened out for t < tmin.

  tmin_leader = .FALSE. ! set .TRUE. if leading records are to be assigned tmin

  ztime0 = -tlarge  ! time of first qualified ~update_time; -infinity initially
  ztimep = -tlarge

  do
     irec = irec + 1
     read(ilun1,'(A)',iostat=ios) aline
     if(ios.ne.0) exit  ! assume EOF

     bline = aline      ! save for messages (aline will be manipulated).

     itilda=0
     icmt=0
     ilen=len(trim(aline))
     do ic=1,ilen
        if(.NOT.iswhite(aline(ic:ic))) then
           if(aline(ic:ic).eq.'~') then
              ! looks like an ~update_time
              itilda=ic
              exit
           else
              if(aline(ic:ic).eq.'!') then
                 icmt=ic
              else
                 if((irec_d0.eq.0).and.(inum_skipped.eq.0)) then
                    ! this records presence of non-comment data prior to first
                    ! ~update_time record
                    if(irec_t0.eq.0) irec_d0=irec ! data, not starting with "~"
                 endif
              endif
              exit
           endif
        endif
     enddo

     if(itilda.eq.0) cycle

     ! look for update time  " ~update_time = <time-value>   ! cmt "
     ! screen out comment if it exists

     ieqs = index(aline,'=')
     if(ieqs.le.0) then
        ierr=1
        write(6,*) ' ?splitn_seek_update: "=" assignment syntax not found: '
        write(6,*) '  '//trim(bline)
        exit
     endif

     iexcl = index(aline,'!')
     if(iexcl.gt.0) then
        aline(iexcl:)=' '
        ilen=len(trim(aline))
     endif

     call uupper(aline)

     in1=0  ! name indices
     in2=0

     iv1=0  ! value index (ilen is the 2nd index, as comment field was removed)

     do ic=1,ilen
        if(in1.eq.0) then
           if(.NOT.iswhite(aline(ic:ic))) then
              in1=ic
           endif
           cycle
        endif
        if(in2.eq.0) then
           if(iswhite(aline(ic:ic)).or.(ic.eq.ieqs)) then
              in2=ic-1
           endif
           cycle
        endif

        if(ic.le.ieqs) cycle

        if(.NOT.iswhite(aline(ic:ic))) then
           iv1=ic
           exit
        endif
     enddo

     if((in1.eq.ieqs).or.(in2.eq.0).or.(iv1.eq.0)) then
        ierr=1
        write(6,*) ' ?splitn_seek_update: "=" assignment syntax parse error: '
        write(6,*) '  '//trim(bline)
        exit
     endif

     if(aline(in1:in2).ne.'~UPDATE_TIME') then
        ierr=1
        write(6,*) ' ?splitn_seek_update: "~" variable must be "~update_time":'
        write(6,*) '  '//trim(bline)
        exit
     endif

     ilenv=ilen-iv1+1
     if(ilenv.gt.20) then
        ierr=1
        write(6,*) ' ?splitn_seek_update: "~update_time" value field too long:'
        write(6,*) '  '//trim(bline)
        exit
     endif

     z20=' '
     z20(20-ilenv+1:20) = aline(iv1:ilen)

     read(z20,'(G20.0)',iostat=ierr) ztime
     if(ierr.ne.0) then
        write(6,*) ' ?splitn_seek_update: "~update_time" value field decode failed:'
        write(6,*) '  '//trim(bline)
        exit
     endif

     if(ztime.lt.tmin) then
        inum_skipped = inum_skipped + 1
        if(inum_skipped.eq.1) write(6,*) '  tmin = ',tmin
        write(6,*) ' %splitn_seek_update: skipped update block (t < tmin): ', &
             ztime
        irec_d0 = 0 ! any leading data (prior to skipped time) is ignored;
        irec_t0 = 0 ! still looking for first good time record.
     else if(ztime.le.ztimep) then
        ierr=1
        write(6,*) &
             ' ?splitn_seek_update: ~update_time times not in ascending order:'
        write(6,*) '  current time: ',ztime,'; prior time: ',ztimep
        exit
     else
        inum_found = inum_found + 1
        if(irec_t0.eq.0) then
           ztime0=ztime
           irec_t0=irec
        endif
        ztimep = ztime
     endif
  enddo
  if(ierr.ne.0) then
     close(unit=ilun1)
     return
  endif

  if((irec_d0.gt.0).and.(irec_d0.lt.irec_t0)) then
     if((ztime0-tmin).gt.zdtmin) then
        inum_found = inum_found + 1
        tmin_leader = .TRUE.
        write(6,*) &
             ' %splitn_seek_update: data records prior to first ~update_time'
        write(6,*) &
             '  are assigned to ~update_time (tmin) = ',tmin
     else
        write(6,*) ' %splitn_seek_update: leading data records skipped:'
        write(6,*) '  The first ~update_time = ',ztime0,' is too soon after'
        write(6,*) '  the minimum time  tmin = ',tmin
     endif
  endif

  write(6,*) ' %splitn_seek_update: ',inum_found,' update blocks found in ', &
       trim(filename)

  if(inum_found.eq.0) then
     ! no updates, but, be sure to signal that a file was read...
     inum_skipped=max(1,inum_skipped)
     close(unit=ilun1)
     return
  endif

  rewind(unit=ilun1)

  call tmpfile_d('splitn_seek_update',tmp_filename,ilz)

  open(unit=ilun2,file=tmp_filename,status='unknown',iostat=ierr)
  if(ierr.ne.0) then
     write(6,*) ' ?splitn_seek_update: could not open temporary file: '
     write(6,*) '  '//tmp_filename(1:ilz)

     close(unit=ilun1)

  else
     ! OK, copy eligible records

     if(tmin_leader) then
        z20=' '
        write(z20,'(1pd19.12)') tmin
        call splitn_fput_clean(z20,'d')

        write(ilun2,'(" ~UPDATE_TIME = ",a)') trim(z20)
     endif

     irec = 0

     irec0=1
     if(inum_skipped.gt.0) then
        irec0 = irec_t0
     endif

     do
        irec = irec + 1
        read(ilun1,'(A)',iostat=ios) aline
        if(ios.ne.0) exit  ! assume EOF

        if(irec.ge.irec0) then
           write(ilun2,'(A)') trim(aline)
        endif
     enddo

     close(unit=ilun1)

     rewind(unit=ilun2)

     if(nupdate.gt.0) then
        if(tmin.le.tup(nupdate)) then
           write(6,*) ' %splitn_seek_update: old updates @ t.ge.tmin removed.'
           call splitn_update_remove(tmin,ierr)
           if(ierr.ne.0) then
              write(6,*) ' ?splitn_seek_update: unexpected splitn_update_remove error (ignored).'
           endif
        else
           kupdate=0
           nxblock=1
        endif
     else
        kupdate=0
        nxblock=1
     endif

     efitin_flag=.FALSE.

     ndblasg1 = 0
     ndblasg2 = 0
     nlackdp = 0

     do
        read(ilun2,'(A)',iostat=ios) aline
        if(ios.ne.0) exit  ! assume EOF

         iblock=nxblock

         newline= nlines+1

         call parse_line(aline,.TRUE.,ierr)
         if(ierr.ne.0) exit  ! this is an actual error

         call add_nl_line(aline)

         if(iblock.lt.nxblock) then
            ! an update block was detected -- mark line position
            linup(iblock)=nlines  ! =newline
         endif

      enddo

      if(ierr.eq.0) then
         call splitn_check_misc(ierr)
      endif

      if(ierr.eq.0) then
         close(unit=ilun2,status='delete')
      else
         write(6,*) ' ?splitn_seek_update: error detected.'
         write(6,*) '  temporary file retained here: ',trim(tmp_filename)
         close(unit=ilun2)
      endif

   endif

CONTAINS
  logical function iswhite(ichar)
    character*1 :: ichar

    ! .TRUE. if blank or TAB

    iswhite = ichar.eq.' ' .OR. ichar.eq.char(9)

  end function iswhite

end subroutine splitn_seek_update
