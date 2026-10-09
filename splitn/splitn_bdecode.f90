subroutine splitn_bdecode1(subname,zinput,zname,iblock,ierr)

  !  decode string <name> -> <name>(uppercase) & iblock=0
  !                <zname>[<num>] -> <name>(uppercase) & iblock=<num>

  !  example: "alpha_e[2]" -> zname="ALPHA_E", iblock=2

  !  error condition: <num> fails to decode as integer, (i8)
  !---------------------------

  character*(*), intent(in) :: subname    ! name of calling subroutine
  character*(*), intent(in) :: zinput     ! string containing: <zname>[<num>]

  character*(*), intent(out) :: zname     ! name only, output in uppercase
  integer, intent(out) :: iblock          ! block number (0 if absent)
  integer, intent(out) :: ierr            ! status code (0=OK)

  !---------------------------
  ! local:
  integer :: ibr1,ibr2   ! brackets
  integer :: ilen
  character*8 zi8
  !---------------------------

  ierr = 0
  iblock = 0
  zname = '??'

  ibr1 = index(zinput,'[')
  ibr2 = index(zinput,']')

  if((ibr1.eq.ibr2).and.(ibr1.le.0)) then

     ! no bracketed number

     iblock = 0
     zname = zinput
     call uupper(zname)

  else

     if(ibr1.eq.1) then
        write(6,*) '?'//trim(subname)//': no name before "[" in input:'
        write(6,*) ' '//trim(zinput)
        ierr = 1
        return
     endif

     zname = zinput(1:ibr1-1)
     call uupper(zname)

     if(ibr2.lt.ibr1) then
        write(6,*) '?'//trim(subname)//': "]" before "[" in input:'
        write(6,*) ' '//trim(zinput)
        ierr = 1
        return
     endif

     if(ibr2.gt.(ibr1+1)) then
        ilen = ibr2 - ibr1 - 1
        if(ilen.gt.8) then
           write(6,*) '?'//trim(subname)//': argument inside "[,]" too long:'
           write(6,*) ' '//trim(zinput)
           ierr = 1
           return
        endif

        zi8=' '
        zi8(8-ilen+1:8) = zinput(ibr1+1:ibr2-1)
        read(zi8,'(I8)',iostat=ierr) iblock
        if(ierr.ne.0) then
           write(6,*) '?'//trim(subname)//': integer decode failed inside "[,]".'
           write(6,*) ' '//trim(zinput)
           ierr = 1
           iblock = 0
           return
        endif

     endif

  endif

end subroutine splitn_bdecode1

subroutine splitn_bdecode2(subname,zinput,zname,iblock1,iblock2,ierr)

  !  decode string <name> -> <name>(uppercase) & iblock1=iblock2=0
  !                <zname>[<num>] -> <name>(uppercase) & 
  !                                     iblock1=iblock2=<num>
  !                <zname>[<num1>:<num2>] -> <name>(uppercase) & 
  !                                     iblock1=<num1>, iblock2=<num2>
  !                <zname>[*] -> <name>(uppercase) & iblock1=0 & iblock2=N
  !                                     where N=nupdate = #updates that exist.

  !  examples: "alpha_e[2]" -> zname="ALPHA_E", iblock1=iblock2=2
  !            "alpha_e"    -> zname="ALPHA_E", iblock1=iblock2=0
  !            "alpha_e[1:3]" -> zname="ALPHA_E", iblock1=1, iblock2=3
  !            "alpha_e[*]" -> zname="ALPHA_E", iblock1=0, iblock2=nupdate

  !  error condition: <num> fails to decode as integer, (i8)
  !---------------------------

  character*(*), intent(in) :: subname    ! name of calling subroutine
  character*(*), intent(in) :: zinput     ! string containing: <zname>[<num>]

  character*(*), intent(out) :: zname     ! name only, output in uppercase
  integer, intent(out) :: iblock1,iblock2 ! range of block numbers out
  integer, intent(out) :: ierr            ! status code (0=OK)

  !---------------------------
  ! local:
  integer :: ibr1,ibr2   ! brackets
  integer :: icolon      ! colon
  integer :: ilen
  integer :: nupdate
  character*8 zi8
  !---------------------------

  ierr = 0
  iblock1 = 0
  iblock2 = 0
  zname = '??'

  call splitn_update_ntimes(nupdate)

  ibr1 = index(zinput,'[')
  ibr2 = index(zinput,']')
  icolon = index(zinput,':')

  if((ibr1.eq.ibr2).and.(ibr1.le.0)) then

     ! no bracketed number(s)

     zname = zinput
     call uupper(zname)

  else

     if(ibr1.eq.1) then
        call errmsg('no name before "[" in input:')
        return
     endif

     zname = zinput(1:ibr1-1)
     call uupper(zname)

     if(ibr2.lt.ibr1) then
        call errmsg('"]" before "[" in input:')
        return
     endif

     ! ibr2 > ibr1 now established.

     if(zinput(ibr1:ibr2).eq.'[*]') then
        iblock1=0
        iblock2=nupdate
        ierr = 0
        return
     endif

     if(icolon.gt.0) then
        if((icolon.lt.ibr1).or.(icolon.gt.ibr2)) then
           call errmsg('misplaced colon.')
           return
        endif
     endif

     ! if icolon gt 0:  ibr1 < icolon < ibr2 established.

     if(zinput(ibr1:ibr2).eq.'[:]') then
        iblock1=0
        iblock2=nupdate
        ierr = 0
        return
     endif

     if(icolon.eq.0) then
        if(ibr2.gt.(ibr1+1)) then
           ilen = ibr2 - ibr1 - 1
           call chklen; if(ierr.ne.0) return

           zi8=' '
           zi8(8-ilen+1:8) = zinput(ibr1+1:ibr2-1)
           read(zi8,'(I8)',iostat=ierr) iblock1
           if(ierr.ne.0) then
              call errmsg('integer decode failed inside "[,]".')
              iblock1 = 0
              return
           endif
           iblock2 = iblock1
        else
           ! "[]"
           iblock1 = 0
           iblock2 = 0
        endif

     else
        if(icolon.eq.ibr1+1) then
           iblock1=0
        else
           ilen = icolon - ibr1 - 1
           call chklen; if(ierr.ne.0) return

           zi8=' '
           zi8(8-ilen+1:8) = zinput(ibr1+1:icolon-1)
           read(zi8,'(I8)',iostat=ierr) iblock1
           if(ierr.ne.0) then
              call errmsg('integer decode failed inside "[,]".')
              iblock1 = 0
              return
           else
              iblock1=max(0,min(nupdate,iblock1))
           endif
        endif

        if(icolon.eq.ibr2-1) then
           iblock2=nupdate
        else
           ilen = ibr2 - icolon - 1
           call chklen; if(ierr.ne.0) return

           zi8=' '
           zi8(8-ilen+1:8) = zinput(icolon+1:ibr2-1)
           read(zi8,'(I8)',iostat=ierr) iblock2
           if(ierr.ne.0) then
              call errmsg('integer decode failed inside "[,]".')
              iblock2 = iblock1
              return
           else
              iblock2=max(0,min(nupdate,iblock2))
           endif

           if(iblock1.gt.iblock2) then
              call errmsg('integer indices out of order.')
           endif

        endif


     endif

  endif

CONTAINS

  subroutine chklen
    if(ilen.gt.8) then
       write(6,*) '?'//trim(subname)//': argument inside "[,]" too long:'
       write(6,*) ' '//trim(zinput)
       ierr = 1
    endif
  end subroutine chklen

  subroutine errmsg(msg)
    character*(*), intent(in) :: msg

    write(6,*) '?'//trim(subname)//': '//trim(msg)
    write(6,*) ' '//trim(zinput)
    ierr = 1

  end subroutine errmsg

end subroutine splitn_bdecode2
