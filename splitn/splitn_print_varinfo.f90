subroutine splitn_print_varinfo(iluno,runid,zname,ierr)

  implicit NONE

  !  print out formatted information on a variable in a namelist

  integer, intent(in) :: iluno        ! Fortran LUN to use for output
  character*(*), intent(in) :: runid  ! namelist identification
  character*(*), intent(in) :: zname  ! variable name

  integer, intent(out) :: ierr        ! status code returned (0=OK)

  !--------------------------------------------------------
  ! local:

  character*5 :: ztype
  integer :: ichsize,irank,idims(2,10),idim1,isize,ir,iadr
  integer :: j3a,j3b,j2a,j2b,j3,j2,ii,ic,nupdate
  integer :: iupdatable,idefault,inum_update,inmax,isum,iup
  integer, dimension(:), allocatable :: iupdate_vector
  character*80 zwk
  character*32 zsave
  character*10 zpiece
  character*128 zdflt
  character*132 zchunk

  character*6 zi6
  character*8 zblk
  character*40 zcombo

  character*128, dimension(:), allocatable :: zcdata,zcdat0
  integer, dimension(:), allocatable :: idata,idat0
  logical, dimension(:), allocatable :: ildata,ildat0
  real, dimension(:), allocatable :: zrdata,zrdat0
  real*8, dimension(:), allocatable :: zddata,zddat0

  real*8, dimension(:), allocatable :: tup

  logical, dimension(:), allocatable :: dflt_flags
  logical :: all_dflt

  !--------------------------------------------------------
  !  get update count & times first...
  call splitn_update_ntimes(nupdate)
  if(nupdate.gt.0) then
     allocate(tup(nupdate))
     call splitn_update_times(nupdate,tup,ierr)
  endif

  !  OK get item data type

  call splitn_get_type(zname,ztype,ichsize,ierr)
  if(ierr.ne.0) then
     write(iluno,*) ' ?variable type could not be found:'
     write(iluno,*) trim(runid)//' namelist variable "'//trim(zname)
     return
  endif

  write(iluno,*) ' '
  write(iluno,*) trim(runid)//' namelist variable "'//trim(zname)//'":'

  if(ztype(1:1).eq.'R') then
     write(iluno,*) '  data type:  R (4 byte REAL)'//'; * denotes default value'
  else if(ztype(1:1).eq.'D') then
     write(iluno,*) '  data type:  D (8 byte REAL)'//'; * denotes default value'
  else if(ztype(1:1).eq.'I') then
     write(iluno,*) '  data type:  I (INTEGER)'//'; * denotes default value'
  else if(ztype(1:1).eq.'L') then
     write(iluno,*) '  data type:  L (LOGICAL)'//'; * denotes default value'
  else if(ztype(1:1).eq.'C') then
     write(iluno,*) '  data type:  ',ztype, &
          ' -- string length up to ',ichsize,'; * denotes default value'
  endif

  ! get and display shape of object (scalar/array)

  idim1=0
  call splitn_getdims(zname,irank,10,idims,isize,ierr)
  if(ierr.ne.0) then
     write(iluno,*) ' ?splitn_getdims error (unexpected).'
     return
  endif

  ! if rank > 3, only the rank is displayed (future upgrade possible but
  ! no rank > 3 namelist items exist at present...)

  if(irank.eq.0) then
     write(iluno,*) '  Scalar'
  else if(irank.gt.3) then
     write(iluno,*) '  Array, rank: ',irank,' (rank>3, values not displayed).'
     return
  else

     !  print the dimensions...

     idim1=idims(2,1)-idims(1,1)+1
     if(irank.eq.1) zwk='  Vector, dimension: ('
     if(irank.eq.2) zwk='  2d Array, dimension: ('
     if(irank.eq.3) zwk='  3d Array, dimension: ('
     zsave='---> '//trim(zname)//'('
     do ir=1,irank
        zpiece=' '
        write(zpiece,'(I10)') idims(1,ir)
        call addwk(zwk,zpiece,1)
        call addwk(zwk,':',1)
        if(ir.eq.1) then
           call addwk(zsave,zpiece,1)
           call addwk(zsave,':',1)
        endif
        write(zpiece,'(I10)') idims(2,ir)
        call addwk(zwk,zpiece,1)
        if(ir.eq.1) then
           call addwk(zsave,zpiece,1)
        endif
        if(ir.lt.irank) then
           call addwk(zwk,',',1)
        else
           call addwk(zwk,')',1)
        endif
     enddo
     write(iluno,*) trim(zwk),' #items = ',isize
  endif

  ! default specifier string for item (controls value when item is not
  ! explicitly in namelist)

  call splitn_get_dstr(zname,zdflt,ierr)
  if(ierr.ne.0) then
     write(iluno,*) ' ?splitn_get_dstr error (unexpected).'
     return
  endif
  write(iluno,*) '  Default value specifier: ',trim(zdflt)

  ! get update status

  inmax=max(1,nupdate)
  allocate(iupdate_vector(inmax))

  call splitn_get_udetails(zname,iupdatable, &
       idefault,inum_update,inmax,iupdate_vector,ierr)
  if(ierr.ne.0) then
     write(iluno,*) ' ?splitn_get_udetails error (unexpected).'
     return
  endif

  isum=0
  if((iupdatable.gt.1).and.(nupdate.gt.0)) then
     do ii=1,inmax
        if(iupdate_vector(ii).ge.0) isum = isum + 1
     enddo
  else
     inum_update = 0      ! just to be sure...
     iupdate_vector = -1
  endif

  if(iupdatable.eq.0) then
     write(iluno,*) '  Item not updatable.'
  else if(iupdatable.eq.1) then
     write(iluno,*) '  Item updatable via "trdat" rerun.'
  else
     write(iluno,*) '  Item updatable via TRANSP update namelist.'
  endif

  write(iluno,*) ' '
  if(iupdatable.le.1) then
     if(idefault.eq.-1) then
        write(iluno,*) '  -> Item has default value (not referenced in the namelist).'
     else if(idefault.eq.0) then
        write(iluno,*) '  -> Item is referenced in namelist but retains its default value.'
     else if(idefault.eq.1) then
        write(iluno,*) '  -> Item has non-default value(s).'
     endif
  else
     if(idefault.eq.-1) then
        write(iluno,*) '  -> Item initially has default value (not referenced in main namelist).'
     else if(idefault.eq.0) then
        write(iluno,*) '  -> Item is referenced in main namelist; initially retains its default value.'
     else if(idefault.eq.1) then
        write(iluno,*) '  -> Item starts with non-default value(s).'
     endif
  endif

  if(iupdatable.gt.1) then
     write(iluno,*) ' '
     if(isum.eq.0) then
        write(iluno,*) '  Item is not referenced in update namelists.'
     else
        write(iluno,*) '  Item is referenced ',isum,' time(s) in update namelist(s).'
        write(iluno,*) '  Value changes occur ',inum_update,' time(s).'
     endif
  endif

  ! allocate buffer and get the actual values

  allocate(dflt_flags(isize))

  if(ztype(1:1).eq.'R') then
     allocate(zrdata(isize),zrdat0(isize))
     call splitn_rgetd(zname,isize,zrdat0,ierr)
  else if(ztype(1:1).eq.'D') then
     allocate(zddata(isize),zddat0(isize))
     call splitn_dgetd(zname,isize,zddat0,ierr)
  else if(ztype(1:1).eq.'I') then
     allocate(idata(isize),idat0(isize))
     call splitn_igetd(zname,isize,idat0,ierr)
  else if(ztype(1:1).eq.'L') then
     allocate(ildata(isize),ildat0(isize))
     call splitn_lgetd(zname,isize,ildat0,ierr)
  else if(ztype(1:1).eq.'C') then
     allocate(zcdata(isize),zcdat0(isize))
     call splitn_cgetd(zname,len(zcdata(1)),isize,zcdat0,ierr)
  else
     write(iluno,*) ' ?datatype "'//ztype//'" unrecognized.'
     ierr=99
  endif
  if(ierr.ne.0) return

  ! display values

  do iup = 0, nupdate

     if((iup.gt.0).and.(iupdate_vector(iup).le.0)) cycle

     if(iup.gt.0) then
        write(iluno,*) ' '
        write(iluno,*) '  At namelist update #',iup, &
             ' @time=',tup(iup),' seconds: '
     else if(iupdatable.gt.1) then
        write(iluno,*) ' '
        write(iluno,*) '  Initially: '
     endif

     zi6=' '
     write(zi6,'(i6)') iup
     ic=1
     do ii=len(zi6)-1,1,-1
        if(zi6(ii:ii).eq.' ') then
           ic=ii+1
           exit
        endif
     enddo

     zblk = '['//zi6(ic:len(zi6))//']'
     zcombo = trim(zname)//trim(zblk)

     if(ztype(1:1).eq.'R') then
        call splitn_rget(zcombo,isize,zrdata,ierr)
     else if(ztype(1:1).eq.'D') then
        call splitn_dget(zcombo,isize,zddata,ierr)
     else if(ztype(1:1).eq.'I') then
        call splitn_iget(zcombo,isize,idata,ierr)
     else if(ztype(1:1).eq.'L') then
        call splitn_lget(zcombo,isize,ildata,ierr)
     else if(ztype(1:1).eq.'C') then
        call splitn_cget(zcombo,len(zcdata(1)),isize,zcdata,ierr)
     endif

     ! update dflt_flags(...)
     all_dflt=.TRUE.
     do ii=1,isize
        if(ztype(1:1).eq.'R') then
           dflt_flags(ii) = (zrdata(ii).eq.zrdat0(ii))
        else if(ztype(1:1).eq.'D') then
           dflt_flags(ii) = (zddata(ii).eq.zddat0(ii))
        else if(ztype(1:1).eq.'I') then
           dflt_flags(ii) = (idata(ii).eq.idat0(ii))
        else if(ztype(1:1).eq.'L') then
           dflt_flags(ii) = .NOT.( &
                   (ildata(ii).AND.(.not.ildat0(ii))) .OR. &
                   (ildat0(ii).AND.(.not.ildata(ii))) )
        else if(ztype(1:1).eq.'C') then
           dflt_flags(ii) = (zcdata(ii).eq.zcdat0(ii))
        endif
        all_dflt = all_dflt.AND.dflt_flags(ii)
     enddo

     if(ierr.ne.0) then
        write(iluno,*) ' ?splitn_print_varinfo: error reading: ',trim(zcombo)
     else
        if(irank.eq.0) then

           ! display scalar values

           if(dflt_flags(1)) then
              zdflt='*'
           else
              zdflt=' '
           endif
           if(ztype(1:1).eq.'R') then
              write(iluno,1001) trim(zname),zrdata(1),zdflt
1001          format(/' ---> ',a,' = ',1pe12.5,a1)
           else if(ztype(1:1).eq.'D') then
              write(iluno,1002) trim(zname),zddata(1),zdflt
1002          format(/' ---> ',a,' = ',1pd12.5,a1)
           else if(ztype(1:1).eq.'I') then
              write(iluno,1003) trim(zname),idata(1),zdflt
1003          format(/' ---> ',a,' = ',i10,a1)
           else if(ztype(1:1).eq.'L') then
              write(iluno,1004) trim(zname),ildata(1),zdflt
1004          format(/' ---> ',a,' = ',L1,a1)
           else if(ztype(1:1).eq.'C') then
              write(iluno,1005) trim(zname),trim(zcdata(1)),zdflt
1005          format(/' ---> ',a,' = "',a,'"',a1)
           endif
        else

           ! display vector/array values; max rank of 3

           j3a=1; j3b=1
           j2a=1; j2b=1
           if(irank.ge.3) then
              j3a=idims(1,3); j3b=idims(2,3)
           endif
           if(irank.ge.2) then
              j2a=idims(1,2); j2b=idims(2,2)
           endif
           iadr=1-idim1
           do j3=j3a,j3b
              do j2=j2a,j2b
                 iadr=iadr+idim1
                 call prinvec  ! print values (iadr:iadr+idim1-1)
              enddo
           enddo
        endif
     endif
  enddo

  ! deallocate buffers

  deallocate(dflt_flags)
  if(ztype(1:1).eq.'R') then
     deallocate(zrdata,zrdat0)
  else if(ztype(1:1).eq.'D') then
     deallocate(zddata,zddat0)
  else if(ztype(1:1).eq.'I') then
     deallocate(idata,idat0)
  else if(ztype(1:1).eq.'L') then
     deallocate(ildata,ildat0)
  else if(ztype(1:1).eq.'C') then
     deallocate(zcdata,zcdat0)
  endif

  write(iluno,*) ' '

contains

  subroutine prinvec
    integer :: iper,ict,ia,is,id,isqz

    ! print array data

    if(ztype(1:1).eq.'R') then
       iper=14
    else if(ztype(1:1).eq.'D') then
       iper=14
    else if(ztype(1:1).eq.'I') then
       iper=12
    else if(ztype(1:1).eq.'L') then
       iper=3
    else if(ztype(1:1).eq.'C') then
       iper=80
    endif

    zwk = zsave
    if(irank.eq.1) then
       call addwk(zwk,'):',1)
    else if(irank.eq.2) then
       write(zpiece,'(I10)') j2
       call addwk(zwk,','//zpiece//'):',1)
    else if(irank.eq.3) then
       write(zpiece,'(I10)') j2
       call addwk(zwk,','//zpiece,1)
       write(zpiece,'(I10)') j3
       call addwk(zwk,','//zpiece//'):',1)
    endif

    write(iluno,*) ' '
    write(iluno,'(1x,a)') trim(zwk); zwk=' '
    if(all_dflt) then
       write(iluno,*) '   (all elements are set to their default values).'
       return
    endif

    ia=iadr-1
    do ict=1,idim1
       ia=ia+1
       zchunk=' '
       isqz=1
       if(ztype(1:1).eq.'R') then
          write(zchunk,'(1pe12.5)') zrdata(ia)
       else if(ztype(1:1).eq.'D') then
          write(zchunk,'(1pd12.5)') zddata(ia)
       else if(ztype(1:1).eq.'I') then
          write(zchunk,'(i10)') idata(ia)
       else if(ztype(1:1).eq.'L') then
          write(zchunk,'(L1)') ildata(ia)
       else if(ztype(1:1).eq.'C') then
          if(zcdata(ia).eq.' ') then
             zchunk = '" "'
             isqz=0
          else
             zchunk = '"'//trim(zcdata(ia))//'"'
          endif
       endif
       if(dflt_flags(ia)) call addwk(zchunk,'*',1)
       is = len_trim(zwk)
       id = 2-iper
       do
          id=id+iper
          if(id.gt.is) then
             zwk(id:id)='.'
             call addwk(zwk,zchunk,isqz)
             zwk(id:id)=' '
             exit
          endif
       enddo
       if(len_trim(zwk)+2+iper.ge.80) then
          write(iluno,'(1x,a)') trim(zwk); zwk=' '
       endif
    enddo
    
    if(len_trim(zwk).gt.0) then
       write(iluno,'(1x,a)') trim(zwk); zwk=' '
    endif
    
  end subroutine prinvec

  subroutine addwk(zstr,zstr_add,isqueeze)
    character*(*), intent(inout) :: zstr
    character*(*), intent(in) :: zstr_add
    integer :: isqueeze

    ! append to char string, squeezing out blanks if isqueeze=1

    integer :: iii,jjj,jlen

    jlen=len_trim(zstr_add)
    iii=len_trim(zstr)
    do jjj=1,jlen
       if((isqueeze.ne.1).or.(zstr_add(jjj:jjj).ne.' ')) then
          iii=iii+1
          zstr(iii:iii)=zstr_add(jjj:jjj)
       endif
    enddo
  end subroutine addwk

end subroutine splitn_print_varinfo
