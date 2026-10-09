subroutine trx_wr_protimes(ppath,ierr)

  ! write file of profile output times -- 1st record = #of times
  ! subsequently, one time per record

  ! use data recorded in module -- report error if none is found

  use trx_module
  use trxplib_ps_options
  implicit NONE

  character*(*), intent(in) :: ppath
  integer, intent(out) :: ierr

  !------------------------------------------
  character*30 :: prev_file
  integer :: it,inum,ilun
  integer :: lunzer
  !------------------------------------------

  ierr=0
  if(.not.allocated(time_pr)) then
     write(lunzer(0),*) &
          ' ?trx_wr_protimes: no profile time data; connect to run first.'
     ierr=1
     return
  endif

  inum=nprtime

  call find_io_unit(ilun)

  open(unit=ilun,file=trim(ppath),status='unknown',iostat=ierr)
  if(ierr.ne.0) then
     write(lunzer(0),*) &
          ' ?trx_wr_protimes: open failure, filename was: '//trim(ppath)
     return
  endif

  write(ilun, &
       '(1x,i5,"  ! #profile (1x,i5) // times (1x,1pe13.6) one per line")') &
       inum

  prev_file = ' '
  do it=1,nprtime
     call trx_wr_encode_r4(ilun,time_pr(it),ps_prefix,prev_file)
  enddo

  close(unit=ilun)

end subroutine trx_wr_protimes

subroutine trx_wr_encode_r4(ilun,ztime,zprefix,prev_file)

  ! write protimes.dat or sawtimes.dat lines;
  ! append filename substring

  ! in output record, character positions(1:20) : 1x,1pe13.6
  ! in output record, character positions(20:)  : filename piece

  ! time argument in R4

  implicit NONE

  integer, intent(in) :: ilun
  real, intent(in) :: ztime
  character*(*), intent(in) :: zprefix
  character*(*), intent(inout) :: prev_file

  !----------------------------------
  ! LOCAL:
  character*20 time_encod
  character*30 file_encod,wk
  character*10 zint,zmant
  character*1 char1
  integer :: ic,id,ii,ilen,iexp,idot,inb
  integer :: jexp
  !----------------------------------

  if(zprefix.eq.' ') then
     write(ilun,'(1x,1pe13.6)') ztime
     return
  endif

  time_encod=' '
  write(time_encod,'(1x,1pe13.6)') ztime
  file_encod = time_encod

  ilen = len(trim(file_encod))

  ! find exponent character

  do ic=ilen,1,-1
     char1 = file_encod(ic:ic)
     if((char1.eq.'e').OR.(char1.eq.'E').OR. &
          (char1.eq.'d').OR.(char1.eq.'D')) then
        iexp=ic
        exit
     endif
  enddo
  zint = file_encod(iexp+1:)
  zint = adjustR(zint)
  read(zint,'(I10)') jexp

  if((jexp.gt.6).or.(jexp.lt.-2)) then
     ! preserve scientific notation in filename piece
     wk = file_encod
     do ic=1,ilen
        char1 = wk(ic:ic)
        if(char1.eq.'-') then
           wk(ic:ic)='m'
        else if(char1.eq.'+') then
           wk(ic:ic)='p'
        else if(char1.eq.'.') then
           wk(ic:ic)='x'
        else if((char1.eq.'e').OR.(char1.eq.'E').OR. &
          (char1.eq.'d').OR.(char1.eq.'D')) then
           wk(ic:ic)='E'
        endif
     enddo
  else
     ! discard exponent; move decimal point
     file_encod(iexp:)=' '
     ilen = len(trim(file_encod))

     idot=0
     inb=0
     do ic=1,ilen
        if(file_encod(ic:ic).ne.' ') then
           if(inb.eq.0) inb=ic
        endif
        if(file_encod(ic:ic).eq.'.') then
           idot=ic
           exit
        endif
     enddo

     wk = ' '

     if(jexp.lt.0) then
        if(ztime.lt.(0.0)) then
           if(jexp.eq.-1) then
              wk(inb:) = 'm0x'
              id=0
              ii=inb+2
           else if(jexp.eq.-2) then
              wk(inb:) = 'm0x0'
              id=1
              ii=inb+3
           endif
        else
           if(jexp.eq.-1) then
              wk(inb:) = '0x'
              id=0
              ii=inb+1
           else if(jexp.eq.-2) then
              wk(inb:) = '0x0'
              id=1
              ii=inb+2
           endif
        endif

        do ic=1,ilen
           if(file_encod(ic:ic).eq.'+') cycle
           if(file_encod(ic:ic).eq.'-') cycle
           if(file_encod(ic:ic).eq.'.') cycle
           if(file_encod(ic:ic).eq.' ') cycle
           ii = ii + 1
           id = id + 1
           wk(ii:ii)=file_encod(ic:ic)
           if(id.ge.6) exit
        enddo

     else
        if(ztime.lt.(0.0)) then
           wk(inb:) = 'm'
           id=-1
           ii=inb
        else
           id=-1
           ii=inb-1
        endif

        do ic=1,ilen
           if(file_encod(ic:ic).eq.'+') cycle
           if(file_encod(ic:ic).eq.'-') cycle
           if(file_encod(ic:ic).eq.'.') cycle
           if(file_encod(ic:ic).eq.' ') cycle
           ii = ii + 1
           id = id + 1
           wk(ii:ii)=file_encod(ic:ic)
           if(id.ge.6) exit
           jexp = jexp-1
           if(jexp.eq.-1) then
              ii = ii + 1
              wk(ii:ii)='x'
           endif
        enddo
     endif

  endif

  wk = adjustL(wk)
  file_encod = trim(zprefix)//trim(wk)//'.cdf'

  if(file_encod.eq.prev_file) then
     file_encod = trim(file_encod)//'a'
  endif

  prev_file = file_encod

  write(ilun,'(a,a)') time_encod,trim(file_encod)

end subroutine trx_wr_encode_r4
