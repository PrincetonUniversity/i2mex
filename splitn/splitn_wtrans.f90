subroutine splitn_wtrans_test

  ! test driver for splitn_wtrans:  load some random changes into
  ! the write buffer arrays; call splitn_wtrans to transfer these
  ! changes into the main buffers and express these changes in the
  ! file text form.

  ! execute after reading a namelist (when dbuf is copied into dbuf_w
  ! and similarly for all the other buffers)

  use splitn_module
  implicit NONE

  !-------------------------------
  integer i,i2,istart,iinc
  !-------------------------------

  kupdate = 0   ! applies to base namelist only

  iinc=1000

  if(nreal.gt.0) then

     istart=max(1,min(nreal,iinc)/2)
     do i=istart,nreal,iinc
        i2=min(nreal,i+2)
        rbuf_w(i:i2)=(rbuf(i)-0.5e0)*2
     enddo

  endif

  if(nr8.gt.0) then

     istart=max(1,min(nr8,iinc)/2)
     do i=istart,nr8,iinc
        i2=min(nr8,i+2)
        dbuf_w(i:i2)=(dbuf(i)-0.5d0)*2
     enddo

  endif

  if(nint.gt.0) then

     istart=max(1,min(nint,iinc)/2)
     do i=istart,nint,iinc
        i2=min(nint,i+2)
        intbuf_w(i:i2)=(intbuf(i)-1)*2
     enddo

  endif

  if(nlog.gt.0) then

     istart=max(1,min(nlog,iinc)/2)
     do i=istart,nlog,iinc
        i2=min(nlog,i+2)
        logbuf_w(i:i2)=.TRUE.
     enddo

  endif

  if(nchv.gt.0) then

     istart=max(1,min(nchv,iinc)/2)
     do i=istart,nchv,iinc
        i2=min(nchv,i+2)
        chbuf_w(i:i2)='Hi'
     enddo

  endif

  call splitn_wtrans

end subroutine splitn_wtrans_test

subroutine splitn_wtrans

  !  transcribe write buffer into main buffer, with update to file text

  use splitn_module
  implicit NONE

  !-----------------------------------
  integer i,idim,iadr0,iadr1,iadc0,iadc1,isize,indx(10),inum,ierr
  integer irank

  logical diff
  !-----------------------------------

  !  loop over all namelist variables

  do i=1,nvars
     isize=1
     irank=varlist(i)%rank
     do idim=1,irank
        isize=isize*(varlist(i)%dims(2,idim)-varlist(i)%dims(1,idim)+1)
     enddo
     iadr0=varlist(i)%addr
     iadr1=iadr0+isize-1
     iadc0=iadr0-1

     do
        !  look for start of change block
        iadc0=iadc0+1
        if(iadc0.gt.iadr1) exit
        diff = .not.is_equal(iadc0)
        if(diff) then
           !  look for end of change block with same value
           iadc1=iadc0-1
           do
              iadc1=iadc1+1
              if(iadc1.eq.iadr1) exit
              if(.not.is_equal2(iadc0,iadc1+1)) exit
           enddo

           inum=iadc1-iadc0+1
           
           if(irank.gt.0) then
              call gen_indx
           endif
           call put

           iadc0=iadc1  ! and resume search...
        endif
     enddo
  enddo
  contains
    subroutine gen_indx
      integer kblock,ind,iwk,irat,ii

      !  generate index 
      indx=0
      iwk=iadc0-iadr0
      do ii=irank,1,-1
         kblock=(varlist(i)%dims(2,ii)-varlist(i)%dims(1,ii))+1
         irat=iwk/kblock
         ind = (iwk - irat*kblock) + varlist(i)%dims(1,ii)
         indx(ii)=ind
         iwk=irat
      enddo
    end subroutine gen_indx

    subroutine put

      !  request transfer

      if(varlist(i)%type.eq.'R') then
         if(irank.eq.0) then
            call splitn_rput_sc(varlist(i)%name,rbuf_w(iadc0),ierr)
         else
            call splitn_rput_ar(varlist(i)%name,rbuf_w(iadc0), &
                 irank,indx(1:irank),inum,ierr)
         endif

      else if(varlist(i)%type.eq.'D') then
         if(irank.eq.0) then
            call splitn_dput_sc(varlist(i)%name,dbuf_w(iadc0),ierr)
         else
            call splitn_dput_ar(varlist(i)%name,dbuf_w(iadc0), &
                 irank,indx(1:irank),inum,ierr)
         endif

      else if(varlist(i)%type.eq.'I') then
         if(irank.eq.0) then
            call splitn_iput_sc(varlist(i)%name,intbuf_w(iadc0),ierr)
         else
            call splitn_iput_ar(varlist(i)%name,intbuf_w(iadc0), &
                 irank,indx(1:irank),inum,ierr)
         endif

      else if(varlist(i)%type.eq.'L') then
         if(irank.eq.0) then
            call splitn_lput_sc(varlist(i)%name,logbuf_w(iadc0),ierr)
         else
            call splitn_lput_ar(varlist(i)%name,logbuf_w(iadc0), &
                 irank,indx(1:irank),inum,ierr)
         endif

      else if(varlist(i)%type(1:1).eq.'C') then
         if(irank.eq.0) then
            call splitn_chput_sc(varlist(i)%name,chbuf_w(iadc0),ierr)
         else
            call splitn_chput_ar(varlist(i)%name,chbuf_w(iadc0), &
                 irank,indx(1:irank),inum,ierr)
         endif

      endif

    end subroutine put

    logical function is_equal(iadc)

      !  return TRUE if main buffer and write buffer elements are equal

      integer, intent(in) :: iadc   ! buffer address

      if(varlist(i)%type.eq.'R') then
         is_equal = rbuf(iadc).eq.rbuf_w(iadc)
      else if(varlist(i)%type.eq.'D') then
         is_equal = dbuf(iadc).eq.dbuf_w(iadc)
      else if(varlist(i)%type.eq.'I') then
         is_equal = intbuf(iadc).eq.intbuf_w(iadc)
      else if(varlist(i)%type.eq.'L') then
         is_equal = logbuf(iadc).and.logbuf_w(iadc)
         is_equal = is_equal.or.((.not.logbuf(iadc)).and.(.not.logbuf_w(iadc)))
      else if(varlist(i)%type(1:1).eq.'C') then
         is_equal = chbuf(iadc).eq.chbuf_w(iadc)
      endif

    end function is_equal
           
    logical function is_equal2(iad1,iad2)

      !  return TRUE if two write buffer elements are equal

      integer, intent(in) :: iad1,iad2   ! buffer addresses

      if(varlist(i)%type.eq.'R') then
         is_equal2 = rbuf_w(iad1).eq.rbuf_w(iad2)
      else if(varlist(i)%type.eq.'D') then
         is_equal2 = dbuf_w(iad1).eq.dbuf_w(iad2)
      else if(varlist(i)%type.eq.'I') then
         is_equal2 = intbuf_w(iad1).eq.intbuf_w(iad2)
      else if(varlist(i)%type.eq.'L') then
         is_equal2 = logbuf_w(iad1).and.logbuf_w(iad2)
         is_equal2 = is_equal2.or. &
              ((.not.logbuf_w(iad1)).and.(.not.logbuf_w(iad2)))
      else if(varlist(i)%type(1:1).eq.'C') then
         is_equal2 = chbuf_w(iad1).eq.chbuf_w(iad2)
      endif

    end function is_equal2
           
end subroutine splitn_wtrans
