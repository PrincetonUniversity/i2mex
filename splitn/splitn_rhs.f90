subroutine splitn_rhs(zline,zrhs,jj,indx,inx,ios)

  ! parse RHS of namelist assignment to an array or range of array
  ! elements.

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zline  ! namelist line being analyzed.
  character*(*), intent(in) :: zrhs   ! RHS (no leading/trailing blanks)
  integer, intent(in) :: jj           ! index to namelist item being read
  
  integer, intent(in) :: inx          ! size of indx(...) array
  integer, intent(in) :: indx(inx)    ! indices to first array element to set

  integer, intent(out) :: ios         ! status code, 0=OK, on exit.

  !-------------------------------------------
  ! although namelist arrays can have multiple dimensions (up to 4 as
  ! of Jan. 2003), the storage associated with them in splitn_module is
  ! associated with flat singly dimensioned storage arrays-- a separate
  ! array for each data type.  So, address translation is done:
  !
  !   if a is declared a(i0:i1,j0:j1,k0:k1,m0:m1)
  !           with N1=(i1-i0+1), N2=(j1-j0+1), N3=(k1-k0+1), then,
  !
  !   a(i,j,k,m) is stored at a location offset from a base address 
  !           by this amount:
  !
  !           [[[(m-m0)*N3 + (k-k0)]*N2 + (j-j0)]*N1 + (i-i0)]*N0
  !
  !           with N0=1
  !
  ! arrays of lower dimensionality can still use this formula by 
  ! treating e.g. m0=m1=m=0
  !-------------------------------------------

  integer indims(inx)  ! Nj's
  integer idims(2,inx) ! actual array dimensioning information

  integer iadr_base    ! base address for array object
  integer iadr0        ! address of first element to be set from RHS
  integer iadrmax      ! maximum allowable address settable from RHS
                       ! (determined by size of namelist array object)
  integer iadrvmx      ! maximum address staying within the same array
                       ! row, i.e. without changing any but the first
                       ! of the array element indices.

  integer iadrnlin     ! address (in ilines array) for storing current
                       ! line number (associating data with file line no.)

  integer i,ix,iscan,ilen
  integer iw1,iw2,idelim,ir1,ir2,irepeat,iextra
  integer icura,imaxa,icurl,imaxl,iwarn,iwarn1,iwarn2,ier
  integer :: itmp
  integer :: iadls

  real, dimension(:), allocatable :: rsave
  real*8, dimension(:), allocatable :: dsave
  integer, dimension(:), allocatable :: isave
  logical, dimension(:), allocatable :: lsave
  character*150, dimension(:), allocatable :: chsave
  integer, dimension(:), allocatable :: ichange
  integer ii,iparen,inam2
  !-------------------------------------------

  iwarn=0
  iwarn1=0
  iwarn2=0
  ios=0

  ilen=len(zrhs)
  if(inx.ne.maxrank) then
     ios=-99
     write(6,*) ' ?splitn: internal error, maximum array rank inconsistency.'
     return
  endif

  idims = varlist(jj)%dims

  itmp = varlist(jj)%addr_st
  if(kupdate.eq.0) then
     iadr_base = varlist(jj)%addr
  else
     iadr_base = itmp
     if(varlist(jj)%type.eq.'R') then
        iadr_base = nreal + (kupdate-1)*nreal_st + iadr_base
     else if(varlist(jj)%type.eq.'I') then
        iadr_base = nint + (kupdate-1)*nint_st + iadr_base
     else if(varlist(jj)%type.eq.'L') then
        iadr_base = nlog + (kupdate-1)*nlog_st + iadr_base
     else if(varlist(jj)%type.eq.'D') then
        iadr_base = nr8 + (kupdate-1)*nr8_st + iadr_base
     else if(varlist(jj)%type.eq.'C') then
        iadr_base = nchv + (kupdate-1)*nchv_st + iadr_base
     endif
  endif

  indims(1)=1
  do i=1,inx-1
     indims(i+1)=(idims(2,i)-idims(1,i))+1
  enddo

  ! calculate addresses as offsets first

  iadr0=0
  iadrvmx=0
  iadrmax=0

  do i=inx,1,-1
     iadr0=(iadr0+(indx(i)-idims(1,i)))*indims(i)
     iadrmax=(iadrmax+(idims(2,i)-idims(1,i)))*indims(i)
     if(i.eq.1) then
        ix=idims(2,i)
     else
        ix=indx(i)
     endif
     iadrvmx=(iadrvmx+(ix-idims(1,i)))*indims(i)
  enddo

  ! add in base to get actual addresses

  iadrnlin=iadr0+varlist(jj)%nlinadr   ! file line no. storage address
  if(kupdate.gt.0) then
     iadrnlin = iadr0+varlist(jj)%nlinadr_st
     iadrnlin = nreal+nint+nlog+nr8+nchv + &
          (kupdate-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st) + iadrnlin
  endif

  itmp=iadr0+itmp
  iadr0=iadr0+iadr_base
  iadrmax=iadrmax+iadr_base
  iadrvmx=iadrvmx+iadr_base

  allocate(ichange(iadr0:iadrmax))

  !---------------
  ! OK now scan the RHS

  ichange=0   ! note any real changes in value...

  icura=iadr0

  iscan=0
  do

     call splitn_nxgrp(zrhs,iscan,iw1,iw2,idelim,ir1,ir2,irepeat,iextra)
     if(iextra.gt.0) then
        write(6,*) ' ?splitn_rhs: extraneous characters after quote string.'
        return
     endif

     imaxa=icura+irepeat-1
     if(imaxa.gt.iadrmax) then
        write(6,*) ' ?splitn_rhs: namelist assignment exceeds array size.'
        write(6,*) zline
        ios=1
        return
     else if(imaxa.gt.iadrvmx) then
        if(iwarn.eq.0) then
           iwarn=iwarn+1
           write(6,*) ' ------------- '
           write(6,*) &
                ' %splitn_rhs: warning: namelist assign spans multiple', &
                ' rows of array.'
           write(6,*) zline
        endif
     endif

     ! addresses appear to be OK

     if(iw2.ge.iw1) then
        ! non null value; decode according to type
        if(varlist(jj)%type.eq.'R') then
           allocate(rsave(icura:imaxa))
           rsave=rbuf(icura:imaxa)
           call splitn_rdcod(zrhs(iw1:iw2),rbuf(icura:imaxa),irepeat,ier)
           do ii=icura,imaxa
              if(rbuf(ii).ne.rsave(ii)) ichange(ii)=1
           enddo
           deallocate(rsave)
           if(varlist(jj)%steerable.eq.2) then
              iadls = varlist(jj)%nlinadr_st
              call rset_ucop(rbuf,icura,irepeat,itmp,iadls,kupdate,nupdate_max)
           endif
        else if(varlist(jj)%type.eq.'I') then
           allocate(isave(icura:imaxa))
           isave=intbuf(icura:imaxa)
           call splitn_idcod(zrhs(iw1:iw2),intbuf(icura:imaxa),irepeat,ier)
           do ii=icura,imaxa
              if(intbuf(ii).ne.isave(ii)) ichange(ii)=1
           enddo
           deallocate(isave)
           if(varlist(jj)%steerable.eq.2) then
              iadls = varlist(jj)%nlinadr_st
              call iset_ucop(intbuf,icura,irepeat,itmp,iadls,kupdate,nupdate_max)
           endif
        else if(varlist(jj)%type.eq.'L') then
           allocate(lsave(icura:imaxa))
           lsave=logbuf(icura:imaxa)
           call splitn_ldcod(zrhs(iw1:iw2),logbuf(icura:imaxa),irepeat,ier)
           do ii=icura,imaxa
              if(lsave(ii)) then
                 if(.not.logbuf(ii)) ichange(ii)=1
              else
                 if(logbuf(ii)) ichange(ii)=1
              endif
           enddo
           deallocate(lsave)
           if(varlist(jj)%steerable.eq.2) then
              iadls = varlist(jj)%nlinadr_st
              call lset_ucop(logbuf,icura,irepeat,itmp,iadls,kupdate,nupdate_max)
           endif
        else if(varlist(jj)%type.eq.'D') then
           allocate(dsave(icura:imaxa))
           dsave=dbuf(icura:imaxa)
           call splitn_ddcod(zrhs(iw1:iw2),dbuf(icura:imaxa),irepeat,ier)
           do ii=icura,imaxa
              if(dbuf(ii).ne.dsave(ii)) ichange(ii)=1
           enddo
           deallocate(dsave)
           if(varlist(jj)%steerable.eq.2) then
              iadls = varlist(jj)%nlinadr_st
              call dset_ucop(dbuf,icura,irepeat,itmp,iadls,kupdate,nupdate_max)
           endif
        else if(varlist(jj)%type(1:1).eq.'C') then
           allocate(chsave(icura:imaxa))
           chsave=chbuf(icura:imaxa)
           call splitn_cdcod(zrhs(iw1:iw2), &
                varlist(jj)%chsize,chbuf(icura:imaxa),irepeat,ier)
           do ii=icura,imaxa
              if(chbuf(ii).ne.chsave(ii)) ichange(ii)=1
           enddo
           deallocate(chsave)
           if(varlist(jj)%steerable.eq.2) then
              iadls = varlist(jj)%nlinadr_st
              call chset_ucop(chbuf,icura,irepeat,itmp,iadls,kupdate,nupdate_max)
           endif
        endif
        if(ier.ne.0) then
           write(6,*) ' ?splitn_rhs: error parsing namelist array values:'
           write(6,*) zline
           ios=1
           return
        endif
     else
        write(6,*) ' ?splitn_rhs: null value field.'
        write(6,*) '  line starts: ',trim(zline)
        ier=1
        ios=1
        return
     endif
     iscan = idelim

     ! data values have been assigned.  Now put in links back to the file
     ! line-- set sign bit in case of null value assignment.

     ios=0
     icurl=iadrnlin+(icura-iadr0)
     imaxl=iadrnlin+(imaxa-iadr0)
     do i=icurl,imaxl
        ii=i-iadrnlin+iadr0
        if(ilines(i).ne.0) then
           if(ichange(ii).eq.1) then  ! detect multiply assigned elements
              iwarn2=iwarn2+1         ! WITH change in value
           else
              iwarn1=iwarn1+1         ! NO change in value
           endif
        endif

        iline_nstack = iline_nstack + 1
        iline_stack(1,iline_nstack) = i
        iline_stack(2,iline_nstack) = newline

        !xxxx   ilines(i)=newline  -- deferred, see add_nl_line in module

        ivrange(:,i)=0
        irrange(:,i)=0
        krepeat(i)=0

        if(iw2.ge.iw1) then           ! last assignment wins
           ivrange(1,i)=iw1+val1-1
           ivrange(2,i)=iw2+val1-1
        endif
        if(ir1.gt.0) then
           krepeat(i)=irepeat
           irrange(1,i)=ir1+val1-1
           irrange(2,i)=ir2+val1-1
        endif
     enddo

     if(idelim.ge.ilen) exit
     icura=icura+irepeat
     itmp =itmp +irepeat

  enddo

  if(max(iwarn1,iwarn2).gt.0) then
     iparen=index(zline,'(')
     if(iparen.le.nam1) iparen=nam2+1
     if(iparen.lt.nam2) then
        inam2=iparen-1
     else
        inam2=nam2
     endif
  endif

  if(iwarn1.gt.0) call splitn_addwarn(zline(nam1:inam2),dblasg1,ndblasg1, &
       len(dblasg1(1)),kdblasg1,maxwarn)
  if(iwarn2.gt.0) call splitn_addwarn(zline(nam1:inam2),dblasg2,ndblasg2, &
       len(dblasg2(1)),kdblasg2,maxwarn)

  deallocate(ichange)

end subroutine splitn_rhs
