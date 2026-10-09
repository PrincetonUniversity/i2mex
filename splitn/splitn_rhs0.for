      subroutine splitn_rhs0(zline,zrhs,ivals,jj)
C
      use splitn_module
      implicit NONE
C
C  check RHS of scalar value assignment; verify there is only one value.
C
C  input:
C
      character*(*), intent(in) :: zline  ! input line
      character*(*), intent(in) :: zrhs ! RHS of namelist assignment
C                    **NO leading or trailing blanks**
C
      integer,intent(out) :: ivals      ! #values found (1 expected).
      integer,intent(in) :: jj          ! index to namelist item definition.
C
C  output:
C
CC      integer ivals                     ! =1 if only one RHS value
C                   =0 if type error
C                   .gt.1 if multiple RHS values.
C
C---------------------------------------------
      integer iextra,irepeat,idelim,iw1,iw2,iadr,iadlp,itmp,ier
      integer ir1,ir2
      integer :: iadls
C
      real rsave(1)
      real*8 dsave(1)
      integer isave(1)
      logical lsave(1)
      character*150 chsave(1)
C
      integer ichange                   ! detect change in value
C
C---------------------------------------------
C
      ivals=1
      call splitn_nxgrp(zrhs,0,iw1,iw2,idelim,ir1,ir2,irepeat,iextra)
      if(iextra.gt.0) then
         write(6,*)
     >      ' ?splitn: extraneous characters after quote string.'
         ivals=0
         return
      endif
C
      if(irepeat.ne.1) then
         write(6,*) ' ?splitn_rhs0: repeat count for scalar value.'
         ivals=irepeat
         return
      else
C
         ichange=0
         ivals=1
         if(kupdate.eq.0) then
            iadr=varlist(jj)%addr
            iadlp=varlist(jj)%nlinadr
            itmp=varlist(jj)%addr_st
         else
            ! we got through splitn_lhs; we know this is an updatable quantity
            itmp=varlist(jj)%nlinadr_st
            iadlp = nreal+nint+nlog+nr8+nchv +
     >           (kupdate-1)*(nreal_st+nint_st+nlog_st+nr8_st+nchv_st)
     >           + itmp
            itmp=varlist(jj)%addr_st
         endif
         if(iw2.ge.iw1) then
            if(varlist(jj)%type.eq.'R') then
               if(kupdate.gt.0) then
                  iadr = nreal + (kupdate-1)*nreal_st + itmp
               endif
               rsave=rbuf(iadr:iadr)
               call splitn_rdcod(zrhs(iw1:iw2),rbuf(iadr:iadr),1,ier)
               if(rbuf(iadr).ne.rsave(1)) ichange=1
               if(varlist(jj)%steerable.eq.2) then
                  iadls = varlist(jj)%nlinadr_st
                  call rset_ucop(rbuf,iadr,1,itmp,iadls,kupdate,
     >                 nupdate_max)
               endif
            else if(varlist(jj)%type.eq.'I') then
               if(kupdate.gt.0) then
                  iadr = nint + (kupdate-1)*nint_st + itmp
               endif
               isave=intbuf(iadr:iadr)
               call splitn_idcod(zrhs(iw1:iw2),intbuf(iadr:iadr),1,ier)
               if(intbuf(iadr).ne.isave(1)) ichange=1
               if(varlist(jj)%steerable.eq.2) then
                  iadls = varlist(jj)%nlinadr_st
                  call iset_ucop(intbuf,iadr,1,itmp,iadls,kupdate,
     >                 nupdate_max)
               endif
            else if(varlist(jj)%type.eq.'L') then
               if(kupdate.gt.0) then
                  iadr = nlog + (kupdate-1)*nlog_st + itmp
               endif
               lsave=logbuf(iadr:iadr)
               call splitn_ldcod(zrhs(iw1:iw2),logbuf(iadr:iadr),1,ier)
               if(lsave(1)) then
                  if(.not.logbuf(iadr)) ichange=1
               else
                  if(logbuf(iadr)) ichange=1
               endif
               if(varlist(jj)%steerable.eq.2) then
                  iadls = varlist(jj)%nlinadr_st
                  call lset_ucop(logbuf,iadr,1,itmp,iadls,kupdate,
     >                 nupdate_max)
               endif
            else if(varlist(jj)%type.eq.'D') then
               if(kupdate.gt.0) then
                  iadr = nr8 + (kupdate-1)*nr8_st + itmp
               endif
               dsave=dbuf(iadr:iadr)
               call splitn_ddcod(zrhs(iw1:iw2),dbuf(iadr:iadr),1,ier)
               if(dbuf(iadr).ne.dsave(1)) ichange=1
               if(varlist(jj)%steerable.eq.2) then
                  iadls = varlist(jj)%nlinadr_st
                  call dset_ucop(dbuf,iadr,1,itmp,iadls,kupdate,
     >                 nupdate_max)
               endif
            else if(varlist(jj)%type(1:1).eq.'C') then
               if(kupdate.gt.0) then
                  iadr = nchv + (kupdate-1)*nchv_st + itmp
               endif
               chsave=chbuf(iadr:iadr)
               call splitn_cdcod(zrhs(iw1:iw2),
     >            varlist(jj)%chsize,chbuf(iadr:iadr),1,ier)
               if(chbuf(iadr).ne.chsave(1)) ichange=1
               if(varlist(jj)%steerable.eq.2) then
                  iadls = varlist(jj)%nlinadr_st
                  call chset_ucop(chbuf,iadr,1,itmp,iadls,kupdate,
     >                 nupdate_max)
               endif
            endif
         else
            write(6,*) ' ?splitn_rhs0: null value field.'
            write(6,*) '  line starts: ',trim(zline)
            ier=1                       ! null value
         endif

         if(ier.gt.0) then
            ivals=0                     ! type mismatch
         else
            if(idelim.lt.len(zrhs)) then
               write(6,*) ' ?splitn_rhs0: multiple values for scalar.'
               ivals=2                  ! multiple values
            endif
         endif
      endif
C
      if(ivals.eq.1) then
         if(ilines(iadlp).ne.0) then
            if(ichange.eq.0) then
               call splitn_addwarn(zline(nam1:nam2),dblasg1,ndblasg1,
     >            len(dblasg1(1)),kdblasg1,maxwarn)
            else
               call splitn_addwarn(zline(nam1:nam2),dblasg2,ndblasg2,
     >            len(dblasg2(1)),kdblasg2,maxwarn)
            endif
         endif

         iline_nstack = iline_nstack + 1
         iline_stack(1,iline_nstack) = iadlp
         iline_stack(2,iline_nstack) = newline

         !xxx  ilines(iadlp)=newline   -- deferred, see add_nl_line in module

         ivrange(:,iadlp)=0
         irrange(:,iadlp)=0
         krepeat(iadlp)=0

         if(iw2.ge.iw1) then
            ivrange(1,iadlp)=iw1+val1-1
            ivrange(2,iadlp)=iw2+val1-1
         endif
      endif
C
      return
      end
