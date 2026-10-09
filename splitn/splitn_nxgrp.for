      subroutine splitn_nxgrp(zline,iprev,iw1,iw2,idelim,
     >   ir1,ir2,irepeat,iextra)
C
C find next group of values: could be one value followed by a
C comma or end of line, or just the end of line (a null value),
C or nn*(a value) or nn* & end of line (nn null values).
C
C  input
C
      character*(*) zline               ! namelist line input
      integer iprev                     ! start search at: "=" or prev. ",".
C                                         or whitespace delimiter
C  output
C    iw1:iw2 -- range of characters giving value field (null value: iw2<iw1);
C    ir1:ir2 -- range of characters giving repeat count field, with "*",
C               or 0:0.
C
      integer iw1,iw2                   ! delimit found value
C         iw2.lt.iw1 on output means:  null value
C
      integer idelim                    ! terminating delimiter / EOL
C
      integer ir1,ir2                   ! repeat count found (both 0 if none)
      integer irepeat                   ! repeat count
C
      integer iextra                    ! =0: OK, .gt.0: extraneous chars
C
C  iextra is set if there are extraneous characters after the value
C  but before the delimiter.
C
C--------------------------------------
C
      character*10 zdigs
C
      data zdigs/'0123456789'/
C--------------------------------------
C
C  the following is the answer if a null value is detected
C
      iextra=0
      iw1=iprev+1
      iw2=iprev
      ir1=0
      ir2=0
      irepeat=1
      ilen=len(zline)
      idelim=ilen
C
      ic=iprev
C
C  find first non-white-space character
C
 10   continue
      ic=ic+1
      if(ic.gt.idelim) then
C  EOL, all blank
         return
      endif
      if(zline(ic:ic).eq.' ') go to 10
      if(zline(ic:ic).eq.char(9)) go to 10
C
C  non-blank found...
C
      if(zline(ic:ic).eq.',') then
         idelim=ic
         return                         ! null value
      endif
C
C  non-null... digit?
C
      idig=index(zdigs,zline(ic:ic))
      if(idig.eq.0) go to 50
C
C  yes, digit(s)...
      ic0=ic
      icount=idig-1
C
C  looking at an integer...
C
 20   continue
      ic=ic+1
      if(ic.gt.ilen) then
         iw1=ic0
         iw2=ilen
         idelim=ilen
         irepeat=1
         return                         ! a single integer
      endif
C
      idig=index(zdigs,zline(ic:ic))
      if(idig.gt.0) then
         icount=10*icount+(idig-1)
         go to 20
      endif
C
C  not a digit.  Look for first nonblank.
C
      ic1=ic-1
      ic=ic1
 30   continue
      ic=ic+1
      if(ic.gt.ilen) then
         iw1=ic0
         iw2=ic1
         idelim=ilen
         irepeat=1
         return                         ! a single integer
      endif
      if(zline(ic:ic).eq.' ') go to 30
      if(zline(ic:ic).eq.char(9)) go to 30
C
C  a non-blank.  comma is a delimiter.
C
      if(zline(ic:ic).eq.',') then
         iw1=ic0
         iw2=ic1
         idelim=ic
         irepeat=1
         return                         ! a single integer
      endif
C
C  a "*" means a repeat count
C
      if(zline(ic:ic).eq.'*') then
         ir1=ic0
         ir2=ic
         irepeat=icount
C
C  look for next non-blank
 40      continue
         ic=ic+1
         if(ic.gt.ilen) then
            return                      ! repeated nulls
         endif
         if(zline(ic:ic).eq.' ') go to 40
         if(zline(ic:ic).eq.char(9)) go to 40
C  nonblank
         if(zline(ic:ic).eq.',') then
            idelim=ic
            return                      ! repeated nulls
         endif
C  nonblank, noncomma -- a value
         go to 50
      else
C
C  it was not a repeat count; must have been a value of some kind
C
         ic=ic0
         go to 50
      endif
C
C------------------------------
C  we have a value starting at ic
C
 50   continue
      iw1=ic
      if(zline(ic:ic).eq.'"') go to 80 ! quote string
      if(zline(ic:ic).eq."'") go to 80 ! quote string
C
C  non quoted value ends at next non-blank
C
 60   continue
      ic=ic+1
      if(ic.gt.ilen) then
         iw2=ilen
         idelim=ilen
         return
      endif
      if(zline(ic:ic).eq.',') then
         iw2=ic-1
         idelim=ic
         return
      endif
      if((zline(ic:ic).eq.' ').or.(zline(ic:ic).eq.char(9))) then
         iw2=ic-1
         go to 100
      endif
      go to 60
C
C  quote string
C
 80   continue
C
C  look for terminating quote
C
      ic=ic+1
      if(ic.gt.ilen) then
         iw2=ilen
         idelim=ilen
         return
      endif
      if(zline(ic:ic).eq.'''') then
         if(ic.lt.ilen) then
            if(zline(ic+1:ic+1).eq.'''') then
               ic=ic+1                  ! paired quotes '' -> literal '
               go to 80
            endif
         endif
      else
         go to 80                       ! not a quote
      endif
C
C  have end of quote string-- must be delimited by whitespace or
C  comma or EOL
C
      iw2=ic
      if(iw2.ge.ilen) go to 100
      ii=ic+1
      if((zline(ii:ii).ne.',').and.(zline(ii:ii).ne.' ').and.
     >   (zline(ii:ii).ne.char(9))) then
         iextra=1                       ! bogus extra chars after quote string
      endif
      go to 100
C
 100  continue
C
C  now look for terminating delimiter
C
      ic=ic+1
      if(ic.gt.ilen) then
         idelim=ilen
         return
      endif
C
      if(zline(ic:ic).eq.',') then
         idelim=ic
         return
      endif
      if(zline(ic:ic).eq.' ') go to 100
      if(zline(ic:ic).eq.char(9)) go to 100
C
C  non-blank non comma:
C  assume whitespace delimiter
C
      idelim=ic-1
      return
C
      end
 
