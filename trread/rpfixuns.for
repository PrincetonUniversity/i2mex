      subroutine rpfixuns(units,transform,iok)
C
C  call up an RPLOT units transformation.
C
C  does stuff like #/cm3/sec -> #/sec (e.g. for "VOLINT")
C
C  in/out
      character*(*) units               ! units to be transformed
C
C  in
      character*(*) transform           ! type of transform
C
C  out
      integer iok                       ! =1 if transform type is known.
C
C  if iok.ne.1 on exit, units are unchanged.  no message is generated.
C
C     iok=0 means transform not recognized.
C     iok=2 means transform recognized but units change call did not
C           change the units string
C
C----------------------------
C
      character*20 ztest
C
      character*10 zabbr
C----------------------------
C
      ztest=transform
      call trcaps(ztest)
C
      iok=0
      if(ztest.eq.'VOLINT') then
         iok=1
         iint=1
      else if(ztest.eq.'FLXINT') then
         iok=1
         iint=2
      else if(ztest.eq.'ARINT') then
         iok=1
         iint=3
      else if(ztest.eq.'GRAD') then
         iok=1
         iint=4
      else if(ztest.eq.'LOGDERIV') then
         iok=1
         iint=5
      else if(ztest.eq.'SCALEN') then
         iok=1
         iint=6
      else if(ztest.eq.'TIME_DERIV') then
         iok=1
         iint=9
      endif
C
      if(iok.eq.1) then
         call mglbsw(iint,units,isw)   ! RPLOT units transform
         if(isw.eq.0) iok=2
      endif
C
      return
      end
