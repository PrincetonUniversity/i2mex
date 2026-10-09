      subroutine mgtryadd(ipkg,if)
C
      use cplotr_mod
C
C  plcexec subroutine (MG_ADDFUN and MG_CREATE commands use this).
C  look at argument #if.  If not defaulted, look for signed function id.
C  if valid after checks, add to multigraph package.
C
      character*21 zid
      character*22 zidarg
C
      call plcgarg(if,zidarg)
      if(zidarg.eq.'%empty') return
C
      ifin=len_trim(zidarg)
      istz=1
      isign=1
C
C  check for signed argument
      if(zidarg(1:1).eq.'+') istz=istz+1
      if(zidarg(1:1).eq.'-') then
         istz=istz+1
         isign=-1
      endif
C
      zid=zidarg(istz:ifin)
      if(iintb(ipkg).eq.0) then
C  expect profile
         ifun=ifind_ordr(abr,iordrr,nfxt,zid)
         if(ifun.eq.0) then
            call zermsg(' %plcexec:  not a profile function:  '//zid)
         endif
      else
C  expect scalar
         ifun=ifind_ordr(abt,iordrt,nft,zid)
         if(ifun.eq.0) then
            call zermsg(' %plcexec:  not a profile function:  '//zid)
         endif
      endif
      if(infb(ipkg).eq.naxmgf) then
         call zermsg(' %plcexec: cannot add: '//trim(zid)//
     >        ' multigraph is full.')
         ifun=0
      endif
      if(ifun.gt.0) then
C  check against duplication
         idup=0
         do ifck=1,infb(ipkg)
            ifunck=abs(ifunb(ifck,ipkg))
            if(ifun.eq.ifunck) idup=1
         enddo
         if(idup.gt.0) then
            call zermsg(
     >         ' %plcexec: duplicate, not added to multigraph:  '//zid)
         else
C  add function into the packages
            infb(ipkg)=infb(ipkg)+1
            ifunb(infb(ipkg),ipkg)=isign*ifun
         endif
      endif
C
      return
      end
