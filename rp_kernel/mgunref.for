      subroutine mgunref(zname)
C
      use cplotr_mod
C
      character*(*) zname               ! name of function
C
C  dmc 11 Aug 1999 -- remove any multigraph references to the named
C  function, which is about to be deleted.
C
      character*21 ztest
C-----------------------------------------
C
      ztest=zname
      call uupper(ztest)
C
      do ipkg=1,nbal
         inum=infb(ipkg)
         do ipos=1,inum
            ifun=iabs(ifunb(ipos,ipkg))
            if(iintb(ipkg).eq.1) then
               if(ztest.eq.abt(ifun)) then
                  call mgdelfun(ipkg,ifun,ipos)
                  go to 1000
               endif
            else
               if(ztest.eq.abr(ifun)) then
                  call mgdelfun(ipkg,ifun,ipos)
                  go to 1000
               endif
            endif
         enddo
      enddo
C
 1000 continue
      return
      end
 
