      subroutine mg_fixref(itype,ifnum)
C
      use cplotr_mod
C
C  fix MG references after deletion of function #ifnum
C  --all references to function #s > ifnum have to be decremented
C
      integer, intent(in) :: itype  ! 0=profile, 1=scalar
      integer, intent(in) :: ifnum  ! deleted function #
C
C  dmc 11 Aug 1999 -- remove any multigraph references to the named
C  function, which is about to be deleted.
C
C-----------------------------------------
C
      do ipkg=1,nbal
         if(iintb(ipkg).ne.itype) cycle
         inum=infb(ipkg)
         do ipos=1,inum
            ifun=iabs(ifunb(ipos,ipkg))   ! remove sign
            if(ifun.gt.ifnum) then
               if(ifun.ne.ifunb(ipos,ipkg)) then
                  ! -sign is preserved
                  ifun=ifun-1
                  ifunb(ipos,ipkg)=-ifun
               else
                  ! +sign is preserved
                  ifunb(ipos,ipkg)=ifun-1
               endif
            endif
         enddo
      enddo

      return
      end
 
