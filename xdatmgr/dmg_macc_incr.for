      subroutine dmg_macc_incr

      ! increment MACC -- but, if it much exceeds MAXENT, do a cleanup

      use datmgr_mod

      implicit NONE

      integer, dimension(:), allocatable :: lacc_from, lacc_to
      integer :: ilim,imax,ii,jj,ict
      !--------------------------------------

      MACC = MACC + 1

      ilim = 3*MAXENT
      if(MACC.lt.ilim) return

      ! OK, compress the MACC references
      ! this code takes a range of scattered integers from [1:ilim]
      ! and while maintaining ordering reduces them into the range
      ! 1:[# of non-zero entries]

      allocate(lacc_from(ilim),lacc_to(ilim))
      lacc_from = 0
      lacc_to = 0

      jj = 1
      do
         if(lacc(jj).gt.0) then
            lacc_from(lacc(jj)) = jj
         endif

         if(dmglbl(jj).eq.'%FINI') exit
         jj = lnext(jj)
         if(jj.eq.0) exit
      enddo

      ict=0
      do ii=1,ilim
         if(lacc_from(ii).gt.0) then
            ict=ict+1
            lacc_to(ii)=ict
         endif
      enddo

      do ii=1,ilim
         if(lacc_to(ii).gt.0) then
            jj=lacc_from(ii)
            lacc(jj)=lacc_to(ii)
         endif
      enddo

      MACC = ict + 1

      deallocate(lacc_from,lacc_to)

#ifdef __DEBUG
      call dmprin('dmg_macc_incr',1000)
#endif

      return
      end
