subroutine dmg_wkxpand(idiff,ier)

  ! dmg_wkxpand
  !  (try to) make space in memory for expansion of low-end buffers by
  !  the amount (idiff).  If NO_DELETE is set, try hard.

  use datmgr_mod
  implicit NONE

  integer, intent(in) :: idiff   ! amount to expand
  integer, intent(out) :: ier    ! status code on exit, 0=OK

  !---------------------------------------
  !  %F(T) and at least one other low-end named buffer is expected to exist.
  !  zero or more high end buffers may also be present
  !---------------------------------------

  integer :: iadr_lob  ! address of last word of current low end buffer data
  integer :: iadr_hia  ! address of 1st word of high end buffer data

  integer :: iadr_mida ! address of 1st word of midsection buffer data
  integer :: iadr_midb ! address of last word of midsection buffer data

  integer :: ict_lo    ! number of low end buffer blocks found (.gt.1 expected)

  integer :: j,jj,jsave,ifree,ihilo,ifree_hi

  integer :: iadr1,iadr2
  integer :: istart,istartt,istartlo,iend
  integer :: inum,itot,ia,inew

  integer, dimension(:), allocatable :: idelm,idell

  logical :: idlock
  !---------------------------------------
  
  do
     ! initialize...

     ier = 0

     iadr_lob = 0
     iadr_mida = ndbsiz + 1
     iadr_midb = 0
     iadr_hia = ndbsiz + 1

     ict_lo = 0

     istartlo = ndbsiz + 1

     ! scan the blocks...

     do j=1,ndent
        if(nwds(j).le.0) cycle

        if(mprio(j).eq.lo_end_prio) then
           ict_lo = ict_lo + 1
           iadr_lob = max(iadr_lob, locd(j)+nwds(j)-1 )
           if(dmglbl(j).ne.'%F(T)') then
              istartlo=min(istartlo,locd(j))
           else
              istartt=locd(j)
           endif

        else if(mprio(j).eq.hi_end_prio) then
           iadr_hia = min(iadr_hia, locd(j))

        else
           iadr_mida = min(iadr_mida, locd(j))
           iadr2 = locd(j)+nwds(j)-1
           if(iadr2.gt.iadr_midb) then
              jsave = j
              iadr_midb = iadr2
           endif
        endif
     enddo

     ! sanity checks...

     if(ict_lo.lt.2) then
        write(lundmo,*) &
           ' ?dmg_wkxpand unexpected: .lt.2 low end memory blocks.'
        ier = ier + 1
     endif

     if(iadr_lob.ge.iadr_mida) then
        write(lundmo,*) &
           ' ?dmg_wkxpand unexpected: low and midrange memory blocks overlap.'
        ier = ier + 1
     endif

     if(iadr_midb.ge.iadr_hia) then
        write(lundmo,*) &
           ' ?dmg_wkxpand unexpected: midrange and high memory blocks overlap.'
        ier = ier + 1
     endif

     if(istartt.gt.istartlo) then
        write(lundmo,*) &
           ' ?dmg_wkxpand unexpected: %F(T) not the 1st low end block.'
        ier = ier + 1
     endif

     if(ier.gt.0) then
        write(lundmo,*) &
           ' ?dmg_wkxpand(xdatmgr): error in DATBUF(...) entries.'
        exit
     endif

     !----------------------------------
     ! OK

     ifree = min(iadr_mida,iadr_hia) - iadr_lob - 1

     ! space btw low end and high end blocks...

     ihilo = iadr_hia - iadr_lob - 1
     if(ihilo.lt.idiff) then
        ! more space definitely needed
        call dmg_datbuf_expand(0)
        cycle
     endif

     ! OK if we get this far there may be midrange blocks "in the way..."
     ! here is the space that needs to be cleared:

     ifree_hi = iadr_hia - iadr_midb - 1

     iadr1 = iadr_lob + 1
     iadr2 = iadr_lob + idiff

     ! add up the space in the blocks that need moving

     inum = 0  ! number of blocks
     itot = 0  ! sum of space in blocks

     if(allocated(idell)) deallocate(idell)
     allocate(idell(ict_lo-1))
     ict_lo = 0

     do j=1,ndent
        if(nwds(j).le.0) cycle
        if(mprio(j).eq.lo_end_prio) then
           ict_lo=ict_lo + 1
           if(ict_lo.gt.1) then
              idell(ict_lo-1)=j   ! low end list
           endif
           cycle
        endif
        if(mprio(j).eq.hi_end_prio) cycle

        istart = locd(j)

        if(istart.le.iadr2) then
           inum = inum + 1
           itot = itot + nwds(j)
        endif

     enddo

     if(ifree_hi.lt.itot) then
        ! apparent insufficient space to which to move the (inum) blocks...
        ! more space definitely needed
        call dmg_datbuf_expand(0)
        cycle
     endif

     ! So, move the "in-the-way" blocks...

     idlock = no_delete
     no_delete = .FALSE.

     if(allocated(idelm)) deallocate(idelm)
     allocate(idelm(inum))
     inum = 0

     do j=1,ndent
        if(nwds(j).le.0) cycle
        if(mprio(j).eq.lo_end_prio) cycle
        if(mprio(j).eq.hi_end_prio) cycle

        istart = locd(j)
        iend = locd(j)+nwds(j)-1

        if(istart.le.iadr2) then
           call dminew(jsave,jj,mprio(j))
           nwds(jj)=nwds(j)
           lacc(jj)=lacc(j)

           ! jsave was the highest midrange block, but make
           ! sure it was high enough

           locd(jj)=max(locd(jj),iadr2+1)

           ! copy the data
           do ia=istart,iend
              inew = locd(jj) + (ia-istart)
              datbuf(inew)=datbuf(ia)
           enddo

           dmglbl(jj)=dmglbl(j)
           dmglbl(j)='(deleted)'

           ! add to deletions list

           inum = inum + 1
           idelm(inum) = j

           jsave = jj  ! save the new one, tacked on now...
        endif
     enddo

     ! now delete the copied blocks from the middle section...
     do ia=1,inum
        call dmidel(idelm(ia))
     enddo

     ! now copy the low section blocks' data, adjust addresses
     do ia=1,ict_lo-1
        locd(idell(ia)) = locd(idell(ia))+idiff
     enddo

     do ia=iadr_lob,istartlo,-1
        datbuf(ia+idiff)=datbuf(ia)
     enddo

     no_delete = idlock

     deallocate(idell,idelm)

     exit  ! all done...
  enddo

end subroutine dmg_wkxpand
