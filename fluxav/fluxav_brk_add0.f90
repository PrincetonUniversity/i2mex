subroutine fluxav_brk_addperio(oldlist,n_old,newlist,n_new,mergelist,n_merge)
!
!  merge two lists in ascending order, squeezing out duplicates, with only
!  the subset of newlist that is between oldlist(1) and oldlist(n_old)
!  considered.
!
!  presume a periodic coordinate, i.e. newlist can be rotated, and shifted
!  by integer multiples of 2pi
!
!    **newlist** should cover a range of 2pi
!
!  input lists assumed already in ascending order *without checking*
!
  implicit NONE
!
!   input:
!
  integer n_old
  real*8 oldlist(n_old)         ! old list
!
  integer n_new
  real*8 newlist(n_new)         ! new list
!
!   output:
!
  real*8 mergelist(*)           ! merged list
  integer n_merge               ! length of merged list
!
!--------------------------------------
  real*8 :: snewlist(n_new),rnewlist(n_new)
  integer :: ishift
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  integer :: i,ifound,inf
!--------------------------------------

  if(oldlist(1).lt.newlist(1)) then
     ishift=1 + (newlist(1)-oldlist(1))/c2pi
     snewlist = newlist - ishift*c2pi
  else if(oldlist(1).ge.newlist(n_new)) then
     ishift=1 + (oldlist(1)-newlist(n_new))/c2pi
     snewlist = newlist + ishift*c2pi
  else
     snewlist = newlist
  endif

  ! 1 correction allowed

  if(oldlist(1).lt.snewlist(1)) then
     snewlist = snewlist-c2pi
  else if(oldlist(1).ge.snewlist(n_new)) then
     snewlist = snewlist+c2pi
  endif

  if((oldlist(n_old).le.snewlist(1)).or.(oldlist(1).ge.snewlist(n_new))) then
     n_merge=n_old
     mergelist(1:n_old)=oldlist
     return
  endif

  ! OK range lf lists intersects; find first snewlist element that is
  ! greater than the first oldlist element

  ifound=0
  do i=1,n_new
     if(snewlist(i).gt.oldlist(1)) then
        ifound=i
        exit
     endif
  enddo

  if(ifound.eq.0) then
     n_merge=n_old
     mergelist(1:n_old)=oldlist
     return
  endif

  inf=n_new-ifound+1
  rnewlist(1:inf)=snewlist(ifound:n_new)
  rnewlist(inf+1:n_new-1)=snewlist(2:ifound-1)+c2pi
  rnewlist(n_new)=rnewlist(1)+c2pi

  ! OK now do the merge...

  call fluxav_brk_addsub(oldlist,n_old,rnewlist,n_new,mergelist,n_merge)

end subroutine fluxav_brk_addperio
!-------------------------------------------------------------------
subroutine fluxav_brk_addsub(oldlist,n_old,newlist,n_new,mergelist,n_merge)
!
!  merge two lists in ascending order, squeezing out duplicates, with only
!  the subset of newlist that is between oldlist(1) and oldlist(n_old)
!  considered.
!
!  input lists assumed already in ascending order *without checking*
!
  implicit NONE
!
!   input:
!
  integer n_old
  real*8 oldlist(n_old)         ! old list
!
  integer n_new
  real*8 newlist(n_new)         ! new list
!
!   output:
!
  real*8 mergelist(*)           ! merged list
  integer n_merge               ! length of merged list
!
!--------------------------------------
!
  integer i,i1,i2,inum
!
  i1=0
  i2=0
!
  do i=1,n_new
     if(newlist(i).gt.oldlist(1)) then
        if(i1.eq.0) i1=i
     endif
     if(newlist(i).lt.oldlist(n_old)) then
        i2=i
     endif
  enddo

  if((min(i1,i2).eq.0).or.(i2.lt.i1)) then
     mergelist(1:n_old)=oldlist(1:n_old)
     n_merge = n_old
  else

     inum=i2-i1+1
     call fluxav_brk_add0(oldlist,n_old,newlist(i1:i2),inum, &
          mergelist,n_merge)

  endif
end subroutine fluxav_brk_addsub
!--------------------------------------
!-------------------------------------------------------------------
subroutine fluxav_brk_add0(oldlist,n_old,newlist,n_new,mergelist,n_merge)
!
!  merge two lists in ascending order, squeezing out duplicates
!  input lists assumed already in ascending order *without checking*
!
  implicit NONE
!
!   input:
!
  integer n_old
  real*8 oldlist(n_old)         ! old list
!
  integer n_new
  real*8 newlist(n_new)         ! new list
!
!   output:
!
  real*8 mergelist(*)           ! merged list
  integer n_merge               ! length of merged list
!
!--------------------------------------
!
  integer i_old,i_new
  real*8 znew
!
!--------------------------------------
!
  i_old=0
  i_new=0
!
  n_merge=0
!
!  loop until ascending lists are merged into one ascending list.
!
  do
!
!  loop exit test...
!
     if((i_old.ge.n_old).and.(i_new.ge.n_new)) then
        return
!
     else if(i_old.ge.n_old) then
        i_new=i_new+1
        n_merge=n_merge+1
        mergelist(n_merge)=newlist(i_new)
!
     else if(i_new.ge.n_new) then
        i_old=i_old+1
        n_merge=n_merge+1
        mergelist(n_merge)=oldlist(i_old)
!
     else
!
!   both lists are in play
!
        if(oldlist(i_old+1).lt.newlist(i_new+1)) then
           i_old=i_old+1
           znew=oldlist(i_old)
        else
           i_new=i_new+1
           znew=newlist(i_new)
        endif
!
        if(n_merge.eq.0) then
           n_merge=1
           mergelist(n_merge)=znew
!
!  no duplicates...
!
        else if(znew.gt.mergelist(n_merge)) then
           n_merge=n_merge+1
           mergelist(n_merge)=znew
        endif
!
     endif
  enddo
!
end subroutine fluxav_brk_add0

subroutine fluxav_brktol(list,nlist,ztol)

  ! squeeze out successive entries of list that are within ztol of eachother
  ! list must be in ascending order, **not checked**

  real*8 :: list(*)  ! modified on output
  integer :: nlist   ! modified on output
  real*8,intent(in) :: ztol

  !---------------
  real*8 :: ztoli
  integer :: i,ii
  real*8, parameter :: czero = 0.0d0
  real*8, parameter :: ceps4 = 1.0d-4
  !---------------

  ! sanity check...

  if(ztol.le.czero) return
  ztoli=min(ztol,max(abs(list(1)),abs(list(nlist)))*ceps4)

  !---------------

  ii=1
  do i=2,nlist
     if((list(i)-list(i-1)).gt.ztoli) then
        ii=ii+1
        if(ii.lt.i) list(ii)=list(i)
     endif
  enddo

  if(ii.lt.nlist) list(ii+1:nlist)=0
  nlist=ii

end subroutine fluxav_brktol
