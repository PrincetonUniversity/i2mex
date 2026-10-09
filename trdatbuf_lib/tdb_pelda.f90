subroutine tdb_npel(d,inpel)
  use trdatbuf_obj
  implicit NONE
  !  return number of pellet events
  type (trdatbuf) :: d
  integer, intent(out) :: inpel ! number of pellets this run

  inpel = d%npelda

end subroutine tdb_npel

subroutine tdb_pelda(d,ipel,p,ierr)
  use trdatbuf_obj
  use trdatbuf_aux
  implicit NONE
  !  return data associated with a particular pellet
  type (trdatbuf) :: d
  integer, intent(in) :: ipel   ! index of desired pellet data (ipel'th pellet)
  type (pelget) :: p            ! pellet data returned here
  integer, intent(out) :: ierr  ! return .ne.0 if ipel input out of range

  !-------------------------------------
  integer :: iadr,ifield,ifields
  !-------------------------------------

  if((ipel.le.0).or.(ipel.gt.d%npelda)) then
     ierr=1
     return
  endif

  ierr=0

  iadr=d%lpelda
  ifields = d%datbuf(iadr) + 0.1

  ifield = 0

  iadr=d%lpelda + ipel - d%npelda
  do
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     iadr = iadr + d%npelda
     p%tpel = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%apel = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%apel2 = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%fpel2 = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%nmpela = d%datbuf(iadr) + 0.1

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%plrsta = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%plysta = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%plphia = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%plthea = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%pelrad = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%pelvel = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%kpellet = d%datbuf(iadr) + 0.1

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%tpele = d%datbuf(iadr)

     iadr = iadr + d%npelda
     ifield = ifield+1
     if(ifield.gt.ifields) exit
     p%freqpel = d%datbuf(iadr)

     exit
  enddo

end subroutine tdb_pelda
