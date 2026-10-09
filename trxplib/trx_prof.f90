!  trx_prof -- get conversion factor, read & scale profile
!  trx_kprof -- read & scale profile
!
!--------------------------------------------------------------------
!  get units conversion profile, then, read & scale profile (trx_kprof,
!  below)
!
subroutine trx_prof(zname,zuns_mks,iordr,ident,ierr)
 
  use trx_module
  implicit NONE
 
! ...arguments
 
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zname          ! name of desired profile
  character(*), intent(out) :: zuns_mks      ! mks units label of profile
  integer, intent(in)      :: iordr          ! desired order of fit
  integer, intent(out)     :: ident          ! xplasma object id#
  integer, intent(out)     :: ierr           ! completion code, 0=OK
 
! ...local
 
  real*8                   :: zconv          ! units conversion factor
  integer                  :: iwarn          ! units conversion warning
  integer                  :: ixid           ! x axis id
 
! ------------------------------
 
  call trx_mks_conv(zname,zconv,zuns_mks,iwarn)
  call trx_kprof(zname,iordr,zconv,ident,ixid,ierr)
  call trx_ustore(ierr,ident,ixid,iordr,zuns_mks,iwarn)
 
  return
end subroutine trx_prof
!--------------------------------------------------------------------
!  read a TRANSP profile & interpolate/extrapolate to flux surfaces
!  from zone centers if necessary; scale by units conversion factor
!
subroutine trx_kprof(zname,iordr,zconv,ident,ixid,ierr)
 
  use trx_module
  implicit NONE
 
! ...arguments
 
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zname          ! name of desired profile
  integer, intent(in)      :: iordr          ! desired order of fit
  real*8, intent(in)       :: zconv          ! units conversion factor
  integer, intent(out)     :: ident          ! xplasma object id#
  integer, intent(out)     :: ixid           ! xplasma axis id#
  integer, intent(out)     :: ierr           ! completion code, 0=OK
 
! ...local
 
  real, dimension(nxmax) :: zbuf     ! data read in
  real*8, dimension(nxmax) :: zbuf1  ! converted to real*8
!
  real*8, dimension(nsurf+1) :: zbuf2  ! interp./extrapolated to flux surfaces
  real*8 :: zdum
!
  character(64) :: zlabel            ! label
  character(32) :: zunits            ! units
!
  integer :: istype,ixgot            ! type code, #gotten
  integer :: lunzer
  integer :: iordri,iuser,iclass,idnum_rp
!
! ------------------------------
!
  ierr=0
  ixid=0
  ident=0
!
  call trx_ready('trx_prof',ierr)
  if(ierr.ne.0) return
!
  call eq_gfnum(zname,ident)
  iordri = iordr
  if(ident.ne.0) then
     if(prof_iord(ident).eq.-99) then
        write(lunzer(0),*) '?trx_prof: ',zname,' reserved by TRXPLIB.'
        ierr=1
        return
     endif
     call idchek(zname,idnum_rp,iclass,iuser,0)
     if(iuser.eq.0) then
        return  ! TRREAD immutable object; we have it already...
     else
        iordri=100+iordr  ! replaceable object -- but does type match?
     endif
  endif
!
  call t1profil(zname,zlabel,zunits,time0,delta_t, &
       istype,zbuf,nxmax,ixgot, ierr)
  if(ierr.ne.0) return
!
  zdum=0
!
  if((istype.ne.itype_x).and.(istype.ne.itype_xb).and. &
       (istype.ne.itype_rmajm).and.(istype.ne.itype_rmjsym)) then
     write(lunzer(0),*) '?trx_prof:  ',zname,':  ', &
          'neither a function of x (flux zones or surfaces)'
     write(lunzer(0),*) '            nor a function of major radius.'
     ierr=1
     return
  endif
!
  if((istype.eq.itype_x).or.(istype.eq.itype_xb)) then
     if(ixgot.ne.(nsurf-1)) then
        write(lunzer(0),*) '?trx_prof:  ',zname,':  ixgot+1 .ne. nsurf'
        write(lunzer(0),*) ' (this is unexpected) ixgot=',ixgot,' nsurf=',nsurf
        ierr=1
        return
     endif
  else if(istype.eq.itype_rmajm) then
     if(ixgot.ne.nrmajm) then
        write(lunzer(0),*) '?trx_prof: f(Rmajm) no. of points ixgot=',ixgot
        write(lunzer(0),*) ' expected: NRmajm = ',NRmajm
        ierr=2
        return
     endif
  else if(istype.eq.itype_rmjsym) then
     if(ixgot.ne.nrmjsym) then
        write(lunzer(0),*) '?trx_prof: f(Rmjsym) no. of points ixgot=',ixgot
        write(lunzer(0),*) ' expected: NRmjsym = ',NRmjsym
        ierr=3
        return
     endif
  endif
!
!  OK
!
  zbuf1(1:ixgot) = zbuf(1:ixgot)*zconv
!
  if(istype.eq.itype_xb) then
!
!  copy & extrapolate axial value
!
     zbuf2(2:nsurf)=zbuf1(1:nsurf-1)
     if(ival_axis.eq.1) then
        zbuf2(1)=zval_axis
     else
        if((ibc_axis.eq.1).and.(zbc_axis.eq.0.0E0_R8)) then
!
!  assuming f' -> 0 on axis
!
           zbuf2(1) = zbuf1(1) - (zbuf1(2)-zbuf1(1))*afac0b
        else
!
!  linear extrapolation
!
           zbuf2(1) = zbuf1(1) - (zbuf1(2)-zbuf1(1))*afaclinb
!
        endif
!  extrapolation limit check
        if(isign_check.eq.1) then
           if(zbuf1(1).lt.0.0E0_R8) then
              zbuf2(1)=min(zsign_lim*zbuf1(1),zbuf2(1))
           else if(zbuf1(1).gt.0.0E0_R8) then
              zbuf2(1)=max(zsign_lim*zbuf1(1),zbuf2(1))
           endif
        endif
!
     endif
!
  else if(istype.eq.itype_x) then
!
!  take zone centered data; extrapolate at edges.
!
     zbuf2(2:nsurf)=zbuf1(1:nsurf-1)
!
     if(ival_axis.eq.1) then
        zbuf2(1)=zval_axis
     else
        if((ibc_axis.eq.1).and.(zbc_axis.eq.0.0E0_R8)) then
!
!  assuming f' -> 0 on axis
!
           zbuf2(1) = zbuf1(1) - (zbuf1(2)-zbuf1(1))*afac0
        else
!
!  linear extrapolation
!
           zbuf2(1) = zbuf1(1) - (zbuf1(2)-zbuf1(1))*afaclin
!
        endif
!  extrapolation limit check
        if(isign_check.eq.1) then
           if(zbuf1(1).lt.0.0E0_R8) then
              zbuf2(1)=min(zsign_lim*zbuf1(1),zbuf2(1))
           else if(zbuf1(1).gt.0.0E0_R8) then
              zbuf2(1)=max(zsign_lim*zbuf1(1),zbuf2(1))
           endif
        endif
!
     endif
!
! now... do boundary ... note use of zbuf2(nsurf+1);
!  zbuf2(2:nsurf) will be reset to the original, uninterpolated
!  data values...
!
     if(ival_edge.eq.1) then
        zbuf2(nsurf+1)=zval_edge
     else
        if((ibc_edge.eq.1).and.(zbc_edge.eq.0.0E0_R8)) then
!
!  assuming f' -> 0 on edge
!
           zbuf2(nsurf+1) = &
                zbuf1(nsurf-1) - (zbuf1(nsurf-2)-zbuf1(nsurf-1))*efac0
        else
!
!  linear extrapolation
!
           zbuf2(nsurf+1) = &
                zbuf1(nsurf-1) - (zbuf1(nsurf-2)-zbuf1(nsurf-1))*efaclin
!
        endif
!  extrapolation limit check
        if(isign_check.eq.1) then
           if(zbuf1(nsurf-1).lt.0.0E0_R8) then
              zbuf2(nsurf+1)=min(zsign_lim*zbuf1(nsurf-1),zbuf2(nsurf+1))
           else if(zbuf1(nsurf-1).gt.0.0E0_R8) then
              zbuf2(nsurf+1)=max(zsign_lim*zbuf1(nsurf-1),zbuf2(nsurf+1))
           endif
        endif
     endif
!
  endif       ! zone ctr / zone bdy
!
!  OK... set up the interpolation object
!
  if(istype.eq.itype_xb) then
     ixid=id_rho
     call eqm_rhofun(iordri,id_rho,zname,zbuf2, &
          ibc_axis,zbc_axis,ibc_edge,zbc_edge, &
          ident,ierr)
 
  else if(istype.eq.itype_x) then
     ixid=id_rhozc
     call eqm_rhofun(iordri,id_rhozc,zname,zbuf2, &
          ibc_axis,zbc_axis,ibc_edge,zbc_edge, &
          ident,ierr)
 
  else if(istype.eq.itype_rmajm) then
     ixid=id_rmajm
     call eqm_f1d(iordri,id_rmajm,zname,zbuf1(1:nrmajm), &
          0,zdum,0,zdum, ident,ierr)
 
  else if(istype.eq.itype_rmjsym) then
     ixid=id_rmjsym
     call eqm_f1d(iordri,id_rmjsym,zname,zbuf1(1:nrmjsym), &
          0,zdum,0,zdum, ident,ierr)
 
  endif
!
!  re-init BC flags
!
  call trx_bc_init

end subroutine trx_kprof
