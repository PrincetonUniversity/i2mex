subroutine trx_bxtr(ibccw,ipccw,irzflag,ierr)

  use trx_module
  use trx_bxtr_options
  implicit NONE

!
!  form the B field bicubic splines in the core plasma, and,
!  optionally, extrapolate field to cover (R,Z) rectangular grid.
!
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  integer, intent(inout) :: ibccw ! Btoroidal:  1 for ccw, -1 for cw
  integer, intent(inout) :: ipccw ! Itoroidal:  1 for ccw, -1 for cw

! if the "ccw" variables are 0 on input, try to set them from trdatbuf data...

  integer, intent(in) :: irzflag ! =1:  extrapolate field to (R,Z) region
!   if free boundary data is present, the field is defined over the (R,Z)
!   region, regardless of this switch setting.
!
  integer, intent(out) :: ierr ! completion code, 0=OK
!
!--------------------------------
  external eqm_brz_adhoc,eqm_brzx,eqm_bpsi
!--------------------------------
!
  integer :: lunzer,iwarn,ident=0
  integer :: iccwb,iccwj,ilun,istat
  logical :: ilccw(1),ifreebdy
  integer :: irzflag1
!
!--------------------------------
!
  irzflag1 = irzflag
  if((nRfree.gt.0).and.(nZfree.gt.0)) then
     irzflag1 = 1
     ifreebdy=.TRUE.
  else
     ifreebdy=.FALSE.
  endif

  ierr=0

  if((ibccw.eq.0).or.(ipccw.eq.0)) then
     call trx_trdatbuf_connect(iwarn)
     if(d_data_avail) then
        call tdb_ccwchk(d,iccwb,iccwj)
        if(ibccw.eq.0) then
           ibccw = iccwb
           if(ibccw.eq.0) write(lunzer(0),*) &
                ' ?trx_bxtr: Bphi CCW sign not found in TRANSP archive data.'
        endif
        if(ipccw.eq.0) then
           ipccw = iccwj
           if(ipccw.eq.0) write(lunzer(0),*) &
                ' ?trx_bxtr: Jphi CCW sign not found in TRANSP archive data.'
        endif
     else
        write(lunzer(0),*) ' ?trx_bxtr: TRANSP archive data access failed;'// &
             ' cannot set ibccw or ipccw.'
     endif
  endif

  if((ibccw.ne.1).and.(ibccw.ne.-1)) then
     ierr=ierr+1
     write(lunzer(0),*) ' ?trx_bxtr:  ibccw must be +/-1; ibccw=',ibccw
  endif
  if((ipccw.ne.1).and.(ipccw.ne.-1)) then
     ierr=ierr+1
     write(lunzer(0),*) ' ?trx_bxtr:  ipccw must be +/-1; ipccw=',ipccw
  endif
  if(ierr.gt.0) then
     write(lunzer(0),*) ' %trx_bxtr:  RECOVERY attempt: check namelist.'
     call find_io_unit(ilun)
     call tr_getnl_text(ilun,ierr)
     if(ierr.eq.0) then

        !  NLBCCW
        call tr_getnl_logvec('NLBCCW',ilccw,1,istat)
        if(istat.lt.0) then
           ierr=ierr+1
           write(lunzer(0),*) ' ?trx_bxtr: error accessing NLBCCW.'
        else if(istat.eq.0) then
           write(lunzer(0),*) ' %trx_bxtr: NLBCCW defaulted, .TRUE. assumed.'
           ibccw=1
        else
           if(ilccw(1)) then
              ibccw=1
              write(lunzer(0),*) ' %trx_bxtr: found NLBCCW=.TRUE.'
           else
              ibccw=-1
              write(lunzer(0),*) ' %trx_bxtr: found NLBCCW=.FALSE.'
           endif
        endif

        !  NLJCCW
        call tr_getnl_logvec('NLJCCW',ilccw,1,istat)
        if(istat.lt.0) then
           ierr=ierr+1
           write(lunzer(0),*) ' ?trx_bxtr: error accessing NLJCCW.'
        else if(istat.eq.0) then
           write(lunzer(0),*) ' %trx_bxtr: NLJCCW defaulted, .TRUE. assumed.'
           ipccw=1
        else
           if(ilccw(1)) then
              ipccw=1
              write(lunzer(0),*) ' %trx_bxtr: found NLJCCW=.TRUE.'
           else
              ipccw=-1
              write(lunzer(0),*) ' %trx_bxtr: found NLJCCW=.FALSE.'
           endif
        endif

     endif
  endif
  if(ierr.ne.0) return
!
!  compute interior field
!
  call eqm_bset(ibccw,ipccw)
!
!  computer extrapolated field on (R,Z) grid
!
  if(irzflag1.eq.1) then
     if(ifreebdy) then
        call trx_psirz_load(ierr)
        if(ierr.eq.0) then
           call eqm_brz(eqm_bpsi, edge_smooth ,ierr)
        endif
     else
        if(extrap_method.eq.1) then
           call eqm_brz(eqm_brz_adhoc, edge_smooth ,ierr)
        else
           call eqm_brz(eqm_brzx, edge_smooth ,ierr)
        endif
     endif
  endif
  if(ierr.ne.0) return
!
!  store units
!
  call eq_gfnum('Bmod',ident)
  call trx_ustore(ierr,ident,0,-99,'T',0)
  if(ierr.ne.0) return
!
  call eq_gfnum('BR',ident)
  call trx_ustore(ierr,ident,0,-99,'T',0)
  if(ierr.ne.0) return
!
  call eq_gfnum('BZ',ident)
  call trx_ustore(ierr,ident,0,-99,'T',0)
  if(ierr.ne.0) return
!
  if(irzflag1.eq.1) then
     call eq_gfnum('Bphi',ident)
     call trx_ustore(ierr,ident,0,-99,'T',0)
     if(ierr.ne.0) return
  endif
!
  return
!
end subroutine trx_bxtr
