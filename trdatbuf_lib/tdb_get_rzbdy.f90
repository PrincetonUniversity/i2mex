subroutine tdb_get_rzbdy(d,ztime,nbdy,rbdy,zbdy,ier,zmid,r1,r2)

  !  Get R,Z coords along boundary, either from UFILES or IMAS

  use iso_c_binding, only: rp => c_double
  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE
  integer, external :: lunmsg_tdb

  type (trdatbuf), intent(in) :: d
  real(rp), intent(in) :: ztime       ! time to which to interpolate
  integer, intent(in) :: nbdy         ! size of bdry coord arrays
  real(rp), dimension(nbdy), intent(out) :: rbdy,zbdy ! bndry pt coords
  integer, intent(out) :: ier         ! status code: 0=OK
  real(rp), intent(out), optional :: zmid    ! Midplane height
  real(rp), intent(out), optional :: r1, r2  ! Boundary-midplane intersections
  !------------------------------------
  real(rp), dimension(2) :: rmid
  real(rp) :: zfrac1,zfrac2
  integer  :: nonlin,rzsize,ia,iloct,inumt,it1,jb
  !------------------------------------

  ier = 0

  ! Error check
  nonlin = lunmsg_tdb(0)

  if (nbdy.ne.d%nbdy) then
     write(nonlin,*)'Boundary array size mismatch in tdb_get_rzbdy.'
     ier = 1
     return
  endif

  if (nbdy.le.0) then
     write(nonlin,*) 'tdb_get_rzbdy: no boundary coordinate data available.'
     ier = 2
     return
  endif

  iloct = d%ltime2;  inumt = d%ntime2

  if((iloct.eq.0).or.(inumt.eq.0)) then
     write(nonlin,*) '  R_bdy(t,j) timebase address: ',iloct
     write(nonlin,*) '             timebase size:    ',inumt
     write(nonlin,*) '  --> zero not expected.'
     ier=3
     return
  endif

  ! Find closest time slices to ztime, interpolation coefficients
  call tdbsub_lookup(d%datbuf(iloct:iloct+inumt-1),inumt,ztime,it1,zfrac2)
  zfrac1 = 1.0_rp - zfrac2

  ! Interpolate R_boundary from trdatbuf data
  ia = d%ldrbdy + it1 - 1
  do jb=1,nbdy
     rbdy(jb) = zfrac1*d%datbuf(ia) + zfrac2*d%datbuf(ia+1)
     ia = ia + inumt
  enddo

  ! Interpolate Z_boundary from trdatbuf data
  ia = d%ldzbdy + it1 - 1
  do jb=1,nbdy
     zbdy(jb) = zfrac1*d%datbuf(ia) + zfrac2*d%datbuf(ia+1)
     ia = ia + inumt
  enddo

  if (present(zmid)) then
     ! Find z at midplane
     ia = d%ldatzpl + it1 - 1
     zmid = zfrac1*d%datbuf(ia) + zfrac2*d%datbuf(ia+1)

     if (present(r1).and.present(r2)) then
        ! Find intersections of boundary with midplane
        inumt = 1;  jb = 2
        do while (jb.le.nbdy)
           if ((zbdy(jb)-zmid)*(zbdy(jb-1)-zmid).le.0.0) then
              rmid(inumt) = rbdy(jb-1) + &
                   (rbdy(jb)-rbdy(jb-1))*(zmid-zbdy(jb-1))/(zbdy(jb)-zbdy(jb-1))
              inumt = inumt + 1
              if (inumt.gt.2) exit
              jb = jb + 1
           endif
           jb = jb + 1
        enddo

        if (inumt.lt.3) then
           write(nonlin,*)'%tdb_get_rzbdy: boundary lacks two midplane intersections.'
           ier = 4
           return
        endif

        if (rmid(1).lt.rmid(2)) then
           r1 = rmid(1);  r2 = rmid(2)
        else
           r1 = rmid(2);  r2 = rmid(1)
        endif
     endif
  endif
end subroutine tdb_get_rzbdy

!=======================================================================
subroutine tdb_get_bdysize(d,inbdy)

  !  Return size of boundary coord arrays to go with
  !  rbdy(theta),zbdy(theta) data vs. time

  use trdatbuf_obj
  implicit NONE

  type(trdatbuf), intent(in) :: d
  integer, intent(out)       :: inbdy  ! the array size, or 0 if none found.

  inbdy = d%nbdy

end subroutine tdb_get_bdysize
