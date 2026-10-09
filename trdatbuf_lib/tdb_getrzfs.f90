subroutine tdb_getrzfs_sizes(d,mth,mxi)
 
  !  get th and xi grid sizes for explicit R(theta,xi), Z(theta,xi)
  !  flux surface data

  use trdatbuf_obj
  IMPLICIT NONE

  type (trdatbuf) :: d
  integer, intent(out) :: mth  ! th (poloidal angle) grid size
  integer, intent(out) :: mxi  ! xi (sqrt(phi/philim)) radial grid size
  !-------------------------

  mth = d%nthfs
  mxi = d%nxfs

end subroutine tdb_getrzfs_sizes

subroutine tdb_getrzfs_xi(d,xi,nmax,ngot)

  !  get xi grid for explicit R(theta,xi), Z(theta,xi)
  !  flux surface data

  use trdatbuf_obj
  IMPLICIT NONE

  type (trdatbuf) :: d
  integer, intent(in) :: nmax      ! size of xi array provided
  real*8, intent(out) :: xi(nmax)  ! xi array
  integer, intent(out) :: ngot     ! number of xi values returned

  ! ngot=0 & xi = 0 is returned if there is no such data available or
  ! if nmax is too small.  In either case a warning message is printed.

  !-------------------------
  integer :: lunmsg_tdb
  !-------------------------

  ngot = 0
  xi = 0

  if(d%nxfs.eq.0) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs_xi:  no data available'
     return
  endif

  if(d%nxfs.gt.nmax) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs_xi:  passed xi grid size too small:'
     write(lunmsg_tdb(0),*) '  nmax (passed): ',nmax
     write(lunmsg_tdb(0),*) '  grid size in "d": ',d%nxfs
     return
  endif

  ngot = d%nxfs

  xi(1:ngot) = d%datbuf(d%lxfs:d%lxfs+d%nxfs-1)

end subroutine tdb_getrzfs_xi

subroutine tdb_getrzfs_th(d,th,nmax,ngot)

  !  get th grid for explicit R(theta,xi), Z(theta,xi)
  !  flux surface data

  use trdatbuf_obj
  IMPLICIT NONE

  type (trdatbuf) :: d
  integer, intent(in) :: nmax      ! size of th array provided
  real*8, intent(out) :: th(nmax)  ! th array
  integer, intent(out) :: ngot     ! number of th values returned

  ! ngot=0 & th = 0 is returned if there is no such data available or
  ! if nmax is too small.  In either case a warning message is printed.

  !-------------------------
  integer :: lunmsg_tdb
  !-------------------------

  ngot = 0
  th = 0

  if(d%nthfs.eq.0) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs_th:  no data available'
     return
  endif

  if(d%nthfs.gt.nmax) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs_th:  passed th grid size too small:'
     write(lunmsg_tdb(0),*) '  nmax (passed): ',nmax
     write(lunmsg_tdb(0),*) '  grid size in "d": ',d%nthfs
     return
  endif

  ngot = d%nthfs

  th(1:ngot) = d%datbuf(d%lthfs:d%lthfs+d%nthfs-1)

end subroutine tdb_getrzfs_th

subroutine tdb_getrzfs(d,zt,mth,mxi,rsurf,zsurf,ierr)

  !  get explicit R(theta,xi), Z(theta,xi) flux surface data interpolated
  !  to time "zt".  Array sizes provided must match the available data
  !  precisely, else ierr=1 is set and a message is printed.

  use trdatbuf_obj
  use tdbsub_uts
  IMPLICIT NONE

  type (trdatbuf) :: d
  real*8, intent(in) :: zt   ! time (seconds) to which to interpolate data

  integer, intent(in) :: mth,mxi  ! array dimensions -- must match!
  !  these are for the th (poloidal angle) and xi (radial flux coordinate)
  !  grids respectively

  real*8, intent(out) :: rsurf(mth,mxi)  ! R(th,xi)
  real*8, intent(out) :: zsurf(mth,mxi)  ! Z(th,xi)

  integer, intent(out) :: ierr           ! status code, 0=OK

  !--------------------------
  integer :: lunmsg_tdb
  integer :: ixi,ith,iad,ilt,int,it,ildr,ildz,incr
  real*8 :: zf
  !--------------------------

  ierr=0

  if(mth.ne.d%nthfs) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs: "th" array dimension mismatch:'
     write(lunmsg_tdb(0),*) '  passed value (mth):  ',mth
     write(lunmsg_tdb(0),*) '  stored value in "d": ',d%nthfs
     ierr=1
  endif

  if(mxi.ne.d%nxfs) then
     write(lunmsg_tdb(0),*) ' ?tdb_getrzfs: "xi" array dimension mismatch:'
     write(lunmsg_tdb(0),*) '  passed value (mxi):  ',mxi
     write(lunmsg_tdb(0),*) '  stored value in "d": ',d%nxfs
     ierr=1
  endif

  if(ierr.ne.0) then
     rsurf = 0
     zsurf = 0
     return
  endif

  !--------------------------------------
  !  error checks complete; interpolate the data

  !  get time interpolation factors (it:it+1) (zf=0.0 ->it; zf=1.0 ->it+1)
   
  ilt=d%ltime2
  int=d%ntime2
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt,it,zf)

  ildr=d%lrfs
  ildz=d%lzfs

  do ith=1,mth
     do ixi=1,mxi
        incr = int*((ixi-1)*mth + (ith-1))

        iad = ildr + incr + (it-1)
        rsurf(ith,ixi) = d%datbuf(iad)*(ONE-zf) + d%datbuf(iad+1)*zf

        iad = ildz + incr + (it-1)
        zsurf(ith,ixi) = d%datbuf(iad)*(ONE-zf) + d%datbuf(iad+1)*zf

     enddo
  enddo

end subroutine tdb_getrzfs
