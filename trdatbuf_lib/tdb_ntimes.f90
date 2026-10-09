subroutine tdb_ntimes(d,intime1,intime2)
  !
  !  return sizes of trdatbuf input 1d f(t) timebase and 2d f(t,x) timebase

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(out) :: intime1
  integer, intent(out) :: intime2

  intime1 = d%ntime1
  intime2 = d%ntime2

end subroutine tdb_ntimes

subroutine tdb_range_ntimes1(d,zt1,zt2,intime1)

  !  return no. of pts in tdb input 1d f(t) timebase within the time
  !  range (zt1,zt2)

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(in) :: zt1,zt2  ! time range
  integer, intent(out) :: intime1 ! no. of time points returned

  !----------------
  integer :: it1,it2,ilt,int
  real*8 :: zf1,zf2
  !----------------

  ilt=d%ltime1
  int=d%ntime1
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt1,it1,zf1)
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt2,it2,zf2)

  if(d%datbuf(ilt+it1-1).lt.zt1) it1=it1+1
  if(d%datbuf(ilt+it2).le.zt2) it2=it2+1

  intime1=it2-it1+1

end subroutine tdb_range_ntimes1

subroutine tdb_range_times1(d,zt1,zt2,ztimes,ierr)

  !  return the time points in tdb input 1d f(t) timebase within the time
  !  range (zt1,zt2)

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(in) :: zt1,zt2 ! time range
  real*8, intent(inout) :: ztimes(:) ! time points returned
  integer, intent(out) :: ierr ! error (ierr=1) if ztimes array too small

  !----------------
  integer :: it1,it2,ilt,int,intime1
  real*8 :: zf1,zf2
  !----------------

  ilt=d%ltime1
  int=d%ntime1
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt1,it1,zf1)
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt2,it2,zf2)

  if(d%datbuf(ilt+it1-1).lt.zt1) it1=it1+1
  if(d%datbuf(ilt+it2).le.zt2) it2=it2+1

  intime1=it2-it1+1

  ztimes=0
  if(size(ztimes).lt.intime1) then
     ierr=1
  else
     ierr=0
     ztimes(1:intime1)=d%datbuf(ilt+it1-1:ilt+it2-1)
  endif

end subroutine tdb_range_times1

subroutine tdb_range_ntimes2(d,zt1,zt2,intime2)

  !  return no. of pts in tdb input 2d f(x,t) timebase within the time
  !  range (zt1,zt2)

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(in) :: zt1,zt2  ! time range
  integer, intent(out) :: intime2 ! no. of time points returned

  !----------------
  integer :: it1,it2,ilt,int
  real*8 :: zf1,zf2
  !----------------

  ilt=d%ltime2
  int=d%ntime2
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt1,it1,zf1)
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt2,it2,zf2)

  if(d%datbuf(ilt+it1-1).lt.zt1) it1=it1+1
  if(d%datbuf(ilt+it2).le.zt2) it2=it2+1

  intime2=it2-it1+1

end subroutine tdb_range_ntimes2

subroutine tdb_range_times2(d,zt1,zt2,ztimes,ierr)

  !  return the time points in tdb input 2d f(x,t) timebase within the time
  !  range (zt1,zt2)

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(in) :: zt1,zt2 ! time range
  real*8, intent(inout) :: ztimes(:) ! time points returned
  integer, intent(out) :: ierr ! error (ierr=1) if ztimes array too small

  !----------------
  integer :: it1,it2,ilt,int,intime2
  real*8 :: zf1,zf2
  !----------------

  ilt=d%ltime2
  int=d%ntime2
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt1,it1,zf1)
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt2,it2,zf2)

  if(d%datbuf(ilt+it1-1).lt.zt1) it1=it1+1
  if(d%datbuf(ilt+it2).le.zt2) it2=it2+1

  intime2=it2-it1+1

  ztimes=0
  if(size(ztimes).lt.intime2) then
     ierr=1
  else
     ierr=0
     ztimes(1:intime2)=d%datbuf(ilt+it1-1:ilt+it2-1)
  endif

end subroutine tdb_range_times2
