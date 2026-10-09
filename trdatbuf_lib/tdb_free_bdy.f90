!
! Set of trdatbuf subroutines for accessing free boundary data.
!

!
! --------------- tdb_freebdy_timefac ------------
! get the interpolation factor at ztime.  All free boundary data is
! currently mapped to LTIME2.
!
subroutine tdb_freebdy_timefac(d,ztime,it,zf)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  real*8,  intent(in)  :: ztime   ! time (seconds)
  integer, intent(out) :: it      ! time bin
  real*8,  intent(out) :: zf      ! interpolation factor w/in bin

  integer :: ilt, int

  ilt = d%ltime2
  int = d%ntime2
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,ztime,it,zf)
end subroutine tdb_freebdy_timefac

!
! --------------- tdb_freebdy_pretimefac ------------
! get the interpolation factor at ztime for the pre TINIT data.  All pre TINIT data is
! currently mapped to LTIMEC.
!
subroutine tdb_freebdy_pretimefac(d,ztime,it,zf)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  real*8,  intent(in)  :: ztime   ! time (seconds)
  integer, intent(out) :: it      ! time bin
  real*8,  intent(out) :: zf      ! interpolation factor w/in bin

  integer :: ilt, int

  ilt = d%ltimec
  int = d%ntimec

  if (ilt<=0 .or. int<=0) then
     print *, '?tdb_freebdy_pretimefac: unexpectedly no TINIT_PRECOIL data'
     call bad_exit
  end if

  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,ztime,it,zf)
end subroutine tdb_freebdy_pretimefac

!
! ------------ tdb_freebdy_dims ----------
! Return the dimensions associated with the free boundary data.
!
subroutine tdb_freebdy_dims(d, npfc)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  integer, intent(out) :: npfc    ! number of poloidal field coils or 0 if none

  npfc = d%NXPFC
end subroutine tdb_freebdy_dims


!
! ------------ tdb_freebdy_pfcnames -------------
! Return the poloidal field coil names.  An error is
! returned if there are no coil names or the dimension of the input
! array is incorrect.
!
subroutine tdb_freebdy_pfcnames(d, n, pfnam, ier)
  use trdatbuf_obj
  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  integer, intent(in)  :: n       ! number of coils expected
  character(len=8), dimension(n), intent(out) :: pfnam  ! pf coil names returned
  integer, intent(out) :: ier     ! nonzero on error

  integer :: ldpfn   ! datbuf address of start of pf name data
  integer :: npfc    ! number of pf coils from trdatbuf
  integer :: ix      ! loop variable

  integer,external :: lunmsg_tdb

  ier = 0

  npfc  = d%NXPFC
  ldpfn = d%LDPFN

  if (ldpfn<=0 .or. npfc<=0) then
     write(lunmsg_tdb(0),'(a)') '?tdb_freebdy_pfcnames: there is no poloidal field coilname data'
     ier=1 ; return
  end if

  if (n/=npfc) then
     write(lunmsg_tdb(0),'(a)')'?tdb_freebdy_pfcnames: mismatch in number of PF coilnames'
     write(lunmsg_tdb(0),'(a,i5)')'      expected:  ', n
     write(lunmsg_tdb(0),'(a,i5)')'      available: ', npfc
     ier=1 ; return
  end if

  do ix=1,npfc
     call str2real(pfnam(ix), d%datbuf(ldpfn+8*(ix-1)), -1)
  enddo !ix
end subroutine tdb_freebdy_pfcnames


!
! ------------ tdb_freebdy_pfccurs -------------
! Return the poloidal field coil currents at a time.  An error is
! returned if there are no currents or the dimension of the input
! array is incorrect.
!
subroutine tdb_freebdy_pfccurs(d, ztime, n, cur, ier)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  real*8,  intent(in)  :: ztime   ! time (seconds)
  integer, intent(in)  :: n       ! number of currents expected
  real*8,  intent(out) :: cur(n)  ! poloidal field coil currents returned
  integer, intent(out) :: ier     ! nonzero on error

  integer :: ldpfc   ! datbuf address of start of pfc data
  integer :: npfc    ! number of currents from trdatbuf
  integer :: nt2     ! number of time points for data
  integer :: it      ! time bin
  integer :: ix      ! loop variable
  real*8  :: zf      ! interpolation factor w/in bin

  integer,external :: lunmsg_tdb

  ier=0

  npfc  = d%NXPFC
  ldpfc = d%LDPFC

  if (ldpfc<=0 .or. npfc<=0) then
     write(lunmsg_tdb(0),'(a)') '?tdb_freebdy_pfccurs: there is no poloidal field coil current data'
     ier=1 ; return
  end if

  if (n/=npfc) then
     write(lunmsg_tdb(0),'(a)')    '?tdb_freebdy_pfccurs: mismatch in number of PFC currents'
     write(lunmsg_tdb(0),'(a,i5)') '      expected:  ', n
     write(lunmsg_tdb(0),'(a,i5)') '      available: ', npfc
     ier=1 ; return
  end if

  call tdb_freebdy_timefac(d,ztime,it,zf)

  nt2 = d%NTIME2
  do ix = 1, npfc
     cur(ix) = (1.d0-zf)*d%datbuf(ldpfc+(it-1)+(ix-1)*nt2)+zf*d%datbuf(ldpfc+it+(ix-1)*nt2)
  end do
end subroutine tdb_freebdy_pfccurs

!
! ------------ tdb_freebdy_pre_pfc_pcur -------------
! Return the poloidal field coil currents and plamsa current at a pre TINIT time.  An error is
! returned if there are no currents or the dimension of the input array is incorrect.
!
subroutine tdb_freebdy_pre_pfc_pcur(d, ztime, n, cur, pcur, ier)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  implicit NONE

  type(trdatbuf)       :: d       ! trdatbuf
  real*8,  intent(in)  :: ztime   ! time (seconds)
  integer, intent(in)  :: n       ! number of currents expected
  real*8,  intent(out) :: cur(n)  ! poloidal field coil currents returned
  real*8,  intent(out) :: pcur    ! plasma current returned, if not defined in input this will be 0.
  integer, intent(out) :: ier     ! nonzero on error

  integer :: ldpfc   ! datbuf address of start of pfc data
  integer :: npfc    ! number of currents from trdatbuf
  integer :: nt2     ! number of time points for data
  integer :: it      ! time bin
  integer :: ix      ! loop variable
  real*8  :: zf      ! interpolation factor w/in bin

  integer,external :: lunmsg_tdb

  ier=0

  npfc  = d%NXPFC
  ldpfc = d%LDPFC_PRE

  if (ldpfc<=0 .or. npfc<=0 .or. d%NTIMEC<=0) then
     write(lunmsg_tdb(0),'(a)') '?tdb_freebdy_pre_pfc_pcur: there is no poloidal field coil current data pre TINIT'
     ier=1 ; return
  end if

  if (n/=npfc) then
     write(lunmsg_tdb(0),'(a)')    '?tdb_freebdy_pre_pfc_pcur: mismatch in number of PFC currents'
     write(lunmsg_tdb(0),'(a,i5)') '      expected:  ', n
     write(lunmsg_tdb(0),'(a,i5)') '      available: ', npfc
     ier=1 ; return
  end if

  call tdb_freebdy_pretimefac(d,ztime,it,zf)

  nt2 = d%NTIMEC
  do ix = 1, npfc
     cur(ix) = (1.d0-zf)*d%datbuf(ldpfc+(it-1)+(ix-1)*nt2)+zf*d%datbuf(ldpfc+it+(ix-1)*nt2)
  end do

  pcur = (1.d0-zf)*d%datbuf(ldpfc+(it-1)+npfc*nt2)+zf*d%datbuf(ldpfc+it+npfc*nt2)
end subroutine tdb_freebdy_pre_pfc_pcur
