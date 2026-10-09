logical function tdb_fallback(d,ztri,ztime,zresult)
  !
  ! ** NEW SEP 1984 **  bolometer profile data can "fall back to time
  !     series" in which case interpretation is that the total 
  !     radiated power (Watts) was measured, instead of a power 
  !     density

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: ztri  ! trigraph of candidate data
  ! for "fall back to time series"
  real*8, intent(in) :: ztime        ! time to which to interpolate
  real*8, intent(out) :: zresult     ! result, if data is found
  ! the function value is .TRUE., and the interpolated value is returned,
  ! if "fallback to time series" data is present; otherwise the function
  ! returns .FALSE. and zresult=ZERO.

  !------------------------------------------

  character*3 ztest
  integer :: if,ib0,ib1,ilt,int
  real*8 :: zf

  !------------------------------------------
  !  this could be generated code, if we ever
  !  need a lot of these... (dmc May 2005)

  ztest=ztri(1:min(3,len(ztri)))
  if(ztest.eq.'BOL') then
     if(d%nxbol.eq.1) then
        ilt=d%ltime2
        int=d%ntime2
        tdb_fallback=.TRUE.
        call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,ztime,if,zf)
        IB1=d%LFBOL+IF
        IB0=IB1-1
        zresult=d%DATBUF(IB0)+ZF*(d%DATBUF(IB1)-d%DATBUF(IB0))
     else
        tdb_fallback=.FALSE.
        zresult=ZERO
     endif
  else
     tdb_fallback=.FALSE.
     zresult=ZERO
  endif

end function tdb_fallback
