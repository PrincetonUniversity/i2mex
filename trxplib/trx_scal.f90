!  trx_scal -- get conversion factor, read & scale scalar
!  trx_kscal -- read & scale scalar
!  trxc_scal -- C interface to trx_scal
!  trxc_kscal -- C interface to trxc_scal
!
!--------------------------------------------------------------------
!  get units conversion, then, read & scale scalar (trx_kscal,
!  below)
!
subroutine trx_scal(zname,zuns_mks,zval,ierr)
 
  use trx_module
  implicit NONE
 
! ...arguments
 
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zname          ! name of desired scalar
  character(*), intent(out) :: zuns_mks      ! mks units label of profile
  real*8, intent(out)        :: zval         ! scalar value returned
  integer, intent(out)     :: ierr           ! completion code, 0=OK
 
! ...local
 
  real*8                   :: zconv          ! units conversion factor
  integer                  :: iwarn          ! units conversion warning
 
! ------------------------------
 
  call trx_mks_conv(zname,zconv,zuns_mks,iwarn)
  call trx_kscal(zname,zconv,zval,ierr)
 
  return
end subroutine trx_scal
!--------------------------------------------------------------------
!  read a TRANSP scalar & interpolate/extrapolate to flux surfaces
!  from zone centers if necessary; scale by units conversion factor
!
subroutine trx_kscal(zname,zconv,zval,ierr)
 
  use trx_module
  implicit NONE
 
! ...arguments
 
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zname          ! name of desired scalar
  real*8, intent(in)       :: zconv          ! units conversion factor
  real*8, intent(out)      :: zval           ! scalar value returned
  integer, intent(out)     :: ierr           ! completion code, 0=OK
 
! ...local
!
  character(64) :: zlabel            ! label
  character(32) :: zunits            ! units
!
  integer :: istype                  ! type code
  integer :: lunzer
!
  real zvalr4
! ------------------------------
!
  call trx_ready('trx_scal',ierr)
  if(ierr.ne.0) return
!
  call t1scalar(zname,zlabel,zunits,time0,delta_t, zvalr4, ierr)
  zval=zvalr4*zconv
!
  return
end subroutine trx_kscal
