subroutine trx_mks_conv(zname,zconv,zuns_out,iwarn)
!
!  find MKS conversion factor for named quantity
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  character(*), intent(in) :: zname    ! name of data item
  real*8, intent(out) :: zconv         ! conversion factor for data item
  character(*), intent(out) :: zuns_out ! units label (after mks conversion)
  integer, intent(out) :: iwarn        ! conversion warning flag, =0: OK
!
!  if the conversion cannot be done, zconv=1.0 is returned and a warning
!  message is written, and a warning flag is set.
!
!-------------
!  ...local
!
  character(64) :: zlbl
  character(32) :: zuns,zuns_mks
  integer :: imulti
  integer :: istype
!
  real zconvr4
!-------------
!
  call rplabel(zname,zlbl,zuns,imulti,istype)
  call rp_mks_conv(zname,zuns,zuns_mks,zconvr4,iwarn)
!
  zconv=zconvr4
  zuns_out=zuns_mks
!
  return
  end
