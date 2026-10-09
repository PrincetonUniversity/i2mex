subroutine tdb_nzones(nzones,nrmaj,nrsym)
  !
  !  from basic grid size produce midplane grid sizes used by the trdatbuf
  !  profile interpolation routines

  integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)

  integer, intent(out) :: nrmaj  ! (2*nzones+1) no. of flux surface
  !                                midplane intercepts
  integer, intent(out) :: nrsym  ! (4*nzones+5) size of midplane test
  !                                grid (zone bdys & surfaces + 2
  !                                extrapolation points at each end).

  nrmaj = 2*nzones + 1
  nrsym = 4*nzones + 5

end subroutine tdb_nzones

subroutine tdb_xilmp(nzones,xibdys,xilmp)

  integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)
  real*8, dimension(:), intent(in) :: xibdys  ! (nzones+1) x @ zone bdys

  real*8, dimension(:), intent(out) :: xilmp  ! (2*nzones+1) x @midplane
  !                                flux surface intercepts

  integer :: inzp1,i2nzp1,ierr,i,inc
  integer :: lunmsg_tdb

  !--------------------------------

  inzp1 = nzones+1
  i2nzp1 = inzp1 + nzones

  ierr=0
  if(size(xibdys).ne.inzp1) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xilmp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(xibdys) = nzones + 1, but found:',&
          'size(xibdys) = ',size(xibdys)
     ierr = ierr+1
  endif

  if(size(xilmp).ne.i2nzp1) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xilmp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(xilmp) = 2*nzones + 1, but found:',&
          ' size(xilmp) = ',size(xilmp)
     ierr = ierr+1
  endif
 
  if(ierr.gt.0) then
     xilmp = 0
     return
  endif

  xilmp(inzp1)=xibdys(1)

  inc=0
  do i=1,nzones
     inc=inc+1
     xilmp(inzp1+inc)=xibdys(i+1)
     xilmp(inzp1-inc)=-xibdys(i+1)
  enddo

end subroutine tdb_xilmp


subroutine tdb_xisymp(nzones,xibdys,rmajmp,xirsym,rmjsym)
  
  !  generate test grids for looking at symmetrization of profile data

  use trdatbuf_iface, only: tdb_xilmp
  implicit NONE

  integer, intent(in) :: nzones  ! target grid, no. of flux zones (in)
  real*8, dimension(:), intent(in) :: xibdys  ! (nzones+1) x @ zone bdys
  real*8, dimension(:), intent(in) :: rmajmp  ! (2*nzones+1) major radius
  !         at flux surface midplane intercepts (cm)
  real*8, dimension(:), intent(out) :: xirsym ! (4*nzones+5) x @midplane
  !         (double-density grid with 2 pt extension @ bdys)
  real*8, dimension(:), intent(out) :: rmjsym ! (4*nzones+5) R @midplane
  !         (double-density grid with 2 pt extension @ bdys) (cm)

  integer :: inzp1,i2nzp1,i4nzp5,ierr
  integer :: icen,icenr,inc,incr,j,incrp,incri,ir,i4,i3,i2,i1
  integer :: lunmsg_tdb

  real*8, dimension(:), allocatable :: xilmp
  real*8 :: zdr,zdx

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !--------------------------------

  inzp1 = nzones+1
  i2nzp1 = 2*nzones+1
  i4nzp5 = 4*nzones+5

  ierr=0
  if(size(xibdys).ne.inzp1) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xisymp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(xibdys) = nzones + 1, but found:',&
          'size(xibdys) = ',size(xibdys)
     ierr = ierr+1
  endif

  if(size(rmajmp).ne.i2nzp1) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xisymp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(rmajmp) = 2*nzones + 1, but found:',&
          ' size(rmajmp) = ',size(rmajmp)
     ierr = ierr+1
  endif

  if(size(xirsym).ne.i4nzp5) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xisymp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(xirsym) = 4*nzones + 5, but found:',&
          ' size(xirsym) = ',size(xirsym)
     ierr = ierr+1
  endif

  if(size(rmjsym).ne.i4nzp5) then
     write(lunmsg_tdb(0),*) ' ?trdatbuf_lib/tdb_xisymp: with nzones = ',nzones
     write(lunmsg_tdb(0),*) '  expected size(rmjsym) = 4*nzones + 5, but found:',&
          ' size(rmjsym) = ',size(rmjsym)
     ierr = ierr+1
  endif

  if(ierr.ne.0) then
     xirsym=0
     rmjsym=0
     return
  endif

  allocate(xilmp(i2nzp1))
  call tdb_xilmp(nzones,xibdys,xilmp)

  ICEN=NZONES+1  ! CENTER OF RMAJMP GRID (MAGNETIC AXIS)

  ICENR=2*NZONES+3  ! CENTER OF RMJSYM GRID

  !  at axis

  RMJSYM(ICENR)=RMAJMP(ICEN)
  XIRSYM(ICENR)=XILMP(ICEN)

  !  WITHIN THE PLASMA:
 
  INC=0
  INCR=0
  DO J=1,NZONES
     INC=INC+1
     INCRP=INCR
     INCR=INCR+2

     !  THE FLUX SURFACES:

     RMJSYM(ICENR+INCR)=RMAJMP(ICEN+INC)
     XIRSYM(ICENR+INCR)=XILMP(ICEN+INC)

     RMJSYM(ICENR-INCR)=RMAJMP(ICEN-INC)
     XIRSYM(ICENR-INCR)=XILMP(ICEN-INC)

     !  BTW THE FLUX SURFACES:

     INCRI=INCR-1
     RMJSYM(ICENR+INCRI)= &
          0.5_R8*(RMJSYM(ICENR+INCRP)+RMJSYM(ICENR+INCR))
     XIRSYM(ICENR+INCRI)= &
          0.5_R8*(XIRSYM(ICENR+INCRP)+XIRSYM(ICENR+INCR))
     RMJSYM(ICENR-INCRI)= &
          0.5_R8*(RMJSYM(ICENR-INCRP)+RMJSYM(ICENR-INCR))
     XIRSYM(ICENR-INCRI)= &
          0.5_R8*(XIRSYM(ICENR-INCRP)+XIRSYM(ICENR-INCR))

  enddo

  !  BEYOND THE PLASMA BDY

  !  INSIDE
  ZDR=RMJSYM(4)-RMJSYM(3)
  ZDX=XIRSYM(4)-XIRSYM(3)
  RMJSYM(2)=RMJSYM(3)-ZDR
  RMJSYM(1)=RMJSYM(2)-ZDR
  XIRSYM(2)=XIRSYM(3)-ZDX
  XIRSYM(1)=XIRSYM(2)-ZDX

  !  OUTSIDE
  I4=4*NZONES+2
  I3=4*NZONES+3
  I2=4*NZONES+4
  I1=4*NZONES+5
  ZDR=RMJSYM(I4)-RMJSYM(I3)
  ZDX=XIRSYM(I4)-XIRSYM(I3)
  RMJSYM(I2)=RMJSYM(I3)-ZDR
  RMJSYM(I1)=RMJSYM(I2)-ZDR
  XIRSYM(I2)=XIRSYM(I3)-ZDX
  XIRSYM(I1)=XIRSYM(I2)-ZDX

  deallocate(xilmp)

end subroutine tdb_xisymp
