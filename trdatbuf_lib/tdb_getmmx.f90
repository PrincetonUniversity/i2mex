subroutine tdb_getmmx(d,zt,xib,zrmc2,zymc2,mj,mimom,lcentr,nzp1,zdrshaf)
 
  use trdatbuf_obj
  use tdbsub_uts
  IMPLICIT NONE
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
  type (trdatbuf) :: d
  REAL*8,intent(in) :: zt    ! time at which to fetch moments
 
  integer, intent(in) :: mj  ! flux surface array dimension for zrmc2,zymc2
  integer, intent(in) :: mimom ! moments array dimension for zrmc2,zymc2

  real*8,intent(in) :: xib(mj) ! flux surfaces (xib(j)=sqrt(Phi/Philim) at
  ! surface "j").
  
  REAL*8,intent(out) :: zrmc2(mj,0:mimom,2),zymc2(mj,0:mimom,2)
                           ! asymmetric Fourier moments set at xi bdys
 
  REAL*8,intent(out) :: zdrshaf(mj)
                           ! "Shafranov" shift of interior surfaces
                           ! relative to the boundary surface
 
  integer :: lcentr        ! index to magnetic axis
  integer :: nzp1          ! no. of surfaces including mag. axis

  ! using the MMX data, fetch the current equilibrium geometry.
  ! this routine MUST NOT be called unless the MMX data exists.
 
  !----------------------------------------------------
  ! local:
 
  ! local splines...
 
  REAL*8, dimension(:), allocatable :: zwk
  REAL*8, dimension(:,:), allocatable :: zxpkg
  REAL*8, dimension(:,:,:,:), allocatable :: zsrmc,zsymc
 
  ! for time interpolation:
 
  integer it    ! time index
  REAL*8 zf       ! fraction to next index pt.
 
  integer iadr,iady,ics,im,ix,j,ilt,int,idum,ierr
  REAL*8 zdum,zrin,zrout,zymid
 
  integer iadci   ! indexing function trulib/iadci.for
  integer ict(3)

  integer :: lunmsg_tdb
 
  integer :: lep1,lcp1,nzones,imaxe

  data ict/1,0,0/
 
  !----------------------------------------------------

  lcp1=lcentr+1
  nzones=nzp1-1
  lep1=lcentr+nzones
 
  if(d%LDMMX.eq.0) then
     write(lunmsg_tdb(0),*) ' ??? GETMMX called but there is no MMX data.'
     call bad_exit
  endif
 
  if(d%NXMMX.le.1) then
     !  if d%NXMMX=1 there is only boundary data and this routine should
     !  still not be called.
     write(lunmsg_tdb(0),*) ' ??? GETMMX called but there is no MMX x grid.'
     call bad_exit
  endif
 
  if(d%NMOMD.ne.d%MMAX) then
     write(lunmsg_tdb(0),*) ' ??? GETMMX: mmax=',d%MMAX,' nmomd=',d%NMOMD
     write(lunmsg_tdb(0),*) '     these quantities were expected to be equal.'
     call bad_exit
  endif
 
  imaxe = d%mmx_maxe  ! max scaling exponent for moments
  if(imaxe.eq.0) imaxe = 16  ! the old default set here.

  zrmc2=0; zymc2=0
 
  ilt=d%ltime2
  int=d%ntime2
  call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,zt,it,zf)
 
  allocate(zwk(d%NXMMX))
  allocate(zxpkg(d%NXMMX,4),zsrmc(4,d%NXMMX,0:d%MMAX,2),zsymc(4,d%NXMMX,0:d%MMAX,2))
 
  zdum=0.0_R8
  call genxpkg(d%NXMMX,d%DATBUF(d%LXMMX),zxpkg,0,0,0,zdum,-3,ierr)
  if(ierr.ne.0) then
     call errmsg_exit(' ?getmmx: unexpected spline "genxpkg" error.')
  endif
 
  !  time interpolate scaled moments on orig. data grid
 
  do ics=1,2
     do im=0,d%MMAX
        do ix=1,d%NXMMX
           iadr=iadci(d, d%LDMMX, it, ix, im, ics)
           zsrmc(1,ix,im,ics)=d%DATBUF(iadr)*(1-zf)+d%DATBUF(iadr+1)*zf
           iady=iadci(d, d%LDMMX, it, ix, im, ics+2)
           zsymc(1,ix,im,ics)=d%DATBUF(iady)*(1-zf)+d%DATBUF(iady+1)*zf
        enddo
     enddo
  enddo
 
  !  set up splines
 
  zdum=0.0_R8
  do ics=1,2
     do im=0,d%MMAX
        call cspline(d%DATBUF(d%LXMMX),d%NXMMX,zsrmc(1,1,im,ics),0,zdum,0,zdum,zwk, &
             d%NXMMX,idum,ierr)
        if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "cspline" error.')
        call cspline(d%DATBUF(d%LXMMX),d%NXMMX,zsymc(1,1,im,ics),0,zdum,0,zdum,zwk, &
             d%NXMMX,idum,ierr)
        if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "cspline" error.')
     enddo
  enddo
 
  !  map to TRANSP grid
 
  do ics=1,2
     do im=0,d%MMAX
        if(im.eq.0) then
 
           ! 0'th moment can have non zero value on axis
 
           call spvec(ict,nzp1,xib(lcentr),nzp1, &
                zrmc2(lcentr:lep1,im,ics),d%NXMMX,zxpkg,zsrmc(1,1,im,ics), &
                idum,ierr)
           if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "spvec" (spline evaluation) error 1.')
           call spvec(ict,nzp1,xib(lcentr),nzp1, &
                zymc2(lcentr:lep1,im,ics),d%NXMMX,zxpkg,zsymc(1,1,im,ics), &
                idum,ierr)
           if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "spvec" (spline evaluation) error 2.')
 
        else
 
           ! higher moments always zero on axis
 
           call spvec(ict,nzones,xib(lcp1),nzones, &
                zrmc2(lcp1:lep1,im,ics),d%NXMMX,zxpkg,zsrmc(1,1,im,ics), &
                idum,ierr)
           if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "spvec" (spline evaluation) error 1.')
           call spvec(ict,nzones,xib(lcp1),nzones, &
                zymc2(lcp1:lep1,im,ics),d%NXMMX,zxpkg,zsymc(1,1,im,ics), &
                idum,ierr)
           if(ierr.ne.0) call errmsg_exit(' ?getmmx: unexpected "spvec" (spline evaluation) error 2.')
 
           ! scale the moments
 
           zrmc2(lcp1:lep1,im,ics)=zrmc2(lcp1:lep1,im,ics)* &
                xib(lcp1:lep1)**min(im,imaxe)
           zymc2(lcp1:lep1,im,ics)=zymc2(lcp1:lep1,im,ics)* &
                xib(lcp1:lep1)**min(im,imaxe)
 
        endif
     enddo
  enddo
 
  !  compute axis shift
 
  zdrshaf = 0
  zdrshaf(lcentr)=zrmc2(lcentr,0,1)
  do j=lcp1,lep1
     call r8_getmpa(zrmc2(j,0:mimom,1:2),zymc2(j,0:mimom,1:2),mimom,d%NMOMD, &
          zrin,zrout,zymid)
     zdrshaf(j)=0.5_R8*(zrin+zrout)
  enddo
 
  !  shift is defined relative to outer bdy
 
  zdrshaf(lcentr:lep1) = zdrshaf(lcentr:lep1) - zdrshaf(lep1)
 
end subroutine tdb_getmmx

