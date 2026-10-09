subroutine trx_wall(ilun,inumR,inumZ,ierr)
!  get the limiter locations; form (R,Z) grid to cover just this space...
!  see trx_wall_RZ (below).
!
  use trx_module, only: symZ

  implicit NONE
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
  integer ilun                ! lun for file read (input)
  integer inumR,inumZ         ! cartesian overlay grid sizes (input)
  integer ierr                ! completion code:  0=OK
!
  real*8 Rmin,Rmax            ! Rmin,Rmax of grid (0.0 for automatic grid)
  real*8 Zmin,Zmax            ! Zmin,Zmax of grid (0.0 for automatic grid)
!
  real*8 :: ztol
!
  symZ = .true.

  Rmin=0.0_R8
  Rmax=0.0_R8
  Zmin=0.0_R8
  Zmax=0.0_R8
  call trx_wall_rz(ilun,Rmin,Rmax,inumR,Zmin,Zmax,inumZ,ierr)
  return
end subroutine trx_wall

subroutine trx_wall_freebdy(ilun,ierr)

  use trx_module
  implicit NONE

  !  use the free bdy (R,Z) grid -- consider it an error if there is none

  integer, intent(in) :: ilun
  integer, intent(out) :: ierr

  !------------------------------
  integer :: lunzer
  !------------------------------

  if(min(nRfree,nZfree).eq.0) then
     ierr=1
     write(lunzer(0),*) ' ? trx_wall_freebdy: called w/o free bdy data.'
     return
  endif

  call trx_wall_grid(ilun,Rgrid_free,nRfree,Zgrid_free,nZfree,ierr)

end subroutine trx_wall_freebdy

subroutine trx_wall_grid(ilun,r8grid,inumR,z8grid,inumZ,ierr)
!  get the limiter locations; apply a prescribed (R,Z) grid


!  get the limiter locations; form (R,Z) grid to cover just this space...
!  see trx_wall_RZ (below).
!
  use trx_module, only: symZ

  implicit NONE
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
  integer ilun                ! lun for file read (input)
  integer inumR,inumZ         ! cartesian overlay grid sizes (input)
  real*8 R8grid(inumR),Z8grid(inumZ)  ! the grids themselves (input)
  integer ierr                ! completion code:  0=OK
!
  real*8 Rmin,Rmax            ! Rmin,Rmax of grid (0.0 for automatic grid)
  real*8 Zmin,Zmax            ! Zmin,Zmax of grid (0.0 for automatic grid)
!
  real*8 :: ztol
  integer :: id_R,id_Z
  integer :: lunzer
!
  symZ = .false.

  Rmin=-1.0_R8
  Rmax=-1.0_R8
  Zmin=-1.0_R8
  Zmax=-1.0_R8
  call trx_wall_rz(ilun,Rmin,Rmax,inumR,Zmin,Zmax,inumZ,ierr)

  ztol=0.01d0*max((Rmax-Rmin),(Zmax-Zmin))

  if((Rmin.lt.r8grid(1)-ztol).or.(Rmax.gt.r8grid(inumR)+ztol).or. &
       (Zmin.lt.z8grid(1)-ztol).or.(Zmax.gt.z8grid(inumZ)+ztol)) then
     write(lunzer(0),*) ' ?? coverage of grid is incomplete.'
     write(lunzer(0),*) '    R grid range: ',r8grid(1),r8grid(inumR)
     write(lunzer(0),*) '    need: ',Rmin,Rmax
     write(lunzer(0),*) '    Z grid range: ',z8grid(1),z8grid(inumZ)
     write(lunzer(0),*) '    need: ',Zmin,Zmax
     ierr=1
  else
     !if not enough precision be sure xplasma will not give an error
     if(r8grid(1).gt.Rmin) r8grid(1)=Rmin
     if (r8grid(inumR).lt.Rmax) r8grid(inumR)=Rmax
     if(z8grid(1).gt.Zmin) z8grid(1)=Zmin
     if(z8grid(inumZ).lt.Zmax) z8grid(inumZ)=Zmax
     call eqm_rzgrid(r8grid,z8grid,-1,-1,inumR,inumZ,ztol, &
             id_R,id_Z,ierr)
  endif

  return
end subroutine trx_wall_grid
 
subroutine trx_wall_RZ(ilun,Rmin,Rmax,inumR,Zmin,Zmax,inumZ,ierr)
!
!  read from the TRANSP namelist the TRANSP circles & lines
!  limiter specifications, then create this limiter in xplasma,
!  and create an inumR x inumZ cartesian box that covers a
!  rectangular region enclosing both the core plasma and the
!  scrape-off plasma (btw limiter and core plasma)
!
!  MOD DMC Nov. 2006 -- use output database limiter if available...
!
  use trx_module

  implicit NONE
!
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer ilun                ! lun for file read (input)
  real*8 Rmin,Rmax            ! Rmin,Rmax of grid (0.0 for automatic grid)
  real*8 Zmin,Zmax            ! Zmin,Zmax of grid (0.0 for automatic grid)
  integer inumR,inumZ         ! cartesian overlay grid sizes (input)
  integer ierr                ! completion code:  0=OK
!
! **NOTE** Rmin,Rmax,Zmin,Zmax modified on output
!          to give the actual limits used, i.e. taking limiter
!          into account.  If these quantities are non-zero on
!          input, they define a minimum range to be covered,
!          which may be expanded to reach the limiter given
!          in the TRANSP namelist, if necessary.
!
! if Rmin,Rmax,Zmin and Zmax are  all <0, do not set the grid here;
! but DO return the approximate limiter Rmin,Rmax,Zmin,Zmax in these
! variables
!
! mod DMC Jan 2007 -- force Zmin = -Zmax !!
!   required for some physics codes that use geqdsk files generated
!   from xplasma
!
!--------------------------------------------------------------------
!
  integer ncirlm,nlinlm,itype,i,ilims,irank,idims(10)
  real*8, dimension(:), allocatable :: alnlimR,alnlimZ,alnlimt
  real*8, dimension(:), allocatable :: crlimR,crlimZ,crlimrad
  real*8, dimension(:), allocatable :: rpts,ypts
!
  real*8 Rmin_use,Rmax_use,Zmin_use,Zmax_use,zdist
  real*8 zRgrid(-1:inumR+2),zZgrid(-1:inumZ+2)
  real*8 zero,zone,zzdum
  real*8 ratmin,ratmax,zat
  real*8 rp1,yp1,rp2,yp2
!
! SINGLE precision for t1profil interface:
!
  real :: ztdum
  real, dimension(:), allocatable :: r4vals
!
  character*64 zlabel,zxabb
  character*32 zunits
  integer :: imulti,istype,ineed_grid,isize,idum
!
  logical :: nogrid
!
  integer lunzer
!
  data zero/0.0_R8/
  data zone/1.0_R8/
!
!--------------------------------------------------------------------
!    first see if limiter data is in the run output database
!
  nogrid=.FALSE.
  if(max(rmin,rmax,zmin,zmax).lt.zero) nogrid=.TRUE.

  call rplabel('RLIM',zlabel,zunits,imulti,istype)
  if(istype.gt.0) call rplabel('YLIM',zlabel,zunits,imulti,istype)

  if(istype.gt.0) then
     ineed_grid = 1
     call rpdims(istype,irank,idims,zxabb,ierr)
     if(ierr.ne.0) call bad_exit

     isize=idims(1)
     allocate(rpts(isize),ypts(isize),r4vals(isize))

     ztdum=0

     call t1profil('RLIM',zlabel,zunits,ztdum,ztdum, &
          istype,r4vals,isize,idum,ierr)
     if(ierr.ne.0) call bad_exit
     rpts = r4vals*0.01_R8  ! -> m

     call t1profil('YLIM',zlabel,zunits,ztdum,ztdum, &
          istype,r4vals,isize,idum,ierr)
     if(ierr.ne.0) call bad_exit
     ypts = r4vals*0.01_R8  ! -> m

!
! prevent limiter range from extending too far from plasma boundary unless it would cause a failure
!
     call rational_limiter(inumR, ratmin, ratmax, zat, ierr)
     if (ierr.eq.0) then
        rp1 = min(ratmax,max(ratmin,rpts(1)))
        yp1 = min(zat   ,max(-zat,  ypts(1)))

        do i=2,isize
           rp2 = min(ratmax,max(ratmin,rpts(i)))
           yp2 = min(zat   ,max(-zat,  ypts(i)))

           if((rp1.eq.rp2).and.(yp1.eq.yp2)) goto 150  ! eqm_cbdy will fail -- try to avoid this

           rp1=rp2
           yp1=yp2
        end do

        rpts = min(rpts,ratmax)
        rpts = max(rpts,ratmin)
        ypts = min(ypts,zat)
        ypts = max(ypts,-zat)

150     continue
     else
        write(lunzer(0),*) ' !trxplib: unexpected error finding plasma boundary before limiter'  
     end if

     call eqm_cbdy(isize,rpts,ypts,ierr)
     if(ierr.eq.0) then
!
!  fetch (R,Z) limits needed to cover limiters
!
        if(nogrid) then
           call my_eq_bdlims(itype,Rmin,Rmax,Zmin,Zmax,ierr)
        else
           call my_eq_bdlims(itype,Rmin_use,Rmax_use,Zmin_use,Zmax_use,ierr)
        endif
     endif

  else
     write(lunzer(0),*) ' %trxplib: no contour limiter; use namelist...'
!--------------------------------------------------------------------
!    no limiter data in database so...
!
!    read the TRANSP namelist; get limiter / wall information
!
     call tr_getnl_text(ilun,ierr)  ! try reading the namelist
!
!  namelist read completed
!--------------------------------------------------------------------
!  get the limiter info ...
!
!  no. of TRANSP limiters
     call trx_getnlims(ncirlm,nlinlm,ierr)
     if(ierr.ne.0) return
!
!  circular limiters
!
     allocate(crlimR(max(1,ncirlm)))
     allocate(crlimZ(max(1,ncirlm)))
     allocate(crlimrad(max(1,ncirlm)))
     crlimR=0
     crlimZ=0
     crlimrad=0
     if(ncirlm.gt.0) call trx_getcrlim(ncirlm,crlimR,crlimZ,crlimrad,ierr)
     if(ierr.ne.0) go to 200
!
!  line limiters
!
     allocate(alnlimR(max(1,nlinlm)))
     allocate(alnlimZ(max(1,nlinlm)))
     allocate(alnlimt(max(1,nlinlm)))
     alnlimR=0
     alnlimZ=0
     alnlimt=0
     if(nlinlm.gt.0) call trx_getlnlim(nlinlm,alnlimR,alnlimZ,alnlimt,ierr)
     if(ierr.ne.0) go to 100
!
     ilims=nlinlm+ncirlm
     if(ilims.eq.0) then
        call eq_glimrz(zone,Rmin_use,Rmax_use,Zmin_use,Zmax_use,ierr)
        Zmax_use = max(abs(Zmin_use),abs(Zmax_use))
        Zmin_use = -Zmax_use
        if(ierr.ne.0) call bad_exit
        zdist=0.05*(Rmax_use-Rmin_use)
        write(lunzer(0),*) ' %trx_wall: no limiter spec in namelist.'
        write(lunzer(0),*) '  --assuming distance bdy to wall is:  ',zdist
     endif
!
!  OK -- call xplasma setup routine
!
     if((Rmin.eq.zero).and.(Rmax.eq.zero).and. &
          (Zmin.eq.zero).and.(Zmax.eq.zero)) then

        ineed_grid = 0

        if(ilims.gt.0) then
           call eqm_tbdy_grid(nlinlm,alnlimR,alnlimZ,alnlimt, &
                ncirlm,crlimR,crlimZ,crlimrad, &
                inumR,inumZ, &
                id_R,id_Z, ierr)
           call my_eq_bdlims(itype,Rmin,Rmax,Zmin,Zmax,ierr)
        else
           Rmin=Rmin_use-zdist
           Rmax=Rmax_use+zdist
           Zmin=Zmin_use-zdist
           Zmax=Zmax_use+zdist
           call eqm_dbdy_grid(zdist,Rmin,Rmax,Zmin,Zmax, &
                inumR,inumZ, &
                id_R,id_Z, ierr)
        endif
     else
!
!  set up limiters
!
        ineed_grid = 1

        if(ilims.gt.0) then
           call eqm_tbdy(nlinlm,alnlimR,alnlimZ,alnlimt, &
                ncirlm,crlimR,crlimZ,crlimrad, ierr)
           if(ierr.eq.0) then
!
!  fetch (R,Z) limits needed to cover limiters
!
              if(nogrid) then
                 call my_eq_bdlims(itype,Rmin,Rmax,Zmin,Zmax,ierr)
              else
                 call my_eq_bdlims(itype,Rmin_use,Rmax_use,Zmin_use,Zmax_use, &
                      ierr)
              endif
!
           endif
!
        else
!
!  const. dist. limiter
!
           Rmin_use=Rmin_use-zdist
           Rmax_use=Rmax_use+zdist
           Zmin_use=Zmin_use-zdist
           Zmax_use=Zmax_use+zdist
           if(nogrid) then
              Rmin=Rmin_use
              Rmax=Rmax_use
              Zmin=Zmin_use
              Zmax=Zmax_use
           endif
           call eqm_dbdy(zdist,Rmin_use,Rmax_use,Zmin_use,Zmax_use,ierr)
!
        endif
     endif
  endif
!
!  limiter acquired...
!  expand to user specified limits if necessary
!
  if(nogrid) ineed_grid=0
  if(ineed_grid.eq.1) then
     if(ierr.eq.0) then
        if((Rmin.ne.zero).or.(Rmax.ne.zero)) then
           Rmin_use=min(Rmin_use,Rmin)
           Rmax_use=max(Rmax_use,Rmax)
        endif
        if((Zmin.ne.zero).or.(Zmax.ne.zero)) then
           Zmin_use=min(Zmin_use,Zmin)
           Zmax_use=max(Zmax_use,Zmax)
           if (symZ) then
              Zmax_use = max(abs(Zmin_use),abs(Zmax_use))
              Zmin_use = -Zmax_use
           end if
        endif
!
!  generate grids -- 2 zones extra on each side
!
        do i=-1,inumR+2
           zRgrid(i)=Rmin_use+(i-1)*(Rmax_use-Rmin_use)/(inumR-1)
        enddo
        do i=-1,inumZ+2
           zZgrid(i)=Zmin_use+(i-1)*(Zmax_use-Zmin_use)/(inumZ-1)
        enddo
!
!  create the grids in the xplasma module
!
        call eqm_rzgrid(zRgrid,zZgrid,0,0,inumR+4,inumZ+4,1.0E-5_R8, &
             id_R,id_Z,ierr)
!
        Rmin=Rmin_use
        Rmax=Rmax_use
        Zmin=Zmin_use
        Zmax=Zmax_use
!
     endif
  endif

!  cleanup

100 continue
  if(allocated(alnlimR)) then
     deallocate(alnlimR)
     deallocate(alnlimZ)
     deallocate(alnlimt)
  endif

200 continue
  if(allocated(crlimR)) then
     deallocate(crlimR)
     deallocate(crlimZ)
     deallocate(crlimrad)
  endif

  if(allocated(rpts)) deallocate(rpts,ypts,r4vals)

  return

  contains

    subroutine my_eq_bdlims(itype,Rmini,Rmaxi,Zmini,Zmaxi,ierr)
      integer :: itype
      real*8 :: Rmini,Rmaxi,Zmini,Zmaxi
      integer :: ierr

      !  call eq_bdlims; then enforce Zmini = -Zmaxi

      call eq_bdlims(itype,Rmini,Rmaxi,Zmini,Zmaxi,ierr)
      if (symZ) then
         Zmaxi = max(abs(Zmini),abs(Zmaxi))
         Zmini = -Zmaxi
      end if
    end subroutine my_eq_bdlims

    !
    ! expect a real limiter to be within these limits
    !
    subroutine rational_limiter(ir, rrat_min, rrat_max, zrat, ii)
      use xplasma_obj_instance, only: s
      use xplasma_ctran,        only: xplasma_RZminmax_plasma
      implicit none

      real*8, parameter :: RMIN_FRAC = 0.05d0   ! limiter no closer to 0 then RMIN_FRAC*rbdy_min
      real*8, parameter :: RMAX_FAC  = 0.7d0    ! limiter no further from rbdy_max then RMAX_FRAC*rcenter
      real*8, parameter :: ZR_FAC    = 0.7d0    ! limiter no further from |zbdy| then ZR_FAC*rcenter
      
      integer, intent(in)  :: ir        ! number of R grid points
      real*8,  intent(out) :: rrat_min  ! minimum acceptable R of limiter for gridding
      real*8,  intent(out) :: rrat_max  ! maximum acceptable R of limiter for gridding
      real*8,  intent(out) :: zrat      ! maximum acceptable |Z| of limiter for gridding
      integer, intent(out) :: ii        ! nonzero on error

      real*8 :: rm,rp,zm,zp   ! bdy limits
      real*8 :: rc            ! approximate R center

      ii = 0
      call xplasma_RZminmax_plasma(s, rm,rp, zm,zp, ii)
      
      if (ii==0 .and. rm>0.d0 .and. rp>rm .and. zp>zm) then
         rc = (rm+rp)/2.d0
         
         rrat_max = rp + RMAX_FAC*rc
         rrat_min = RMIN_FRAC*rm
         if (ir>1) rrat_min=max(rrat_min, rrat_max*2.5d0/(ir+1))   ! keep ghost zone>0, rrat_min>2*rrat_max/(ir+1)
         
         zrat = max(abs(zm),abs(zp)) + ZR_FAC*rc
      else
         ii = 1
      end if
    end subroutine rational_limiter
end subroutine trx_wall_RZ


