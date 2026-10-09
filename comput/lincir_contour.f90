!
! construct a contour from line and circle limiters
!
! The algorithm finds the pairwise intersections between all the limiters which 
! are not outside the other limiters.  Segments between these intersections form the 
! the limiter boundary.
!

module lincir_contour_mod

contains

  subroutine lincir_contour(r0, z0,   &
       nlinlm, alnlmr, alnlmy, alnlmt, &
       ncirlm, crlmr1, crlmy1, crlmrd, &
       rout, zout, ier, cseg, debug,   &
       urmin, urmax, uzmin, uzmax)

    use iso_c_binding, only: fp => c_double
    implicit none

    real(fp), parameter :: LARGE_R = 25._fp   ! bound limiters in R to less then this factor of r0
    real(fp), parameter :: SMALL_R = .0001_fp ! bound limiters in R to more then this factor of r0
    real(fp), parameter :: LARGE_Z = 25._fp   ! bound limiters in Z to less then this magnitude of r0

    real(fp), intent(in) :: r0                 ! reference R inside all limiters
    real(fp), intent(in) :: z0                 ! reference Z inside all limiters

    integer,  intent(in) :: nlinlm             ! number of line limiters
    real(fp), intent(in) :: alnlmr(nlinlm)     ! reference R of line limiter
    real(fp), intent(in) :: alnlmy(nlinlm)     ! reference Z of line limiter
    real(fp), intent(in) :: alnlmt(nlinlm)     ! angle in degrees from -R axis upward for line limiter

    integer,  intent(in) :: ncirlm             ! number of circle limiters
    real(fp), intent(in) :: crlmr1(ncirlm)     ! R of center of circle limiter
    real(fp), intent(in) :: crlmy1(ncirlm)     ! Z of center of circle limiter
    real(fp), intent(in) :: crlmrd(ncirlm)     ! radius of circle limiter

    real(fp), pointer    :: rout(:)  ! will be dimensioned and filled with R of contour
    real(fp), pointer    :: zout(:)  ! will be dimensioned and filled with Z of contour
    integer,  intent(out):: ier      ! nonzero on error

    real(fp), intent(in), optional :: cseg    ! maximum segment size on a circle limiter or 0. to use default
    logical,  intent(in), optional :: debug   ! .true. for some debugging output, specifically a gnuplottable file

    real(fp), intent(in), optional :: urmin   ! make bounding min R limiter at this value, overrides LARGE_R
    real(fp), intent(in), optional :: urmax   ! make bounding max R limiter at this value, overrides SMALL_R
    real(fp), intent(in), optional :: uzmin   ! make bounding min Z limiter at this value, overrides LARGE_Z
    real(fp), intent(in), optional :: uzmax   ! make bounding max Z limiter at this value, overrides -LARGE_Z
    !
    ! -------------------------- types --------------------------
    !
    ! identifies an intersection point and one of the limiters
    !
    type interpoint
      integer :: iptr  ! intersection point index
      integer :: ix    ! index 1 or 2 for limiter involved in intersection
    end type interpoint

    !
    ! stupid fortran
    !
    type pinterpoint
      type(interpoint), pointer :: ptr
    end type pinterpoint

    !
    ! at an intersection point this describes the limiter and where
    ! on the limiter it was intersected.  Also contains an interpoint for
    ! forming a linked list
    !
    type interlim
      integer          :: index   ! limiter index >0 for line, <0 for circle
      real(fp)         :: dist    ! distance along limiter 
      type(interpoint) :: next    ! linked list
    end type interlim

    !
    ! describes the intersection between two limiters which lies on the limiter boundary
    !
    type intersect
      real(fp)       :: r,z      ! coordinates of intersection point
      type(interlim) :: lim(2)   ! limiters involved in intersection

      integer :: seg1    ! segment index for segment containing this point
      integer :: seg2    ! segment index for second segment containing this point
    end type intersect

    !
    ! a segment between two intersect points which lies on the limiter boundary
    !
    type segment
      integer  :: a      ! intersect point
      integer  :: b      ! intersect point
      real(fp) :: r      ! if >0. this is the R of the midpoint on a circle 
      real(fp) :: z      ! if >0. this is the Z of the midpoint on a circle 
    end type segment
    !
    ! ------------------------------------------------------------------------
    !
    logical :: isopen               ! .true. when debug file is open
    logical :: gdebug               ! .true. for debugging

    integer :: nlines               ! number of line limiters
    integer :: ncircs               ! number of circle limiters
    integer :: ntot                 ! total number of limiters
    integer :: maxint               ! maximum number of pairwise intersections
    integer :: i,j,k,is,ic,imin     ! temp
    integer :: nc                   ! number of limiter intersections
    integer :: nseg                 ! number of segments
    integer :: ncont                ! number of contour points
    integer :: ifirst,icur,iprev    ! first intersect point, current and last intersect point added to contour
    integer :: i1(1)                ! stupid fortran
    integer :: idir                 ! step along rcont,zcont in ccw direction
    integer :: imaxR,imaxZ,iminZ    ! rcont,zcont indices which achieve R,Z limits

    real(fp) :: csegx               ! maximum segment size on a circle limiter
    real(fp) :: q                   ! temp
    real(fp) :: det                 ! determinant,descriminant
    real(fp) :: dr,dz,delta         ! delta R,Z,distance
    real(fp) :: rx,zx,rc,zc         ! intersection and computed points
    real(fp) :: li,lj               ! distance along line from reference point
    real(fp) :: pi,pj               ! limiter radii
    real(fp) :: n2i,d2i,n2j,d2j     ! numerator,denominator terms
    real(fp) :: olphai,olphaj       ! circle intersection angles relative to center-center vector
    real(fp) :: dalphai,dalphaj     ! angle of delta vector with respect to each circle limiter
    real(fp) :: alphai(2),alphaj(2) ! circle intersection angles relative to R for two intersection points
    real(fp) :: smid,pmid2          ! distance along line to midpoint of intersections, circle radius of midpoint
    real(fp) :: s(2)                ! distance along line to intersections
    real(fp) :: xpi                 ! constant named for a desert
    real(fp) :: twopi               ! ?
    real(fp) :: rmin,rmax,zmin,zmax ! bounding box of final contour


    real(fp), allocatable :: rline(:),zline(:)  ! point on line limiters
    real(fp), allocatable :: uline(:),vline(:)  ! unit direction vector for line limiter
    real(fp), allocatable :: sline(:),scirc(:)  ! +-1, sign of side of r0,z0

    real(fp), allocatable :: rcirc(:),zcirc(:)  ! center of circle limiter
    real(fp), allocatable :: pcirc(:)           ! radius of circle limiter
    integer,  allocatable :: kcirc(:)           ! count number of intersections with other limiters

    real(fp), allocatable :: sdist(:),sindex(:) ! for one limiter, distance to intersect point and index of intersect point
    real(fp), allocatable :: rcont(:),zcont(:)  ! R,Z contour

    type(interpoint) :: sptr                    ! temp pointer to intersection point

    type(interpoint),  allocatable, target :: ptrline(:), ptrcirc(:) ! point to first intersect point in linked list for limiter

    type(pinterpoint), allocatable         :: ltrline(:), ltrcirc(:) ! points to last interpoint in linked list for limiter
    type(intersect),   allocatable, target :: inter(:)               ! list of limiter intersections which are on boundary
    type(segment),     allocatable :: seg(:)                 ! segments between boundary point intercepts


    ier=0
    if (associated(rout)) deallocate(rout)
    if (associated(zout)) deallocate(zout)

    if (r0<=0._fp) then
      print *, '?lincir_contour: r0<=0. not allowed'
      ier=1 ; return
    end if

    xpi   = 4._fp*atan(1._fp)
    twopi = 2._fp*xpi
    isopen = .false.

    csegx = 0.1_fp*r0
    if (present(cseg)) then
      if (cseg>0._fp) csegx=cseg
    end if

    gdebug=.false.
    if (present(debug)) gdebug=debug

    ! -- allocate and add extra bounding line limiters --
    nlines = nlinlm+4
    ncircs = ncirlm

    ntot = nlines+ncircs

    maxint = nlines*(nlines-1)/2 + ncircs*(ncircs-1) + 2*nlines*ncircs

    allocate(rline(nlines),zline(nlines), uline(nlines), vline(nlines))
    allocate(sline(nlines), ptrline(nlines), ltrline(nlines))
    allocate(rcirc(ncircs), zcirc(ncircs), pcirc(ncircs))
    allocate(scirc(ncircs), ptrcirc(ncircs), kcirc(ncircs),ltrcirc(ncircs))
    allocate(inter(maxint), sdist(maxint), sindex(maxint), seg(maxint))
    allocate(rcont(maxint), zcont(maxint))

    rline(1:nlinlm) = alnlmr(1:nlinlm)
    zline(1:nlinlm) = alnlmy(1:nlinlm)
    uline(1:nlinlm) = -cos(alnlmt(1:nlinlm)*xpi/180._fp)
    vline(1:nlinlm) =  sin(alnlmt(1:nlinlm)*xpi/180._fp)

    rline(nlinlm+1) = LARGE_R*r0
    if (present(urmax)) rline(nlinlm+1)=urmax
    zline(nlinlm+1) = z0
    uline(nlinlm+1) = 0._fp
    vline(nlinlm+1) = 1._fp

    rline(nlinlm+2) = SMALL_R*r0
    if (present(urmin)) rline(nlinlm+2)=max(urmin,SMALL_R*r0)  ! prevent <0
    zline(nlinlm+2) = z0
    uline(nlinlm+2) = 0._fp
    vline(nlinlm+2) = 1._fp

    rline(nlinlm+3) = r0
    zline(nlinlm+3) = z0 + LARGE_Z*r0
    if (present(uzmax)) zline(nlinlm+3) = uzmax
    uline(nlinlm+3) = 1._fp
    vline(nlinlm+3) = 0._fp

    rline(nlinlm+4) = r0
    zline(nlinlm+4) = z0 - LARGE_Z*r0
    if (present(uzmin)) zline(nlinlm+4) = uzmin
    uline(nlinlm+4) = 1._fp
    vline(nlinlm+4) = 0._fp

    rcirc = crlmr1(1:ncirlm)
    zcirc = crlmy1(1:ncirlm)
    pcirc = crlmrd(1:ncirlm)
    kcirc = 0

    ! -- pick a side of the limiters --
    do i=1, nlines
      q = -(r0-rline(i))*vline(i) + (z0-zline(i))*uline(i)
      if (q==0._fp) then
        print *, '?lincir_contour: r0,z0 lies on line limiter ',i
        ier=1 ; goto 500
      end if
      sline(i) = sign(1.0_fp,q)
      ptrline(i)%iptr=0
      ptrline(i)%ix  =0
      ltrline(i)%ptr => ptrline(i)
    end do

    do i=1, ncircs
      q = ((r0-rcirc(i))**2 + (z0-zcirc(i))**2) - pcirc(i)**2
      if (q==0._fp) then
        print *, '?lincir_contour: r0,z0 lies on circle limiter ',i
        ier=1 ; goto 500
      end if
      scirc(i) = sign(1.0_fp,q)
      ptrcirc(i)%iptr=0
      ptrcirc(i)%ix  =0
      ltrcirc(i)%ptr => ptrcirc(i)
    end do

    if (gdebug) then
      !
      ! open debug file and write out lines and circle limiters
      !
      open(unit=47,file="lincir.dat",status="replace")
      isopen = .true.

      write(47,'(a)') "# plot 'lincir.dat' index 0 title 'limiters', '' index 2 lt 3 " &
           // "title 'contour', '' index 1 with points pt 4 ps 2 title 'intersects'"

      do i=1, nlines
        write(47,'(a,i3)') '# line limiter ',i
        q = 2._fp*max(r0, rline(i), abs(zline(i)))
        write(47,'(2(es14.6,1x))') rline(i)-q*uline(i), zline(i)-q*vline(i)
        write(47,'(2(es14.6,1x))') rline(i), zline(i)
        write(47,'(2(es14.6,1x))') rline(i)+q*uline(i), zline(i)+q*vline(i)
        write(47,'(a)') ' '
      end do

      do i=1, ncircs
        write(47,'(a,i3)') '# circle limiter ',i
        is = max(6,int(twopi*pcirc(i)/csegx))
        do k=0,is
          dr = rcirc(i) + pcirc(i)*cos((twopi*k)/is)
          dz = zcirc(i) + pcirc(i)*sin((twopi*k)/is)
          write(47,'(2(es14.6,1x))') dr, dz
        end do
        write(47,'(a)') ' '
      end do
    end if

    ! -- find intersections on boundary --
    nc = 0                      ! count number of intersections

    do i = 1, nlines-1          ! line-line
      do j = i+1, nlines
        det = vline(i)*uline(j)-uline(i)*vline(j)
        if (abs(det)<1.e-10_fp) cycle   ! parallel lines

        dr  = rline(j) - rline(i)
        dz  = zline(j) - zline(i)
        li  = (-vline(j)*dr + uline(j)*dz)/det
        lj  = (-vline(i)*dr + uline(i)*dz)/det

        rx = rline(i) + li*uline(i)
        zx = zline(i) + li*vline(i)

        rc = rline(j) + lj*uline(j)
        zc = zline(j) + lj*vline(j)

        if (max(abs(rx-rc),abs(zx-zc))>1.e-8_fp*r0) then
          print *, '?lincir_contour: logic error in line intersection'
          ier=1 ; goto 500
        end if

        if (.not. inside_limiter(rx,zx,i,j)) cycle   ! only the boundary intersections

        nc = nc+1
        if (nc>maxint) goto 490

        inter(nc)%r = rx
        inter(nc)%z = zx
        inter(nc)%seg1 = 0
        inter(nc)%seg2 = 0
        inter(nc)%lim(1)%index     = i
        inter(nc)%lim(1)%dist      = li
        inter(nc)%lim(1)%next%iptr = 0
        inter(nc)%lim(1)%next%ix   = 0
        inter(nc)%lim(2)%index     = j
        inter(nc)%lim(2)%dist      = lj
        inter(nc)%lim(2)%next%iptr = 0
        inter(nc)%lim(2)%next%ix   = 0

        ltrline(i)%ptr%iptr = nc
        ltrline(i)%ptr%ix   = 1
        ltrline(i)%ptr      => inter(nc)%lim(1)%next

        ltrline(j)%ptr%iptr = nc
        ltrline(j)%ptr%ix   = 2
        ltrline(j)%ptr      => inter(nc)%lim(2)%next
      end do
    end do

    do i = 1, ncircs-1          ! circle-circle
      do j = i+1, ncircs
        dr = rcirc(j) - rcirc(i)
        dz = zcirc(j) - zcirc(i)

        delta = sqrt(dr**2+dz**2)
        pi = pcirc(i)
        pj = pcirc(j)

        if ((pi+pj)<delta .or. delta<abs(pi-pj)) cycle  ! circles do not intersect

        n2i = pj**2 - (delta-pi)**2
        d2i = (delta+pi)**2 - pj**2

        n2j = pi**2 - (delta-pj)**2
        d2j = (delta+pj)**2 - pi**2

        if (n2i<0._fp .or. d2i<0._fp .or. n2j<0._fp .or. d2j<0._fp) cycle  ! do not intersect numerically

        olphai = 2._fp*atan2(sqrt(n2i),sqrt(d2i))   ! angles relative to delta
        olphaj = 2._fp*atan2(sqrt(n2j),sqrt(d2j))

        dalphai = atan2(dz,dr)      ! theta of delta
        dalphaj = dalphai + xpi     ! atan2(-dz,-dr)
        if (dalphaj>=twopi) dalphaj=dalphaj-twopi

        alphai(1) = dalphai - olphai
        alphaj(1) = dalphaj + olphaj

        alphai(2) = dalphai + olphai
        alphaj(2) = dalphaj - olphaj

        do k = 1, 2                 ! loop over intersection points
          if (k==2 .and. min(abs(olphai),abs(olphaj))<1.e-8_fp) then
            print *, '!lincir_contour: tangential circle-circle limiters'
            cycle 
          end if

          rx = rcirc(i) + pi*cos(alphai(k))
          zx = zcirc(i) + pi*sin(alphai(k))

          rc = rcirc(j) + pj*cos(alphaj(k))
          zc = zcirc(j) + pj*sin(alphaj(k))

          if (max(abs(rx-rc),abs(zx-zc))>1.e-8_fp*r0) then
            print *, '?lincir_contour: logic error in circle-circle intersection'
            ier=1 ; goto 500
          end if

          if (.not. inside_limiter(rx,zx,-i,-j)) cycle   ! only the boundary intersections

          nc = nc+1
          if (nc>maxint) goto 490

          inter(nc)%r    = rx
          inter(nc)%z    = zx
          inter(nc)%seg1 = 0
          inter(nc)%seg2 = 0
          inter(nc)%lim(1)%index     = -i
          inter(nc)%lim(1)%dist      = alphai(k)
          inter(nc)%lim(1)%next%iptr = 0
          inter(nc)%lim(1)%next%ix   = 0
          inter(nc)%lim(2)%index     = -j
          inter(nc)%lim(2)%dist      = alphaj(k)
          inter(nc)%lim(2)%next%iptr = 0
          inter(nc)%lim(2)%next%ix   = 0

          ltrcirc(i)%ptr%iptr = nc
          ltrcirc(i)%ptr%ix   = 1
          ltrcirc(i)%ptr      => inter(nc)%lim(1)%next

          ltrcirc(j)%ptr%iptr = nc
          ltrcirc(j)%ptr%ix   = 2
          ltrcirc(j)%ptr      => inter(nc)%lim(2)%next

          kcirc(j) = kcirc(j)+1
          kcirc(i) = kcirc(i)+1
        end do
      end do
    end do

    do i = 1, nlines          ! line-circle
      do j = 1, ncircs
        dr = rcirc(j) - rline(i)
        dz = zcirc(j) - zline(i)

        smid  = uline(i)*dr + vline(i)*dz
        pmid2 = max(0._fp,dr**2+dz**2-smid**2)    ! should be >0.
        det   = pcirc(j)**2 - pmid2

        if (det<0._fp) cycle   ! does not intersect
        det = sqrt(det)

        s(1) = smid - det
        s(2) = smid + det

        do k = 1, 2            ! loop over two intersection points
          if (k==2 .and. abs(s(2)-s(1))<1.e-8_fp*r0) then
            print *, '!lincir_contour: tangential line-circle limiters'
            cycle 
          end if

          rx = rline(i) + s(k)*uline(i)
          zx = zline(i) + s(k)*vline(i)

          delta = (rx-rcirc(j))**2 + (zx-zcirc(j))**2 - pcirc(j)**2

          if (abs(delta)>1.e-8_fp*r0) then
            print *, '?lincir_contour: logic error in line-circle intersection'
            ier=1 ; goto 500
          end if

          if (.not. inside_limiter(rx,zx,i,-j)) cycle   ! only the boundary intersections

          nc = nc+1
          if (nc>maxint) goto 490

          inter(nc)%r    = rx
          inter(nc)%z    = zx
          inter(nc)%seg1 = 0
          inter(nc)%seg2 = 0
          inter(nc)%lim(1)%index     = i
          inter(nc)%lim(1)%dist      = s(k)
          inter(nc)%lim(1)%next%iptr = 0
          inter(nc)%lim(1)%next%ix   = 0
          inter(nc)%lim(2)%index     = -j
          inter(nc)%lim(2)%dist      = atan2(zx-zcirc(j),rx-rcirc(j))
          inter(nc)%lim(2)%next%iptr = 0
          inter(nc)%lim(2)%next%ix   = 0

          ltrline(i)%ptr%iptr = nc
          ltrline(i)%ptr%ix   = 1
          ltrline(i)%ptr      => inter(nc)%lim(1)%next

          ltrcirc(j)%ptr%iptr = nc
          ltrcirc(j)%ptr%ix   = 2
          ltrcirc(j)%ptr      => inter(nc)%lim(2)%next

          kcirc(j) = kcirc(j)+1
        end do
      end do
    end do

    if (gdebug) then
      if (nc>0) write(47,'(/a)') '# intersect points '
      do i = 1, nc
        write(47,'(2(es14.6,1x))') inter(i)%r, inter(i)%z
      end do
      write(47,'(a)') ' '
    end if

    !
    ! --- find segments ----
    ! A segment is on the limter boundary if its midpoint is inside all the other limiters
    !
    nseg = 0

    do i = 1, nlines
      sptr = ptrline(i)
      do is = 1, maxint
        if (sptr%iptr==0) goto 100
        if ((sptr%ix/=1 .and. sptr%ix/=2) .or. sptr%iptr<0 .or. sptr%iptr>nc) then
          print *, '?lincir_contour: logic error, bad limiter or point index'
          ier=1 ; goto 500
        end if
        if (inter(sptr%iptr)%lim(sptr%ix)%index /= i) then
          print *, '?lincir_contour: logic error, mismatch in line linked list'
          ier=1 ; goto 500
        end if

        sdist(is)  = inter(sptr%iptr)%lim(sptr%ix)%dist
        sindex(is) = sptr%iptr
        sptr       = inter(sptr%iptr)%lim(sptr%ix)%next
      end do
      print *, '?lincir_contour: too many intersections'
      ier=1 ; goto 500

100   continue
      is = is-1
      if (is>1) then
        call tr_dsort(sdist,sindex,is,2)
        do k=1,is-1      ! loop over segments on this limiter
          smid = (sdist(k)+sdist(k+1))/2._fp

          rx = rline(i) + smid*uline(i)   ! midpoint, will test to see if segment is inside or outside all limiters
          zx = zline(i) + smid*vline(i)

          if (inside_limiter(rx,zx,i,0)) then
            call addseg(int(sindex(k)+0.1_fp), int(sindex(k+1)+0.1_fp), 0._fp, 0._fp)  ! it's a boundary segment
            if (ier/=0) goto 500
          end if
        end do
      end if
    end do


    do i = 1, ncircs
      sptr = ptrcirc(i)
      do is = 1, maxint-1             ! will add one more point later
        if (sptr%iptr==0) goto 120
        if ((sptr%ix/=1 .and. sptr%ix/=2) .or. sptr%iptr<0 .or. sptr%iptr>nc) then
          print *, '?lincir_contour: logic error, bad limiter or point index'
          ier=1 ; goto 500
        end if
        if (inter(sptr%iptr)%lim(sptr%ix)%index /= -i) then
          print *, '?lincir_contour: logic error, mismatch in circle linked list'
          ier=1 ; goto 500
        end if

        sdist(is)  = inter(sptr%iptr)%lim(sptr%ix)%dist
        sindex(is) = sptr%iptr
        sptr       = inter(sptr%iptr)%lim(sptr%ix)%next
      end do
      print *, '?lincir_contour: too many intersections'
      ier=1 ; goto 500

120   continue
      is = is-1
      if (is>0) then
        imin = 1
        do k = 1, is
          do while (sdist(k)>=twopi) 
            sdist(k) = sdist(k)-twopi
          end do
          do while(sdist(k)<0._fp)
            sdist(k) = sdist(k)+twopi
          end do
          if (sdist(k)<sdist(imin)) imin=k
        end do

        is = is+1
        sdist(is)  = sdist(imin)+twopi   ! closes last segment on circle
        sindex(is) = sindex(imin)

        call tr_dsort(sdist,sindex,is,2)
        do k=1,is-1      ! loop over segments on this limiter
          smid = (sdist(k)+sdist(k+1))/2._fp

          rx = rcirc(i) + pcirc(i)*cos(smid)   ! midpoint
          zx = zcirc(i) + pcirc(i)*sin(smid)

          if (inside_limiter(rx,zx,-i,0)) then
            call addseg(int(sindex(k)+0.1_fp), int(sindex(k+1)+0.1_fp), rx, zx)
            if (ier/=0) goto 500
          end if
        end do
      end if
    end do

    ! --- build contour ---

    ncont = 0
    rcont = 0._fp
    zcont = 0._fp

    if (nseg>0) then
      is     = 1
      ifirst = seg(is)%a
      iprev  = ifirst
      icur   = 0
      call addcont(ifirst) ; if (ier/=0) goto 500
      do 
        if (seg(is)%a==iprev) then        ! current point is other point in segment
          icur = seg(is)%b
        else if (seg(is)%b==iprev) then
          icur = seg(is)%a
        else
          print *, '?lincir_contour: logic error in segment -> point order'
          ier=1 ; goto 500
        end if
        if (icur<=0 .or. icur>nc) then
          print *, '?lincir_contour: unset point in segment'
          ier=1 ; goto 500
        end if

        if (seg(is)%r>0._fp) then
          ! add circular segment which goes from iprev to icur through seg%r,seg%z
          !call addcont(r=seg(is)%r,z=seg(is)%z)
          call addcirseg(                      &
               inter(iprev)%r, inter(iprev)%z, &
               seg(is)%r,      seg(is)%z,      &
               inter(icur)%r,  inter(icur)%z)
          if (ier/=0) goto 500
        end if

        if (icur==ifirst) exit           ! yeah, it worked!

        call addcont(icur) ; if (ier/=0) goto 500

        if (inter(icur)%seg1==is) then   ! next segment is other segment attached to point
          is = inter(icur)%seg2
        else if (inter(icur)%seg2==is) then
          is = inter(icur)%seg1
        else
          print *, '?lincir_contour: logic error in point -> segment order'
          ier=1 ; goto 500
        end if
        if (is<=0 .or. is>nseg) then
          print *, '?lincir_contour: unset segment in point'
          ier=1 ; goto 500
        end if
        iprev = icur
      end do
    else
      !
      ! no intersections so find innermost circle containing r0,z0
      !
      if (ncircs==0) then
        print *, '?lincir_contour: no segments and no circles?'
        ier=1 ; goto 500
      end if

      imin = 1        ! innermost circle
      det  = 1.e20_fp
      do i=1,ncircs
        q = pcirc(i) - sqrt((r0-rcirc(i))**2+(z0-zcirc(i))**2) ! distance to edge
        if (q>=0._fp .and. q<det) then
          imin = i
          det  = q
        end if
      end do

      k = max(6,int(twopi*pcirc(imin)/csegx))
      do j=0,k
        rx = rcirc(imin) + pcirc(imin)*cos((twopi*j)/k)
        zx = zcirc(imin) + pcirc(imin)*sin((twopi*j)/k)

        call addcont(r=rx,z=zx) ; if (ier/=0) goto 500
      end do

      kcirc(imin) = kcirc(imin)+1     ! so this circle is not checked against contour
    end if

    ! -- copy to output arrays --
    i1 = maxloc(rcont(1:ncont)) ; imaxR = i1(1)
    i1 = maxloc(zcont(1:ncont)) ; imaxZ = i1(1)
    i1 = minloc(zcont(1:ncont)) ; iminZ = i1(1)

    if (imaxZ<imaxR) imaxZ=imaxZ+ncont
    if (iminZ<imaxR) iminZ=iminZ+ncont

    if (iminZ<imaxZ) then
      idir = -1    ! cw so change to ccw
    else
      idir =  1    ! ccw
    end if

    allocate(rout(ncont+1),zout(ncont+1))            ! one extra point to close the contour
    do i = 1, ncont+1
      ic = mod(ncont+(i-1)*idir+imaxR-1,ncont)+1    ! start from larges R and go ccw looping around array ends
      rout(i) = rcont(ic)
      zout(i) = zcont(ic)
    end do

    if (gdebug) then
      write(47,'(/a)') '# contour'
      do i = 1, ncont+1
        write(47,'(2(es14.6,1x))') rout(i),zout(i)
      end do
      write(47,'(a)') ' '
    end if

    rmax = maxval(rout)
    rmin = minval(rout)
    zmax = maxval(zout)
    zmin = minval(zout)

    do i = 1, ncircs
      if (kcirc(i)==0) then
        !
        ! this circle does not have any intersections.  If it is outside all the other limiters
        ! then it can be ignored.  If it is inside the other limiters and contains r0,z0 then it should 
        ! have been caught as the limiter boundary.  If it does not contain r0,z0 but is
        ! inside all the other limiters then it is a floater.
        !
        k = 50     ! check a few points to see if they are inside the other limiters
        ic = 0     ! number of circle points which are inside the other limiters
        do j = 0,k
          rx = rcirc(i)+pcirc(i)*cos((twopi*j)/k)
          zx = zcirc(i)+pcirc(i)*sin((twopi*j)/k)
          if (inside_limiter(rx,zx,-i,0)) ic = ic+1
        end do

        if (ic==0) then
          cycle                  ! probably outside limiter boundary
        else if (ic==(k+1)) then
          if (scirc(i)<0._fp) then
            print *, '?lincir_contour: logic error, inner circle without intersections containing r0,z0 not the limiter?'
          else
            print *, '?lincir_contour: circle limiter ',i,' does not intersect any other limiters'
            print *, '                 but floats inside other limiters and does not contain r0,z0'
          end if
          ier=1 ; goto 500
        else if (ic>0) then
          print *, '?lincir_contour: logic error, circle limiter without intersections has points inside and outside contour'
          ier=1 ; goto 500
        end if
      end if
    end do

    goto 500     ! ok exit

490 continue
    print *, '?lincir_contour: logic error, exceeded intersection storage'
    ier=1 ; goto 500

500 continue
    if (isopen) then
      if (ier==0) then
        write(47,'(a)') '# Ok'
      else
        write(47,'(a)') '# ERROR'
      end if
      close(47)
    end if

    deallocate(rline, zline, uline, vline)
    deallocate(sline, ptrline, ltrline)
    deallocate(rcirc, zcirc, pcirc)
    deallocate(scirc, ptrcirc, kcirc, ltrcirc)
    deallocate(inter, sdist, sindex, seg)
    deallocate(rcont, zcont)

  contains
    !
    ! ------ addcirseg -------
    ! add a circular segment to the contour list
    !
    subroutine addcirseg(r0,z0,r1,z1,r2,z2)
      real(fp),intent(in) :: r0,z0,r1,z1,r2,z2  ! threee points on segment in order of traversal

      integer  :: idir              ! direction of theta for arc from 0 to 1 to 2
      integer  :: kt                ! number of segments to put on an arc
      integer  :: k                 ! temp

      real(fp) :: rm,zm,rp,zp       ! midpoints of line segments connecting first two and last two points
      real(fp) :: drm,dzm,drp,dzp   ! unit vectors orthogonl to line segments
      real(fp) :: dm,dp             ! length of line segments
      real(fp) :: det               ! temp
      real(fp) :: lm,lp             ! distance from midpoints to circle center
      real(fp) :: rc,zc             ! center of circle containing three points
      real(fp) :: dr0,dr1,dr2       ! delta R from center to points
      real(fp) :: dz0,dz1,dz2       ! delta Z from center to points
      real(fp) :: p0,p1,p2          ! radius from center to points
      real(fp) :: t0,t1,t2          ! thetas of three points
      real(fp) :: dt                ! delta theta

      ! -- find center --
      rm = (r0+r1)/2._fp  ; drm = -(z1-z0)
      zm = (z0+z1)/2._fp  ; dzm =  (r1-r0)
      rp = (r1+r2)/2._fp  ; drp = -(z2-z1)
      zp = (z1+z2)/2._fp  ; dzp =  (r2-r1)

      dm = sqrt(drm**2+dzm**2)
      dp = sqrt(drp**2+dzp**2)
      if (dm<1.e-8_fp*r0 .or. dp<1.e-8_fp*r0) then
        print *, '?lincir_contour: coincident points on circular segment'
        ier=1 ; return
      end if

      drm = drm/dm ; dzm = dzm/dm
      drp = drp/dp ; dzp = dzp/dp

      det = dzm*drp-drm*dzp
      if (abs(det)<1.e-8_fp) then
        print *, '?lincir_contour: linear points on circular segment'
        ier=1 ; return
      end if

      lm = (-dzp*(rp-rm)+drp*(zp-zm))/det
      lp = (-dzm*(rp-rm)+drm*(zp-zm))/det

      rc = rm + lm*drm
      zc = zm + lm*dzm

      ! -- establish directions --
      dr0 = r0-rc ; dz0 = z0-zc ; p0 = sqrt(dr0**2+dz0**2)
      dr1 = r1-rc ; dz1 = z1-zc ; p1 = sqrt(dr1**2+dz1**2)
      dr2 = r2-rc ; dz2 = z2-zc ; p2 = sqrt(dr2**2+dz2**2)

      if (abs(p0-p1)>1.e-8*r0 .or. abs(p0-p2)>1.e-8*r0) then
        print *, '?lincir_contour: logic error, bad center'
        ier=1 ; return
      end if

      det = dr0*dz1-dr1*dz0    ! (r0,z0) X (r1,z1)  establishes arc direction
      if (det>0._fp) then
        idir=1
      else
        idir=-1
      end if

      ! -- arc from r0,z0 to r1,z1 --
      t0 = atan2(dz0,dr0)
      t1 = atan2(dz1,dr1)
      if (idir*t1<idir*t0) t1=t1+idir*twopi
      kt = max(1,int(p0*abs(t1-t0)/csegx))
      dt = idir*abs(t1-t0)/kt
      do k = 1, kt-1
        call addcont(r=rc+p0*cos(t0+k*dt),z=zc+p0*sin(t0+k*dt))
        if (ier/=0) return
      end do

      call addcont(r=r1,z=z1) ; if (ier/=0) return

      ! -- arc from r1,z1 to r2,z2 --
      t1 = atan2(dz1,dr1)
      t2 = atan2(dz2,dr2)
      if (idir*t2<idir*t1) t2=t2+idir*twopi
      kt = max(1,int(p0*abs(t2-t1)/csegx))
      dt = idir*abs(t2-t1)/kt
      do k = 1, kt-1
        call addcont(r=rc+p1*cos(t1+k*dt),z=zc+p1*sin(t1+k*dt)) 
        if (ier/=0) return
      end do
    end subroutine addcirseg


    !
    ! ------ addcont -------
    ! add a point to the contour either by intersection point index or by R,Z
    subroutine addcont(i, r, z)
      integer, intent(in),optional :: i    ! intersection point index
      real(fp),intent(in),optional :: r,z  ! point to add

      real(fp) :: rs,zs    ! actual point to add
      real(fp),allocatable :: t(:)  ! temp

      if (present(i) .and. (present(r).or.present(z))) then
        print *, '?lincir_contour: too many arguments to addcont'
        ier=1 ; return
      end if

      if (present(i)) then
        rs = inter(i)%r
        zs = inter(i)%z
      else if (present(r) .and. present(z)) then
        rs = r
        zs = z
      else
        print *, '?lincir_contour: only r or z given to addcont'
        ier=1 ; return
      end if

      ncont = ncont+1
      if (ncont>size(rcont)) then
        allocate(t(size(rcont)))        ! expand rcont,zcont
        t = rcont
        deallocate(rcont)
        allocate(rcont((13*ncont)/10))
        rcont(1:size(t)) = t

        t = zcont
        deallocate(zcont)
        allocate(zcont(size(rcont)))
        zcont(1:size(t)) = t
        deallocate(t)
      end if

      rcont(ncont) = rs
      zcont(ncont) = zs
    end subroutine addcont


    !
    ! ----- addseg -----
    ! create a new segment on the limiter boundary
    !
    subroutine addseg(a,b,r,z)
      integer, intent(in) :: a     ! intersect index of one point on segment
      integer, intent(in) :: b     ! intersect index of second point on segment
      real(fp),intent(in) :: r,z   ! >0. then this is the circle midpoint

      nseg = nseg+1
      if (nseg>size(seg)) then
        print *, '?lincir_contour: too many segments'
        ier=1 ; return
      end if

      seg(nseg)%a = a
      seg(nseg)%b = b
      seg(nseg)%r = r
      seg(nseg)%z = z

      if (inter(a)%seg1==0) then
        inter(a)%seg1=nseg
      else if (inter(a)%seg2==0) then
        inter(a)%seg2=nseg
      else
        goto 25
      end if

      if (inter(b)%seg1==0) then
        inter(b)%seg1=nseg
      else if (inter(b)%seg2==0) then
        inter(b)%seg2=nseg
      else
        goto 25
      end if
      return

25    continue
      print *, '?lincir_contour: logic error, too many segments attach to one point'
      ier=1 
    end subroutine addseg

    !
    ! ------------ inside_limiter -----------
    ! return true if the r,z point is inside or on the limiter
    ! iex,jex >0 indices of line limiters to skip
    ! iex,jex <0 -indices of circle limiters to skip
    !
    function inside_limiter(r,z,iex,jex) result (isin)
      real(fp),  intent(in) :: r,z      ! coordinates of point
      integer, intent(in) :: iex,jex  ! exclude these limiters 
      logical             :: isin     ! .true. if inside all tested limiters

      integer  :: i
      real(fp) :: q

      isin = .false.
      do i=1, nlines
        if (i==iex .or. i==jex) cycle
        q = -(r-rline(i))*vline(i) + (z-zline(i))*uline(i)
        if (sline(i)*q<0._fp) return
      end do

      do i=1, ncircs
        if (i==-iex .or. i==-jex) cycle
        q = ((r-rcirc(i))**2 + (z-zcirc(i))**2) - pcirc(i)**2
        if (scirc(i)*q<0._fp) return
      end do

      isin = .true.
    end function inside_limiter

  end subroutine lincir_contour

end module lincir_contour_mod

