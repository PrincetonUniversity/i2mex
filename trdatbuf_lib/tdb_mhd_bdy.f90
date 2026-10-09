! this collection of routines returns information on the plasma boundary
! based on equilibrium boundary (or entire equilibrium) contained in the
! trdatbuf object passed...

subroutine tdb_bdy_type(d,ibtype)
  use trdatbuf_obj

  ! return the type of boundary data available:
  !   ibtype = 1 -- updown symmetric
  !   ibtype = 0 -- updown asymmetric
  !   ibtype = -1 -- no data

  implicit NONE
  type (trdatbuf) :: d
  integer, intent(out) :: ibtype

  if(d%datbuf(5).gt.0.0d0) then
     ibtype=1
  else if(d%datbuf(5).lt.0.0d0) then
     ibtype=-1
  else
     ibtype=0
  endif

end subroutine tdb_bdy_type

subroutine tdb_bdy_timefac(d,ztime,it,zf)
  use trdatbuf_obj
  use tdbsub_uts  ! private

  !  get time interpolation factor for boundary Fourier moments

  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(in) :: ztime    ! time (seconds)
  integer, intent(out) :: it     ! time bin
  real*8, intent(out) :: zf      ! interpolation factor w/in bin

  !  interpolation of f(t) item stored at datbuf(istart)
  !  will be:  result = d%datbuf(istart+it-1) + &
  !                        zf*d%(datbuf(istart+it)-d%datbuf(istart+it-1))

  integer :: iamoms,ilt,int
  !------------------------------------

  iamoms = max(d%ldmmx,d%lmomd)
  if(iamoms.eq.0) then
     !  must be using R0(t) and a(t) only -- on ltime1 timebase
     ilt = d%ltime1
     int = d%ntime1
     call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,ztime,it,zf)

  else
     !  general moments data on ltime2 timebase
     ilt = d%ltime2
     int= d%ntime2
     call tdbsub_lookup(d%datbuf(ilt:ilt+int-1),int,ztime,it,zf)

  endif

end subroutine tdb_bdy_timefac

real*8 function tdb_bdy_getmom(d,it,zf,imom,itype)
  !
  !  get Fourier moment of plasma boundary
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(in) :: it     ! time bin index
  real*8, intent(in) :: zf      ! time interpolation factor
  integer, intent(in) :: imom   ! moment index
  integer, intent(in) :: itype  ! moment type

  ! TDB_MOMS_RCOS or TDB_MOMS_RSIN or TDB_MOMS_ZCOS or TDB_MOMS_ZSIN
  ! ...function value -- the Fourier moment -- in cm

  integer :: iamoms,jr0,jr1,jr2,jy1,jy2,iadd,ja,iadci,ibtype
  logical :: ilsym

  !------------------------

  iamoms = max(d%ldmmx,d%lmomd)

  call tdb_bdy_type(d,ibtype)
  ilsym = ( ibtype .gt. 0) ! .TRUE. => updown symmetric data stored

  tdb_bdy_getmom = ZERO

  if(ibtype.lt.0) return   ! (no data -- always return ZERO)

  if(iamoms.eq.0) then
     !  must be using R0(t) and a(t) only -- on ltime1 timebase
     if(imom.eq.0) then
        if(itype.eq.TDB_MOMS_RCOS) then
           if(d%ldatpos.gt.0) then
              tdb_bdy_getmom = d%datbuf(d%ldatpos+it-1) + &
                   zf*(d%datbuf(d%ldatpos+it)-d%datbuf(d%ldatpos+it-1))
           endif
        endif
     else if(imom.eq.1) then
        if((itype.eq.TDB_MOMS_RCOS).or.(itype.eq.TDB_MOMS_ZSIN)) then
           if(d%ldatrmn.gt.0) then
              tdb_bdy_getmom = d%datbuf(d%ldatrmn+it-1) + &
                   zf*(d%datbuf(d%ldatrmn+it)-d%datbuf(d%ldatrmn+it-1))
           endif
        endif
     endif

  else if(d%ldmmx.gt.0) then
     !
     ! extracting boundary from full equilibrium dataset
     !
     if(imom.le.d%nmomd) then
        ja = iadci(d, d%ldmmx, it, d%nxmmx, imom, itype)
        tdb_bdy_getmom = d%datbuf(ja) + &
             zf*(d%datbuf(ja+1)-d%datbuf(ja))
     endif

  else
     ! d%ldmom.gt.0 -- general boundary dataset
     if(ilsym) then
        !
        ! only symmetric boundary data was stored
        !
        call momind(d,it,jr0,jr1,jy1)
        if(imom.eq.0) then
           if(itype.eq.TDB_MOMS_RCOS) then
              tdb_bdy_getmom = d%datbuf(jr0) + &
                   zf*(d%datbuf(jr0+1)-d%datbuf(jr0))
           endif
        else if(imom.le.d%nmomd) then
           if(itype.eq.TDB_MOMS_RCOS) then
              iadd=(imom-1)*d%ntime2
              jr1=jr1+iadd
              tdb_bdy_getmom = d%datbuf(jr1) + &
                   zf*(d%datbuf(jr1+1)-d%datbuf(jr1))
           else if(itype.eq.TDB_MOMS_ZSIN) then
              iadd=(imom-1)*d%ntime2
              jy1=jy1+iadd
              tdb_bdy_getmom = d%datbuf(jy1) + &
                   zf*(d%datbuf(jy1+1)-d%datbuf(jy1))
           endif
        endif
     else
        !
        ! asymmetric boundary data was stored
        !
        if(imom.le.d%nmomd) then
           call momind3(d,it,jr1,jr2,jy1,jy2)
           iadd = imom*d%ntime2
           if(itype.eq.TDB_MOMS_RCOS) then
              ja=jr1+iadd
           else if(itype.eq.TDB_MOMS_RSIN) then
              ja=jr2+iadd
           else if(itype.eq.TDB_MOMS_ZCOS) then
              ja=jy1+iadd
           else if(itype.eq.TDB_MOMS_ZSIN) then
              ja=jy2+iadd
           endif
           tdb_bdy_getmom = d%datbuf(ja) + &
                zf*(d%datbuf(ja+1)-d%datbuf(ja))
        endif
     endif
  endif
end function tdb_bdy_getmom



subroutine tdb_bdy_getmoms(d,it,zf,rmcb,ymcb,nmom)
  !
  !  get ALL Fourier moments of plasma boundary
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(in) :: it     ! time bin index
  real*8, intent(in) :: zf      ! time interpolation factor
  integer, intent(in) :: nmom      ! number of moments
  real*8, dimension(0:nmom,2) :: rmcb,ymcb  ! boundary harmonics 0:nmom, cos:sin
  real*8 :: tdb_bdy_getmom

  integer :: im

  do im=0,nmom
     rmcb(im,1)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_RCOS) ! Rcos mom.
     rmcb(im,2)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_RSIN)
     ymcb(im,1)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_ZCOS)
     ymcb(im,2)=tdb_bdy_getmom(d,it,zf,im,TDB_MOMS_ZSIN) ! Zsin mom.
  enddo

  return
end subroutine tdb_bdy_getmoms
