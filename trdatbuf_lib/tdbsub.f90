! collection of private auxilliar subroutines for trdatbuf data lookup,
! averaging and interpolation

real*8 function tdbsub_i1(zt1,zt2,zt,int,zf)

  use tdbsub_uts

  ! time integrate piecewise linear interpolant zf(1:int) at times zt(1:int)
  ! over time range zt1 to zt2

  real*8, intent(in) :: zt1,zt2
  real*8, intent(in) :: zt(int)
  real*8, intent(in) :: zf(int)

  !  zt is assumed to be monotonic increasing.  this is NOT checked.
  !  zt1 is assumed .ge. zt2

  integer :: it1,it2
  real*8 :: zfrac1,zfrac2,zfrac,zdt

  tdbsub_i1=0
  if (zt1.ge.zt(int)) return
  if (zt2.le.zt(1)) return
  if (zt2.le.zt1) return

  if(zt1.le.zt(1)) then
     it1=1
     zfrac1=0
  else
     call tdbsub_lookup(zt,int,zt1,it1,zfrac1)
  endif

  if(zt2.ge.zt(int)) then
     it2=int-1
     zfrac2=1
  else
     call tdbsub_lookup(zt,int,zt2,it2,zfrac2)
  endif

  if(it1.eq.it2) then
     zfrac = (zfrac1+zfrac2)*HALF
     tdbsub_i1 = (zf(it1)+zfrac*(zf(it1+1)-zf(it1)))*(zt2-zt1)
  else
     tdbsub_i1 = ZERO
     do it=it1,it2
        if(it.eq.it1) then
           zdt=zt(it1+1)-zt1
           zfrac=(ONE+zfrac1)*HALF
        else if(it.eq.it2) then
           zdt=zt2-zt(it2)
           zfrac=zfrac2*HALF
        else
           zdt=zt(it+1)-zt(it)
           zfrac=HALF
        endif
        tdbsub_i1=tdbsub_i1 + zdt*(zf(it)+zfrac*(zf(it+1)-zf(it)))
     enddo
  endif

end function tdbsub_i1

real*8 function tdbsub_i2(zt1,zt2,zt,int,zf1,zf2)

  use tdbsub_uts

  ! time integrate piecewise linear interpolant zf1(1:int)*zf2(1:int)
  ! at times zt(1:int) over time range zt1 to zt2

  real*8, intent(in) :: zt1,zt2
  real*8, intent(in) :: zt(int)
  real*8, intent(in) :: zf1(int),zf2(int)

  !  zt is assumed to be monotonic increasing.  this is NOT checked.
  !  zt1 is assumed .ge. zt2

  integer :: it1,it2
  real*8 :: zfrac1,zfrac2,zfrac,zdt,zf12

  tdbsub_i2=0
  if (zt1.ge.zt(int)) return
  if (zt2.le.zt(1)) return
  if (zt2.le.zt1) return

  if(zt1.le.zt(1)) then
     it1=1
     zfrac1=0
  else
     call tdbsub_lookup(zt,int,zt1,it1,zfrac1)
  endif

  if(zt2.ge.zt(int)) then
     it2=int-1
     zfrac2=1
  else
     call tdbsub_lookup(zt,int,zt2,it2,zfrac2)
  endif

  if(it1.eq.it2) then
     zfrac = (zfrac1+zfrac2)*HALF
     zf12 = (zf1(it1)*zf2(it1)+ &
          zfrac*(zf1(it1+1)*zf2(it1+1)-zf1(it1)*zf2(it1)))
     tdbsub_i2 = zf12*(zt2-zt1)
  else
     tdbsub_i2 = ZERO
     do it=it1,it2
        if(it.eq.it1) then
           zdt=zt(it1+1)-zt1
           zfrac=(ONE+zfrac1)*HALF
        else if(it.eq.it2) then
           zdt=zt2-zt(it2)
           zfrac=zfrac2*HALF
        else
           zdt=zt(it+1)-zt(it)
           zfrac=HALF
        endif
        zf12 = zf1(it)*zf2(it)+ &
             zfrac*(zf1(it+1)*zf2(it+1)-zf1(it)*zf2(it))
        tdbsub_i2=tdbsub_i2 + zdt*zf12
     enddo
  endif

end function tdbsub_i2

