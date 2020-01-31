subroutine i2mex_getOriNt1(nt1, ier)
  
  ! Return original number of poloidal rays + 1 in res

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1
  integer, intent(out) :: ier

  ier = 0
  call eq_ngrid(i2mex_o%id_t, nt1)

end subroutine i2mex_getOriNt1

subroutine i2mex_getOriNs(ns, ier)
  
  ! Return original number of radial points in res

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  integer, intent(out) :: ier

  ier = 0
  call eq_ngrid(i2mex_o%id_s, ns)

end subroutine i2mex_getOriNs

subroutine i2mex_getOriT(nt1, res, ier)
  
  ! Return original poloidal grid in res
  ! *NOTE* use i2mex_getoriNt1 to get nt1

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1
  real(i2mex_r8), intent(out) :: res(nt1)
  integer, intent(out) :: ier

  real(i2mex_r8) zt(nt1)
  integer iok, idum, i

  ier = 0
  call eq_grid(i2mex_o%id_t, zt, nt1, idum, iok)

  ! invert order if necessary
  if(zt(nt1) > zt(1)) then
     res(1:nt1) = zt(1:nt1)
  else
     do i = 1, nt1
        RES(I) = ZT(NT1+1-I)
     enddo
  endif

  if(idum > nt1) then
     ier = 65
  endif

end subroutine i2mex_getOriT
  
subroutine i2mex_getOriPsi(ns, res, ier)
  
  ! Return original radial poloidal flux grid in res
  ! *NOTE* use i2mex_getoriNs to get ns

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok, idum

  ier = 0
  call eq_grid(i2mex_o%id_s, s, ns, idum, iok)
  call eq_rgetf(idum, s, i2mex_o%id_psi, 0, res, iok)

  if(idum > ns) then
     ier = 66
  endif

end subroutine i2mex_getOriPsi
  
subroutine i2mex_getOriP(ns, res, ier)
  
  ! Return pressure profile on original grid
  ! *NOTE* use i2mex_getoriNs to get ns

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok, idum

  ier = 0
  call eq_grid(i2mex_o%id_p, res, ns, idum, iok)
  if(idum > ns) then
     ier = 67
  endif

end subroutine i2mex_getOriP
  
subroutine i2mex_getOriG(ns, res, ier)
  
  ! Return g profile on original grid
  ! *NOTE* use i2mex_getoriNs to get ns

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  real(i2mex_r8) :: s(ns)
  integer iok, idum

  ier = 0
  call eq_grid(i2mex_o%id_g, res, ns, idum, iok)
  if(idum > ns) then
     ier = 68
  endif

end subroutine i2mex_getOriG
  
subroutine i2mex_getOriQ(ns, res, ier)
  
  ! Return q profile on original grid
  ! *NOTE* use i2mex_getoriNs to get ns

  use i2mex_mod
  implicit none
  integer, intent(in) :: ns
  real(i2mex_r8), intent(out) :: res(ns)
  integer, intent(out) :: ier

  integer iok, idum

  ier = 0
  call eq_grid(i2mex_o%id_q, res, ns, idum, iok)
  if(idum > ns) then
     ier = 69
  endif

end subroutine i2mex_getOriQ
  
subroutine i2mex_getOriX(nt1, ns, res, ier)

  ! Return X in res on original grid
  ! *NOTE* use i2mex_getOriNt1 and  i2mex_getOriNs to obtain NT1 and NS
  ! respectively.

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(out) :: res(nt1, ns)
  integer, intent(out) :: ier
  
  real(i2mex_r8) t(nt1), psi(ns)
  integer nti1, nsi, iok

  ier=0
  call i2mex_getOriNt1(nti1, iok)
  if(nti1 /= nt1) then
     ier = 70 
     return
  endif

  call i2mex_getOriNs(nsi, iok)
  if(nsi /= ns) then
     ier = 71 
     return
  endif

  call i2mex_getOriT(nt1, t, ier)
  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_getX(nt1, ns, t, psi, res, iok)
  if(iok/=0) ier = 72

end subroutine i2mex_getOriX
  
subroutine i2mex_getOriZ(nt1, ns, res, ier)

  ! Return Z in res on original grid
  ! *NOTE* use i2mex_getOriNt1 and  i2mex_getOriNs to obtain NT1 and NS
  ! respectively.

  use i2mex_mod
  implicit none
  integer, intent(in) :: nt1, ns
  real(i2mex_r8), intent(out) :: res(nt1, ns)
  integer, intent(out) :: ier
  
  real(i2mex_r8) t(nt1), psi(ns)
  integer nti1, nsi, iok

  ier=0
  call i2mex_getOriNt1(nti1, iok)
  if(nti1 /= nt1) then
     ier = 73
     return
  endif

  call i2mex_getOriNs(nsi, iok)
  if(nsi /= ns) then
     ier = 74 
     return
  endif

  call i2mex_getOriT(nt1, t, ier)
  call i2mex_getOriPsi(ns, psi, ier)
  call i2mex_getZ(nt1, ns, t, psi, res, iok)
  if(iok/=0) ier = 75

end subroutine i2mex_getOriZ
  
  
