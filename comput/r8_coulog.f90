subroutine r8_coulog(et,at,zt,denb,tempb,ab,zb,cln,nb,Bmag,llstix)
  !
  ! (see the warning below.)
  !
  ! calculates the coulomb logarithm for a test particle of energy et
  ! (kev), atomic mass at and atomic charge zt, colliding with NB
  ! maxwellian background species with densities denb (/m**3),
  ! temperatures tempb (keV), atomic masses ab, atomic charges zb.  The
  ! electrons should be included as one of the background species, with
  ! Z=-1 and A=(1./1836.1).  The list of coulomb logarithms is returned
  ! in cln.  The magnetic field should be input as bmag (Tesla).
  !
  ! Note that the test particle loglambda's produced by this subroutine
  ! are not symmetric, i.e. it is not true that
  ! loglam_alpha_beta=loglam_beta_alpha.  If such symmetry is desired,
  ! one should set Et=(3/2)T_t and call this routine twice with species
  ! alpha and beta switched, and take the average value of log_lambda
  ! which is returned.
  !
  ! Written by Greg Hammett, 16-May-1989.
  !
  ! The exact details of the coulomb logarithm are probably not worth
  ! worrying about too much.  The recipe I use here is similiar to the
  ! one given in the NRL plasma formulary.  So it includes the quantum
  ! mechanics correction which is usually important for collisions with
  ! electrons.  Rather than completely ignoring the Debye shielding by
  ! particles with vthermal < v, it contains a smooth transition which I
  ! have rigorously shown from the Balescu-Lenard operator is correct in
  ! the limits vthermal>>v and vthermal<<v.  However, the Balescu-Lenard
  ! operator predicts a substantial amount of Cerenkov radiation of
  ! plasma waves if the particle speed exceeds the electron thermal
  ! speed.  I am ignoring this Cerenkov radiation (and its associated
  ! drag).  To look at it properly, Nat Fisch claims that one would have
  ! to consider the reabsorption of the Cerenkov radiation as well.
  !
  ! Warning:  GWH 3/1/91:  the formulas I am using here provide very
  ! little Debye shielding if the test particles are much faster than the
  ! thermal electrons.  In that case, one would need to worry about whether or
  ! not this is correct.  There appear to be fundamental problems with
  ! the Balescu-Lenard operator.
  !
  ! The main references for this stuff are D.V. Sivukhin, Reviews of
  ! Plasma Physics, VOl. 4 (1966) (especially good on quantum corrections,
  ! physical insight), Trubnikov, Reviews of Plasma Physics, Vol. 1
  ! (1963) (simple, complete derivation of Fokker-Planck equation,
  ! complete with justification for the Coulomb Logarithm in terms of
  ! deriving the collision operator from the dielectric response of the
  ! plasma), Rob Goldston's Ph.D. thesis (references on strong magnetic
  ! field corrections, collisions with multi-electron impurities),
  ! Krommes' class notes on the derivation of the Landau operator from the
  ! Balescu-Lenard operator, and my notes on Krommes' notes relating to
  ! the problem when the test particle is much faster than the thermal
  ! electron velocity.
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, one, two
  implicit none
  real(fp), parameter :: C1000 = 1000.0d0
  !
  integer :: nb,i
  !============
  real(fp) :: sum,omega2,vrel2,rmax,rmincl,rminqu,rmin
  !============
  real(fp) :: et,at,zt,Bmag
  real(fp), dimension(nb) :: denb,tempb,ab,zb,cln

  logical :: llstix       ! tbt 7/10/95

  !
  ! first calculate the maximum impact parameter, rmax:
  !
  ! rmax=( sum_over_species omega2/vrel2 )**(-1/2)
  !
  ! where,
  ! omega2 = omega_p**2 + omega_c**2
  ! vrel2=T/m+2.*E_t/m_t
  !
  sum=ZERO
  do i=1,nb
    omega2=1.74d0*zb(i)**2/ab(i)*denb(i) &
         +9.18d15*zb(i)**2/ab(i)**2*Bmag**2
    vrel2=9.58d10*(tempb(i)/ab(i) + TWO*et/at)
    sum=sum+omega2/vrel2
  end do
  rmax=sqrt(ONE/sum)

  ! next calculate rmin, including quantum corrections.  The classical
  ! rmin is:
  !
  ! rmincl = e_alpha e_beta / (m_ab vrel**2)
  !
  ! where m_ab = m_a m_b / (m_a+m_b) is the reduced mass.
  ! vrel**2 = 3 T_b/m_b + 2 E_a / m_a
  ! (Note:  the two different definitions of vrel2 used in this code
  ! are each correct for their application.)
  !
  ! The quantum rmin is:
  !
  ! rminqu = hbar/( 2 exp(0.5) m_ab vrel)
  !
  ! and the proper rmin is the larger of rmincl and rminqu
  !
  do i=1,nb
    vrel2=9.58d10*(3*tempb(i)/ab(i)+2*et/at)
    rmincl=0.13793d0*abs(zb(i)*zt)*(ab(i)+at)/ab(i)/at/vrel2
    rminqu=1.9121d-8*(ab(i)+at)/ab(i)/at/sqrt(vrel2)
    rmin=max(rmincl,rminqu)
    cln(i)=log(rmax/rmin)
    if(cln(i) .lt. ONE) then
      write(6,*) 'warning from COULOG: coulomb logarithms < 1!'
      cln(i)=ONE
    end if
  end do

  ! debug section:  set all log(lambda)'s to a simple formula for
  ! log(lambda_e) to benchmark with Stix's solutions:

  if(llstix) then
    do i=1,nb
      cln(i)=24.0d0-log(sqrt(denb(1)/1.d6)/(tempb(1)*C1000)) !MG not sure it is correct 
    end do
  end if

  return
end subroutine r8_coulog

