subroutine coulog(et,at,zt,denb,tempb,ab,zb,cln,nb,Bmag, &
     llstix)
  use iso_c_binding, only: fp => c_double
!
! (See the warnings below about some (usually small) uncertainties in
! certain parameter regimes.)
!
! calculates the Coulomb logarithm for a test particle of energy et
! (kev), atomic mass at and atomic charge zt, colliding with NB
! maxwellian background species with densities denb (/m**3),
! temperatures tempb (keV), atomic masses ab, atomic charges zb.  The
! electrons should be included as one of the background species, with
! Z=-1 and A=(1./1836.1).  The list of Coulomb logarithms is returned
! in cln.  The magnetic field should be input as bmag (Tesla).
!
! Note that the test particle loglambda's produced by this subroutine
! are not symmetric, i.e. it is not true that
! loglam_alpha_beta=loglam_beta_alpha.  If such symmetry is desired,
! one should set Et=(3/2)T_t and call this routine twice with species
! alpha and beta switched, and take the average value of the log_lambda's
! that are returned.
!
! Written by Greg Hammett, 16-May-1989.
!
! The exact details of the Coulomb logarithm are probably not worth
! worrying about too much.  It has the form
!
!    log(Lambda) = log(r_max/r_min)
!
! For typical values of log(Lambda) ~ 18, a factor of 2 difference in
! r_max/r_min gives only a 4% correction to log(Lambda).  Since the whole
! assumption of small-angle scattering in the Fokker-Planck approximation
! to the collision operator is only accurate to ~O(1/log(Lambda)), this is
! good enough for most purposes.  Often the most important correction is
! the difference between the classical and quantum r_min, which can give about
! an 18% difference in log(Lambda) at T_e=10 keV and n_e=1.e20/m^3 (using the
! simplified thermal expressions in the NRL formulary).
!
! The recipe I use here is similiar to the
! one given in the NRL plasma formulary.  So it includes the quantum
! mechanics correction which is usually important for collisions with
! electrons.  Rather than completely ignoring the Debye shielding by
! particles with vthermal < v, it contains a smooth transition which I
! have rigorously shown (so I have thought) from the Balescu-Lenard
! to be correct in the limits vthermal>>v and vthermal<<v.
!
! In 1989 I wrote the following, which I now think overstates an issue:
! (However, the Balescu-Lenard
! operator predicts a substantial amount of Cerenkov radiation of
! plasma waves if the particle speed exceeds the electron thermal
! speed.  I am ignoring this Cerenkov radiation (and its associated
! drag).  To look at it properly, Nat Fisch claims that one would have
! to consider the reabsorption of the Cerenkov radiation as well.)
!
! I think it was first in 1991 that I found more accurate ways of
! evaluating the integrals in the Balescu-Lenard operator for energetic
! particles, leading to the comment below:
!
! Warning (usually small): GWH 3/1/91 (updated 4/7/2016):
!
! The formulas I am using here reduce Debye shielding somewhat and thus
! increase the drag rate somewhat (via a slow logarithmic factor) if the
! test particle is much faster than the *thermal electrons* (like an
! electron tail).  This corresponds to a (usually small) enhancement of
! the effective log(Lambda) factor by Cerenkov radiation.  This fast
! particle limit is not discussed much in the literature or in textbooks,
! but I think the formulas here are at least roughly okay (usually good to
! a few percent), and are fairly consistent with what one gets from the
! Balescu-Lenard collision operator for high velocity particles colliding
! with thermal particles.  At one time I thought that the Balescu-Lenard
! operator blew up in this limit, but I later found how to do the
! complicated asymptotics of the integrals in this regime properly, and
! it asymptotically has the form given here (if I did the calculations
! right).  Comparing with the Perkins 1965 article cited below also
! supports the formulas here, for the case of very fast ions (faster
! than thermal electrons) slowing down on thermal electrons.  (Compare
! with Eqs. 26 and 28 in Perkins 1965.)  The formulas here also agree
! with Perkins 1965 on the energy drag rate for fast electrons slowing
! down on thermal electrons, but there are some differences for the
! electron drag rate (the rate of loss of directed momentum, sometimes
! called dynamical friction), and for slowing down of fast ions or
! electrons on thermal ions (but for high velocity ions, drag on thermal
! ions is weak compared to drag on thermal electrons anyway).
!
! The biggest difference between the formulas here and in Perkins might
! be for the dynamical friction drag rate for runaway electrons in a
! cold plasma.  The formula here says
!
!    b_max ~ (v/v_t) b_max_thermal,
!
! while Perkins Eq. 28 appears to have
!
!    b_max ~ (v/v_t)^(1/2) b_max_thermal.
!
! So the Perkins formula would reduce the log(Lambda) here by
! log((v_t/v)^(1/2)).  For a 100 keV electron colliding with 10 eV
! thermal electrons, that would reduce ln(Lambda) by about 2.3, or
! about 20% at ln(Lambda) ~ 10.
!
! End of 2016 comment/warning.

!
! The main references for this stuff are D.V. Sivukhin, Reviews of
! Plasma Physics, Vol. 4 (1966) (especially good on quantum corrections,
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
	real et,at,zt,Bmag
	real denb(nb),tempb(nb),ab(nb),zb(nb),cln(nb)

	Logical llstix       ! tbt 7/10/95

!
! first calculate the maximum impact parameter, rmax:
!
! rmax=( sum_over_species omega2/vrel2 )**(-1/2)
!
! where,
! omega2 = omega_p**2 + omega_c**2
! vrel2=T/m+2.*E_t/m_t
!
 sum=0.0
 do i=1,nb
    omega2=1.74*zb(i)**2/ab(i)*denb(i)&
         +9.18e15*zb(i)**2/ab(i)**2*Bmag**2
    vrel2=9.58e10*(tempb(i)/ab(i) + 2.*et/at)
    sum=sum+omega2/vrel2
 end do
 rmax=sqrt(1./sum)
 
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
    vrel2=9.58e10_fp*(3.0_fp*tempb(i)/ab(i)+2.0_fp*et/at)
    rmincl=0.13793_fp*abs(zb(i)*zt)*(ab(i)+at)/ab(i)/at/vrel2
    rminqu=1.9121e-8_fp*(ab(i)+at)/ab(i)/at/sqrt(vrel2)
    rmin=max(rmincl,rminqu)
    cln(i)=log(rmax/rmin)
    if(cln(i) .lt. 1.0_fp) then
       write(6,*) 'warning from COULOG: coulomb logarithms < 1!'
       cln(i)=1.0_fp
    end if
 enddo

! Another 2016 comment:
!
! Comparing with Perkins, 1965,  Phys. Fluids 8, 1361, for v>>v_te,
! and focussing on just the log(Lambda) for the fast particle energy
! loss rate on thermal electrons (the g_1(v) term in his Eq. 26)
! suggests that the above should be modified from:
!
!    ln(Lambda) = log(r_max/r_min)
!
! to
!
!    ln(Lambda) = ln(exp(0.5)/1.467*r_max/r_min)
!               = ln(1.12*r_max/r_min)
!
! This is only a very minor correction!!! (If I did the calculation right.)
! One could eventually also check other regimes and thus the Pade
! interpolation formulas.
!
! Could further modify the above formula to be
!
!    ln(Lambda) => ln(sqrt(1+(1.12*r_max/r_min)**2))
!
! This can be derived from the calculation for the dynamical friction
! drag rate, and gives a positive value for ln(Lambda) for all possible
! r_max/r_min, but it is not correct compared to molecular dynamics
! simulations for small r_max/r_min, where certain correlation effects
! become important.
!
! Expanding for large Lambda, we get
!
!    ln(sqrt(1+(C*Lambda)^2)) = ln(C*Lambda)*ln(sqrt(1+1/(C*Lambda)^2))
!                             ~ ln(C*Lambda) + 0.5/(C*Lambda)^2
!
! and it is known that there are other corrections of order
! 1/sqrt(Lambda_eff) (see Perkins 1965, if I understand him right).
! Perkins calculates the first order correction C, but no one (I think)
! has tried to go beyond that.  These other terms would dominate over
! the above 0.5/Lambda_eff^2 term.  Also, although one could in
! principle calculate the dynamical drag rate or the mean energy loss
! rate to arbitrary order in 1/Lambda, the whole assumption that
! small-angle scattering dominates breaks down so that a Fokker-Planck
! form for the collision operator is no longer valid and one should
! use a Boltzmann collision operator and include large-angle
! scattering.
!
! See papers by Baalrud et al. for a better transition to the strongly coupled
! regime.  An improved fit is if r_max is limited to be no smaller than
! 1/n^(1/3), the average interparticle spacing.  But I think even that does
! not do as well as empirical fits to molecular dynamics (MD) simulations of
! the strongly-coupled regime.

! debug section:  set all log(lambda)'s to a simple formula for
! log(lambda_e) to benchmark with Stix's solutions:

 if(llstix) then
    do i=1,nb
       cln(i)=24.0_fp-log(sqrt(denb(1)/1.e6_fp)/(tempb(1)*1000._fp))
    end do
 end if

 return
end subroutine coulog
