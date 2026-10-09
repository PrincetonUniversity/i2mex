! freeware integrator adapted for TRANSP: added auxilliary argument
! for integrand (dmc Apr 2003)
!
subroutine tr_r8qags(f,aux,a,b,epsabs,epsrel, &
     reslt,abserr,neval,ier, &
     limit,lenw,last,iwork,work)
  !***begin prologue  dqags
  !***date written   800101   (yymmdd)
  !***revision date  830518   (yymmdd)
  !***category no.  h2a1a1
  !***keywords  automatic integrator, general-purpose,
  !             (end-point) singularities, extrapolation,
  !             globally adaptive
  !***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !           de doncker,elise,appl. math. & prog. div. - k.u.leuven
  !***purpose  the routine calculates an approximation result to a given
  !            definite integral  i = integral of f over (a,b),
  !            hopefully satisfying following claim for accuracy
  !            abs(i-result).le.max(epsabs,epsrel*abs(i)).
  !***description
  !
  !        computation of a definite integral
  !        standard fortran subroutine
  !        double precision version
  !
  !
  !        parameters
  !         on entry
  !            f      - double precision
  !                     function subprogram defining the integrand
  !                     function f(x,aux). the actual name for f needs to be
  !                     declared e x t e r n a l in the driver program.
  !
  !            aux(*) - auxilliary arguments for f, passed thru
  !
  !            a      - double precision
  !                     lower limit of integration
  !
  !            b      - double precision
  !                     upper limit of integration
  !
  !            epsabs - double precision
  !                     absolute accuracy requested
  !            epsrel - double precision
  !                     relative accuracy requested
  !                     if  epsabs.le.0
  !                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !                     the routine will end with ier = 6.
  !
  !         on return
  !            reslt - double precision
  !                     approximation to the integral
  !
  !            abserr - double precision
  !                     estimate of the modulus of the absolute error,
  !                     which should equal or exceed abs(i-reslt)
  !
  !            neval  - integer
  !                     number of integrand evaluations
  !
  !            ier    - integer
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier.gt.0 abnormal termination of the routine
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !            error messages
  !                     ier = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more sub-
  !                             divisions by increasing the value of limit
  !                             (and taking the according dimension
  !                             adjustments into account. however, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties. if
  !                             the position of a local difficulty can be
  !                             determined (e.g. singularity,
  !                             discontinuity within the interval) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling the
  !                             integrator on the subranges. if possible,
  !                             an appropriate special-purpose integrator
  !                             should be used, which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 the occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                             the error may be under-estimated.
  !                         = 3 extremely bad integrand behaviour
  !                             occurs at some points of the integration
  !                             interval.
  !                         = 4 the algorithm does not converge.
  !                             roundoff error is detected in the
  !                             extrapolation table. it is presumed that
  !                             the requested tolerance cannot be
  !                             achieved, and that the returned result is
  !                             the best which can be obtained.
  !                         = 5 the integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of ier.
  !                         = 6 the input is invalid, because
  !                             (epsabs.le.0 and
  !                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
  !                             or limit.lt.1 or lenw.lt.limit*4.
  !                             reslt, abserr, neval, last are set to
  !                             zero.except when limit or lenw is invalid,
  !                             iwork(1), work(limit*2+1) and
  !                             work(limit*3+1) are set to zero, work(1)
  !                             is set to a and work(limit+1) to b.
  !
  !         dimensioning parameters
  !            limit - integer
  !                    dimensioning parameter for iwork
  !                    limit determines the maximum number of subintervals
  !                    in the partition of the given integration interval
  !                    (a,b), limit.ge.1.
  !                    if limit.lt.1, the routine will end with ier = 6.
  !
  !            lenw  - integer
  !                    dimensioning parameter for work
  !                    lenw must be at least limit*4.
  !                    if lenw.lt.limit*4, the routine will end
  !                    with ier = 6.
  !
  !            last  - integer
  !                    on return, last equals the number of subintervals
  !                    produced in the subdivision process, detemines the
  !                    number of significant elements actually in the work
  !                    arrays.
  !
  !         work arrays
  !            iwork - integer
  !                    vector of dimension at least limit, the first k
  !                    elements of which contain pointers
  !                    to the error estimates over the subintervals
  !                    such that work(limit*3+iwork(1)),... ,
  !                    work(limit*3+iwork(k)) form a decreasing
  !                    sequence, with k = last if last.le.(limit/2+2),
  !                    and k = limit+1-last otherwise
  !
  !            work  - double precision
  !                    vector of dimension at least lenw
  !                    on return
  !                    work(1), ..., work(last) contain the left
  !                     end-points of the subintervals in the
  !                     partition of (a,b),
  !                    work(limit+1), ..., work(limit+last) contain
  !                     the right end-points,
  !                    work(limit*2+1), ..., work(limit*2+last) contain
  !                     the integral approximations over the subintervals,
  !                    work(limit*3+1), ..., work(limit*3+last)
  !                     contain the error estimates.
  !
  !***references  (none)
  !***routines called  dqagse,xerror
  !***end prologue  dqags
  !
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: a,abserr,b,epsabs,epsrel,f,reslt,aux(*)
  integer :: ier,last,lenw,limit,lvl,l1,l2,l3,neval
  !
  integer, dimension(limit) :: iwork
  real(fp), dimension(lenw) :: work
  !
  external f
  !
  !         check validity of limit and lenw.
  !
  !***first executable statement  dqags
  ier = 6
  neval = 0
  last = 0
  reslt = 0.0_fp
  abserr = 0.0_fp
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
  !
  !         prepare call for dqagse.
  !
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2
  !
  call tr_r8qagse(f,aux,a,b,epsabs,epsrel,limit,reslt,abserr,neval, &
       ier,work(1),work(l1),work(l2),work(l3),iwork,last)
  !
  !         call error handler if necessary.
  !
  lvl = 0
10 if(ier.eq.6) lvl = 1
  !xx      if(ier.ne.0) call xerror(26habnormal return from dqags,26,ier,lvl)
  !xx      return
  return
end subroutine tr_r8qags


subroutine tr_r8qagse( &
     f,aux,a,b,epsabs,epsrel,limit,reslt,abserr,neval, &
     ier,alist,blist,rlist,elist,iord,last)
  !***begin prologue  dqagse
  !***date written   800101   (yymmdd)
  !***revision date  830518   (yymmdd)
  !***category no.  h2a1a1
  !***keywords  automatic integrator, general-purpose,
  !             (end point) singularities, extrapolation,
  !             globally adaptive
  !***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !***purpose  the routine calculates an approximation result to a given
  !            definite integral i = integral of f over (a,b),
  !            hopefully satisfying following claim for accuracy
  !            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
  !***description
  !
  !        computation of a definite integral
  !        standard fortran subroutine
  !        double precision version
  !
  !        parameters
  !         on entry
  !            f      - double precision
  !                     function subprogram defining the integrand
  !                     function f(x,aux). the actual name for f needs to be
  !                     declared e x t e r n a l in the driver program.
  !
  !            aux(*) - auxilliary arguments for f, passed thru
  !
  !            a      - double precision
  !                     lower limit of integration
  !
  !            b      - double precision
  !                     upper limit of integration
  !
  !            epsabs - double precision
  !                     absolute accuracy requested
  !            epsrel - double precision
  !                     relative accuracy requested
  !                     if  epsabs.le.0
  !                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
  !                     the routine will end with ier = 6.
  !
  !            limit  - integer
  !                     gives an upperbound on the number of subintervals
  !                     in the partition of (a,b)
  !
  !         on return
  !            reslt - double precision
  !                     approximation to the integral
  !
  !            abserr - double precision
  !                     estimate of the modulus of the absolute error,
  !                     which should equal or exceed abs(i-reslt)
  !
  !            neval  - integer
  !                     number of integrand evaluations
  !
  !            ier    - integer
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier.gt.0 abnormal termination of the routine
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !            error messages
  !                         = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more sub-
  !                             divisions by increasing the value of limit
  !                             (and taking the according dimension
  !                             adjustments into account). however, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties. if
  !                             the position of a local difficulty can be
  !                             determined (e.g. singularity,
  !                             discontinuity within the interval) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling the
  !                             integrator on the subranges. if possible,
  !                             an appropriate special-purpose integrator
  !                             should be used, which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 the occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                             the error may be under-estimated.
  !                         = 3 extremely bad integrand behaviour
  !                             occurs at some points of the integration
  !                             interval.
  !                         = 4 the algorithm does not converge.
  !                             roundoff error is detected in the
  !                             extrapolation table.
  !                             it is presumed that the requested
  !                             tolerance cannot be achieved, and that the
  !                             returned result is the best which can be
  !                             obtained.
  !                         = 5 the integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of ier.
  !                         = 6 the input is invalid, because
  !                             epsabs.le.0 and
  !                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
  !                             reslt, abserr, neval, last, rlist(1),
  !                             iord(1) and elist(1) are set to zero.
  !                             alist(1) and blist(1) are set to a and b
  !                             respectively.
  !
  !            alist  - double precision
  !                     vector of dimension at least limit, the first
  !                      last  elements of which are the left end points
  !                     of the subintervals in the partition of the
  !                     given integration range (a,b)
  !
  !            blist  - double precision
  !                     vector of dimension at least limit, the first
  !                      last  elements of which are the right end points
  !                     of the subintervals in the partition of the given
  !                     integration range (a,b)
  !
  !            rlist  - double precision
  !                     vector of dimension at least limit, the first
  !                      last  elements of which are the integral
  !                     approximations on the subintervals
  !
  !            elist  - double precision
  !                     vector of dimension at least limit, the first
  !                      last  elements of which are the moduli of the
  !                     absolute error estimates on the subintervals
  !
  !            iord   - integer
  !                     vector of dimension at least limit, the first k
  !                     elements of which are pointers to the
  !                     error estimates over the subintervals,
  !                     such that elist(iord(1)), ..., elist(iord(k))
  !                     form a decreasing sequence, with k = last
  !                     if last.le.(limit/2+2), and k = limit+1-last
  !                     otherwise
  !
  !            last   - integer
  !                     number of subintervals actually produced in the
  !                     subdivision process
  !
  !***references  (none)
  !***routines called  dqelg,dqk21,dqpsrt
  !***end prologue  dqagse
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: a,abseps,abserr,area,area1,area12,area2,a1, &
       a2,b,b1,b2,correc,abs,defabs,defab1,defab2,max, &
       dres,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax, &
       error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,reslt, &
       small,uflow,aux(*)
  integer id,ier,ierro,iroff1,iroff2,iroff3,jupbnd,k,ksgn, &
       ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
  logical extrap,noext
  !
  real(fp), dimension(limit) :: alist, blist, elist, rlist
  real(fp), dimension(3) :: res3la
  real(fp), dimension(52) :: rlist2
  integer, dimension(limit) :: iord
  !
  external f
  !
  !            the dimension of rlist2 is determined by the value of
  !            limexp in subroutine dqelg (rlist2 should be of dimension
  !            (limexp+2) at least).
  !
  !            list of major variables
  !            -----------------------
  !
  !           alist     - list of left end points of all subintervals
  !                       considered up to now
  !           blist     - list of right end points of all subintervals
  !                       considered up to now
  !           rlist(i)  - approximation to the integral over
  !                       (alist(i),blist(i))
  !           rlist2    - array of dimension at least limexp+2 containing
  !                       the part of the epsilon table which is still
  !                       needed for further computations
  !           elist(i)  - error estimate applying to rlist(i)
  !           maxerr    - pointer to the interval with largest error
  !                       estimate
  !           errmax    - elist(maxerr)
  !           erlast    - error on the interval currently subdivided
  !                       (before that subdivision has taken place)
  !           area      - sum of the integrals over the subintervals
  !           errsum    - sum of the errors over the subintervals
  !           errbnd    - requested accuracy max(epsabs,epsrel*
  !                       abs(reslt))
  !           *****1    - variable for the left interval
  !           *****2    - variable for the right interval
  !           last      - index for subdivision
  !           nres      - number of calls to the extrapolation routine
  !           numrl2    - number of elements currently in rlist2. if an
  !                       appropriate approximation to the compounded
  !                       integral has been obtained it is put in
  !                       rlist2(numrl2) after numrl2 has been increased
  !                       by one.
  !           small     - length of the smallest interval considered up
  !                       to now, multiplied by 1.5
  !           erlarg    - sum of the errors over the intervals larger
  !                       than the smallest interval considered up to now
  !           extrap    - logical variable denoting that the routine is
  !                       attempting to perform extrapolation i.e. before
  !                       subdividing the smallest interval we try to
  !                       decrease the value of erlarg.
  !           noext     - logical variable denoting that extrapolation
  !                       is no longer allowed (true value)
  !
  !            machine dependent constants
  !            ---------------------------
  !
  !           epmach is the largest relative spacing.
  !           uflow is the smallest positive magnitude.
  !           oflow is the largest positive magnitude.
  !
  !***first executable statement  dqagse
  epmach = epsilon(0.0d0) !d1mach_tr(4)
  !
  !            test on validity of parameters
  !            ------------------------------
  ier = 0
  neval = 0
  last = 0
  reslt = 0.0_fp
  abserr = 0.0_fp
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0_fp
  elist(1) = 0.0_fp
  if(epsabs.le.0.0_fp.and. &
       epsrel.lt.max(0.5E+02_fp*epmach,0.5E-28_fp)) &
       ier = 6
  if(ier.eq.6) go to 999
  !
  !           first approximation to the integral
  !           -----------------------------------
  !
  uflow = tiny(0.0d0) !d1mach_tr(1)
  oflow = huge(0.0d0) !d1mach_tr(2)
  ierro = 0
  call tr_r8qk21(f,aux,a,b,reslt,abserr,defabs,resabs)
  !
  !           test on accuracy.
  !
  dres = abs(reslt)
  errbnd = max(epsabs,epsrel*dres)
  last = 1
  rlist(1) = reslt
  elist(1) = abserr
  iord(1) = 1
  if(abserr.le.1.0E+02_fp*epmach*defabs.and.abserr.gt.errbnd) &
       ier = 2
  if(limit.eq.1) ier = 1
  if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or. &
       abserr.eq.0.0_fp) go to 140
  !
  !           initialization
  !           --------------
  !
  rlist2(1) = reslt
  errmax = abserr
  maxerr = 1
  area = reslt
  errsum = abserr
  abserr = oflow
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  if(dres.ge.(0.1E+01_fp-0.5E+02_fp*epmach)*defabs) ksgn = 1
  !
  !           main do-loop
  !           ------------
  !
  do last = 2,limit
    !
    !           bisect the subinterval with the nrmax-th largest error
    !           estimate.
    !
    a1 = alist(maxerr)
    b1 = 0.5_fp*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call tr_r8qk21(f,aux,a1,b1,area1,error1,resabs,defab1)
    call tr_r8qk21(f,aux,a2,b2,area2,error2,resabs,defab2)
    !
    !           improve previous approximations to integral
    !           and error and test for accuracy.
    !
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2) go to 15
    if(abs(rlist(maxerr)-area12).gt.0.1E-04_fp*abs(area12) &
         .or.erro12.lt.0.99_fp*errmax) go to 10
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
10  if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
15  rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max(epsabs,epsrel*abs(area))
    !
    !           test for roundoff error and eventually set error flag.
    !
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
    !
    !           set error flag in the case that the number of subintervals
    !           equals limit.
    !
    if(last.eq.limit) ier = 1
    !
    !           set error flag in the case of bad integrand behaviour
    !           at a point of the integration range.
    !
    if(max(abs(a1),abs(b2)).le.(0.1E+01_fp+0.1E+03_fp*epmach)* &
         (abs(a2)+0.1E+04_fp*uflow)) ier = 4
    !
    !           append the newly-created intervals to the list.
    !
    if(error2.gt.error1) go to 20
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 30
20  alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
    !
    !           call subroutine dqpsrt to maintain the descending ordering
    !           in the list of error estimates and select the subinterval
    !           with nrmax-th largest error estimate (to be bisected next).
    !
30  call tr_r8qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    ! ***jump out of do-loop
    if(errsum.le.errbnd) go to 115
    ! ***jump out of do-loop
    if(ier.ne.0) go to 100
    if(last.eq.2) go to 80
    if(noext) go to 90
    erlarg = erlarg-erlast
    if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
    if(extrap) go to 40
    !
    !           test whether the interval to be bisected next is the
    !           smallest interval.
    !
    if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
    extrap = .true.
    nrmax = 2
40  if(ierro.eq.3.or.erlarg.le.ertest) go to 60
    !
    !           the smallest interval has the largest error.
    !           before bisecting decrease the sum of the errors over the
    !           larger intervals (erlarg) and perform extrapolation.
    !
    id = nrmax
    jupbnd = last
    if(last.gt.(2+limit/2)) jupbnd = limit+3-last
    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      ! ***jump out of do-loop
      if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
      nrmax = nrmax+1
    end do
    !
    !           perform extrapolation.
    !
60  numrl2 = numrl2+1
    rlist2(numrl2) = area
    call tr_r8qelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1E-02_fp*errsum) ier = 5
    if(abseps.ge.abserr) go to 70
    ktmin = 0
    abserr = abseps
    reslt = reseps
    correc = erlarg
    ertest = max(epsabs,epsrel*abs(reseps))
    ! ***jump out of do-loop
    if(abserr.le.ertest) go to 100
    !
    !           prepare bisection of the smallest interval.
    !
70  if(numrl2.eq.1) noext = .true.
    if(ier.eq.5) go to 100
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5_fp
    erlarg = errsum
    go to 90
80  small = abs(b-a)*0.375_fp
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
90  continue
  end do
  !
  !           set final result and error estimate.
  !           ------------------------------------
  !
100 if(abserr.eq.oflow) go to 115
  if(ier+ierro.eq.0) go to 110
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(reslt.ne.0.0_fp.and.area.ne.0.0_fp) go to 105
  if(abserr.gt.errsum) go to 115
  if(area.eq.0.0_fp) go to 130
  go to 110
105 if(abserr/abs(reslt).gt.errsum/abs(area)) go to 115
  !
  !           test on divergence.
  !
110 if(ksgn.eq.(-1).and.max(abs(reslt),abs(area)).le. &
         defabs*0.1E-01_fp) go to 130
  if(0.1E-01_fp.gt.(reslt/area).or.(reslt/area).gt.0.1E+03_fp &
       .or.errsum.gt.abs(area)) ier = 6
  go to 130
  !
  !           compute global integral sum.
  !
115 reslt = 0.0_fp
  do k = 1,last
    reslt = reslt+rlist(k)
  end do
  abserr = errsum
130 if(ier.gt.2) ier = ier-1
140 neval = 42*last-21
999 return
end subroutine tr_r8qagse


subroutine tr_r8qk21(f,aux,a,b,reslt,abserr,resabs,resasc)
  !***begin prologue  dqk21
  !***date written   800101   (yymmdd)
  !***revision date  830518   (yymmdd)
  !***category no.  h2a1a2
  !***keywords  21-point gauss-kronrod rules
  !***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !           de doncker,elise,appl. math. & progr. div. - k.u.leuven
  !***purpose  to compute i = integral of f over (a,b), with error
  !                           estimate
  !                       j = integral of abs(f) over (a,b)
  !***description
  !
  !           integration rules
  !           standard fortran subroutine
  !           double precision version
  !
  !           parameters
  !            on entry
  !              f      - double precision
  !                       function subprogram defining the integrand
  !                       function f(x,aux). the actual name for f needs to be
  !                       declared e x t e r n a l in the driver program.
  !
  !            aux(*) - auxilliary arguments for f, passed thru
  !
  !              a      - double precision
  !                       lower limit of integration
  !
  !              b      - double precision
  !                       upper limit of integration
  !
  !            on return
  !              reslt - double precision
  !                       approximation to the integral i
  !                       result is computed by applying the 21-point
  !                       kronrod rule (resk) obtained by optimal addition
  !                       of abscissae to the 10-point gauss rule (resg).
  !
  !              abserr - double precision
  !                       estimate of the modulus of the absolute error,
  !                       which should not exceed abs(i-result)
  !
  !              resabs - double precision
  !                       approximation to the integral j
  !
  !              resasc - double precision
  !                       approximation to the integral of abs(f-i/(b-a))
  !                       over (a,b)
  !
  !***references  (none)
  !***end prologue  dqk21
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) a,absc,abserr,b,centr,abs,dhlgth,max,min, &
       epmach,f,fc,fsum,fval1,fval2,hlgth,resabs,resasc, &
       resg,resk,reskh,reslt,uflow,aux(*)
  integer j,jtw,jtwm1
  external f
  !
  real(fp), dimension(10) :: fv1, fv2
  real(fp), dimension(5) :: wg
  real(fp), dimension(11) :: wgk, xgk
  !
  !           the abscissae and weights are given for the interval (-1,1).
  !           because of symmetry only the positive abscissae and their
  !           corresponding weights are given.
  !
  !           xgk    - abscissae of the 21-point kronrod rule
  !                    xgk(2), xgk(4), ...  abscissae of the 10-point
  !                    gauss rule
  !                    xgk(1), xgk(3), ...  abscissae which are optimally
  !                    added to the 10-point gauss rule
  !
  !           wgk    - weights of the 21-point kronrod rule
  !
  !           wg     - weights of the 10-point gauss rule
  !
  !
  ! gauss quadrature weights and kronron quadrature abscissae and weights
  ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
  ! bell labs, nov. 1981.
  !
  data wg  (  1) / 0.066671344308688137593568809893332_fp/
  data wg  (  2) / 0.149451349150580593145776339657697_fp/
  data wg  (  3) / 0.219086362515982043995534934228163_fp/
  data wg  (  4) / 0.269266719309996355091226921569469_fp/
  data wg  (  5) / 0.295524224714752870173892994651338_fp/
  !
  data xgk (  1) / 0.995657163025808080735527280689003_fp/
  data xgk (  2) / 0.973906528517171720077964012084452_fp/
  data xgk (  3) / 0.930157491355708226001207180059508_fp/
  data xgk (  4) / 0.865063366688984510732096688423493_fp/
  data xgk (  5) / 0.780817726586416897063717578345042_fp/
  data xgk (  6) / 0.679409568299024406234327365114874_fp/
  data xgk (  7) / 0.562757134668604683339000099272694_fp/
  data xgk (  8) / 0.433395394129247190799265943165784_fp/
  data xgk (  9) / 0.294392862701460198131126603103866_fp/
  data xgk ( 10) / 0.148874338981631210884826001129720_fp/
  data xgk ( 11) / 0.000000000000000000000000000000000_fp/
  !
  data wgk (  1) / 0.011694638867371874278064396062192_fp/
  data wgk (  2) / 0.032558162307964727478818972459390_fp/
  data wgk (  3) / 0.054755896574351996031381300244580_fp/
  data wgk (  4) / 0.075039674810919952767043140916190_fp/
  data wgk (  5) / 0.093125454583697605535065465083366_fp/
  data wgk (  6) / 0.109387158802297641899210590325805_fp/
  data wgk (  7) / 0.123491976262065851077958109831074_fp/
  data wgk (  8) / 0.134709217311473325928054001771707_fp/
  data wgk (  9) / 0.142775938577060080797094273138717_fp/
  data wgk ( 10) / 0.147739104901338491374841515972068_fp/
  data wgk ( 11) / 0.149445554002916905664936468389821_fp/
  !
  !
  !           list of major variables
  !           -----------------------
  !
  !           centr  - mid point of the interval
  !           hlgth  - half-length of the interval
  !           absc   - abscissa
  !           fval*  - function value
  !           resg   - result of the 10-point gauss formula
  !           resk   - result of the 21-point kronrod formula
  !           reskh  - approximation to the mean value of f over (a,b),
  !                    i.e. to i/(b-a)
  !
  !
  !           machine dependent constants
  !           ---------------------------
  !
  !           epmach is the largest relative spacing.
  !           uflow is the smallest positive magnitude.
  !
  !***first executable statement  dqk21
  epmach = epsilon(0.0d0) !d1mach_tr(4)
  uflow = tiny(0.0d0) !d1mach_tr(1)
  !
  centr = 0.5_fp*(a+b)
  hlgth = 0.5_fp*(b-a)
  dhlgth = abs(hlgth)
  !
  !           compute the 21-point kronrod approximation to
  !           the integral, and estimate the absolute error.
  !
  resg = 0.0_fp
  fc = f(centr,aux)
  resk = wgk(11)*fc
  resabs = abs(resk)
  do j=1,5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc,aux)
    fval2 = f(centr+absc,aux)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do
  do j = 1,5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc,aux)
    fval2 = f(centr+absc,aux)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do
  reskh = resk*0.5_fp
  resasc = wgk(11)*abs(fc-reskh)
  do j=1,10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do
  reslt = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)
  if(resasc.ne.0.0_fp.and.abserr.ne.0.0_fp) &
       abserr = resasc* &
       min(0.1E+01_fp,(0.2E+03_fp*abserr/resasc)**1.5_fp)
  if(resabs.gt.uflow/(0.5E+02_fp*epmach)) abserr = max &
       ((epmach*0.5E+02_fp)*resabs,abserr)
  return
end subroutine tr_r8qk21
