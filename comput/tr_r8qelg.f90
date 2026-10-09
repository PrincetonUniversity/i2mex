subroutine tr_r8qelg(n,epstab,reslt,abserr,res3la,nres)
  !***begin prologue  dqelg
  !***refer to  dqagie,dqagoe,dqagpe,dqagse
  !***revision date  830518   (yymmdd)
  !***keywords  epsilon algorithm, convergence acceleration,
  !             extrapolation
  !***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
  !           de doncker,elise,appl. math & progr. div. - k.u.leuven
  !***purpose  the routine determines the limit of a given sequence of
  !            approximations, by means of the epsilon algorithm of
  !            p.wynn. an estimate of the absolute error is also given.
  !            the condensed epsilon table is computed. only those
  !            elements needed for the computation of the next diagonal
  !            are preserved.
  !***description
  !
  !           epsilon algorithm
  !           standard fortran subroutine
  !           double precision version
  !
  !           parameters
  !              n      - integer
  !                       epstab(n) contains the new element in the
  !                       first column of the epsilon table.
  !
  !              epstab - double precision
  !                       vector of dimension 52 containing the elements
  !                       of the two lower diagonals of the triangular
  !                       epsilon table. the elements are numbered
  !                       starting at the right-hand corner of the
  !                       triangle.
  !
  !              reslt - double precision
  !                       resulting approximation to the integral
  !
  !              abserr - double precision
  !                       estimate of the absolute error computed from
  !                       result and the 3 previous results
  !
  !              res3la - double precision
  !                       vector of dimension 3 containing the last 3
  !                       results
  !
  !              nres   - integer
  !                       number of calls to the routine
  !                       (should be zero at first call)
  !
  !***end prologue  dqelg
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: abserr,abs,delta1,delta2,delta3,max,           &
       epmach,epsinf,error,err1,err2,err3,e0,e1,e1abs,e2,e3, &
       oflow,res,reslt,ss,tol1,tol2,tol3
  integer :: i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
  real(fp), dimension(52) :: epstab
  real(fp), dimension(3) :: res3la
  !
  !           list of major variables
  !           -----------------------
  !
  !           e0     - the 4 elements on which the computation of a new
  !           e1       element in the epsilon table is based
  !           e2
  !           e3                 e0
  !                        e3    e1    new
  !                              e2
  !           newelm - number of elements to be computed in the new
  !                    diagonal
  !           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
  !           reslt - the element in the new diagonal with least value
  !                    of error
  !
  !           machine dependent constants
  !           ---------------------------
  !
  !           epmach is the largest relative spacing.
  !           oflow is the largest positive magnitude.
  !           limexp is the maximum number of elements the epsilon
  !           table can contain. if this number is reached, the upper
  !           diagonal of the epsilon table is deleted.
  !
  !***first executable statement  dqelg
  epmach = epsilon(0.0d0) !d1mach_tr(4)
  oflow =  huge(0.0d0) !d1mach_tr(2)
  nres = nres+1
  abserr = oflow
  reslt = epstab(n)
  if(n.lt.3) go to 100
  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = oflow
  num = n
  k1 = n
  do i = 1,newelm
    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs = abs(e1)
    delta2 = e2-e1
    err2 = abs(delta2)
    tol2 = max(abs(e2),e1abs)*epmach
    delta3 = e1-e0
    err3 = abs(delta3)
    tol3 = max(e1abs,abs(e0))*epmach
    if(err2.gt.tol2.or.err3.gt.tol3) go to 10
    !
    !           if e0, e1 and e2 are equal to within machine
    !           accuracy, convergence is assumed.
    !           reslt = e2
    !           abserr = abs(e1-e0)+abs(e2-e1)
    !
    reslt = res
    abserr = err2+err3
    ! ***jump out of do-loop
    go to 100
10  e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 = abs(delta1)
    tol1 = max(e1abs,abs(e3))*epmach
    !
    !           if two elements are very close to each other, omit
    !           a part of the table by adjusting the value of n
    !
    if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
    ss = 0.1E+01_fp/delta1+0.1E+01_fp/delta2-0.1E+01_fp/delta3
    epsinf = abs(ss*e1)
    !
    !           test to detect irregular behaviour in the table, and
    !           eventually omit a part of the table adjusting the value
    !           of n.
    !
    if(epsinf.gt.0.1E-03_fp) go to 30
20  n = i+i-1
    ! ***jump out of do-loop
    go to 50
    !
    !           compute a new element and eventually adjust
    !           the value of result.
    !
30  res = e1+0.1E+01_fp/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2+abs(res-e2)+err3
    if(error.gt.abserr) go to 40
    abserr = error
    reslt = res
40  continue
  end do
  !
  !           shift the table.
  !
50 if(n.eq.limexp) n = 2*(limexp/2)-1
  ib = 1
  if((num/2)*2.eq.num) ib = 2
  ie = newelm+1
  do i=1,ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do
  if(num.eq.n) go to 80
  indx = num-n+1
  do i = 1,n
    epstab(i)= epstab(indx)
    indx = indx+1
  end do
80 if(nres.ge.4) go to 90
  res3la(nres) = reslt
  abserr = oflow
  go to 100
  !
  !           compute error estimate
  !
90 abserr = abs(reslt-res3la(3))+abs(reslt-res3la(2)) &
        +abs(reslt-res3la(1))
  res3la(1) = res3la(2)
  res3la(2) = res3la(3)
  res3la(3) = reslt
100 abserr = max(abserr,0.5E+01_fp*epmach*abs(reslt))
  return
end subroutine tr_r8qelg
