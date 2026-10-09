subroutine trx_calc(zname,zdescr,zunits,zexpr,iwarn,ierr)
!
!  invoke the rplot calculator -- return no data
!  the calculator can create named items, to be retrieved by
!  subsequent calls...
!
  implicit NONE
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!  input:
  character*(*) zname  ! name of item to define in computation
  character*(*) zdescr ! description of item
  character*(*) zunits ! physical units
  character*(*) zexpr  ! calculator expression
!
!  output:
  integer iwarn,ierr   ! rpcalc warning and error flags
!
!  local...
!
  integer istype
  real zdum(2)
  integer idum
!
  character*500 zbig_expr
!
!----------------------
!
  zbig_expr = zname//','//zdescr//','//zunits//' = '//zexpr
  call rpcalc(zbig_expr,zdum,0,idum,istype,iwarn,ierr)
!
  return
end subroutine trx_calc
