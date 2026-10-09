/*
  Pair of subroutines for returning the elapsed wall time using the
  gettimeofday() system subroutine.  Usage from fortran

    integer :: istore(8)  ! large enough to hold timeval structure with possible
                          ! future expansion
    real*8  :: dtime      ! elapsed wall time

    call wall_mark_time(istore)           ! ignoring return value
    ... do stuff ...
    call wall_elapsed_time(istore,dtime)

    print *, 'Elapsed wall time = ',dtime
*/

#include <unistd.h>
#include <sys/time.h>
#include <stdio.h>
#include "f77name.h"

/*
  --------------- wall_mark_time ---------------
  Store the current time in an integer array.  The array must be large enough
  to hold a timeval structure.  On a 64 bit RHEL 5 machine this was 16 bytes.
  So to accomodate possible future expansion I would suggest using a longer
  array.  No effort was made to align the timeval data with the integer array.

  Returns 0 if ok, 1 on failure.
*/
int F77NAME(wall_mark_time)(int* itv) {
  struct timeval* ptv ;

  ptv = (struct timeval*) itv ;

  if (gettimeofday(ptv,NULL)!=0) {
    printf("?wall_mark_time: unexpected failure of gettimeofday()") ;
    return 1 ;
  } ;

  return 0 ;
}

/*
  --------------- wall_elapsed_time ---------------
  Return the elapsed wall time in seconds since wall_mark_time() was called with 
  the integer array in the first argument.

  Returns 0 if ok, 1 on failure.
*/
int F77NAME(wall_elapsed_time)(int* itv, double* dtime) {
  struct timeval* pv1 ;
  struct timeval  tv2 ;

  pv1 = (struct timeval*) itv ;

  if (gettimeofday(&tv2,NULL)!=0) {
    printf("?wall_elapsed_time: unexpected failure of gettimeofday()") ;
    return 1 ;
  } ;

  *dtime = (tv2.tv_sec-pv1->tv_sec) + 1.e-6*(tv2.tv_usec-pv1->tv_usec) ;
  return 0 ;
}
