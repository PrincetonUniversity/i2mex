/*  dmc mod -- name must vary depending on machine type & f77 / c linkage */
/*  for portability, use cpp lines */
 
#include "f77name.h"
 
/*  Function WHATBIN
 
      Determine parameters required for linear interpolation
      Find which "bin" an interpolant falls in and return index
      of RHS of that bin
	
      Input:
	  ini	=	1 to begin search from last result
			0 to begin search from lhs of the grid
	  xi(*) =       x grid
	  nxi	=	dimension of x grid
	  x     =       interpolant
	
      Output:
	  fnc	=	index of rhs of enclosing bin
	                0 if out of range
 
      Side Effects:
	  ini	=       1 ; initialization has been established
	
      Use:
	  y = y(j-1) + ( y(j)-y(j-1) ) / ( xi(j)-xi(j-1) ) * dx
	  where dx = x - xi(j-1)
          where j = WHATBIN
 
      Notes:
          Index returned is in FORTRAN notation
 
      By: Dick Wieland
          Jan, 1996
*/
 
int  F77NAME(whatbin) (int *ini, float *xi, int *nxi, float *x)
{
  int iz2 ;
  static int i ;
 
  if (*ini==0) {
    i=1;          /* reset to RHS of 1st bin */
    *ini=1; }
 
  iz2=-1;
 
/* Deal with Special Cases */
 
/* Out of Range */
  if ((*x < xi[i-1]) || (*x > xi[*nxi-1])) return (iz2+1) ;
 
/* At 1st Point on LHS */
  if (*x == xi[0]) {
    iz2=1;
    return (iz2+1); }  /* Fortran Index */
 
/* Find the RHS Index for the bin */
  while (1) {
    if (xi[i] >= *x) {
      iz2=i;
      return (iz2+1); }  /* Fortran Index */
    else {
      i++;
    }
  }
}
/*  dmc mod -- name must vary depending on machine type & f77 / c linkage */
/*  for portability, use cpp lines */
 
#include "f77name.h"
 
/*  Function WHATBIN
 
      Determine parameters required for linear interpolation
      Find which "bin" an interpolant falls in and return index
      of RHS of that bin
	
      Input:
	  ini	=	1 to begin search from last result
			0 to begin search from lhs of the grid
	  xi(*) =       x grid
	  nxi	=	dimension of x grid
	  x     =       interpolant
	
      Output:
	  fnc	=	index of rhs of enclosing bin
	                0 if out of range
 
      Side Effects:
	  ini	=       1 ; initialization has been established
	
      Use:
	  y = y(j-1) + ( y(j)-y(j-1) ) / ( xi(j)-xi(j-1) ) * dx
	  where dx = x - xi(j-1)
          where j = WHATBIN
 
      Notes:
          Index returned is in FORTRAN notation
 
      By: Dick Wieland
          Jan, 1996
*/
 
int  F77NAME(whatbin_r8) (int *ini, double *xi, int *nxi, double *x)
{
  int iz2 ;
  static int i ;
 
  if (*ini==0) {
    i=1;          /* reset to RHS of 1st bin */
    *ini=1; }
 
  iz2=-1;
 
/* Deal with Special Cases */
 
/* Out of Range */
  if ((*x < xi[i-1]) || (*x > xi[*nxi-1])) return (iz2+1) ;
 
/* At 1st Point on LHS */
  if (*x == xi[0]) {
    iz2=1;
    return (iz2+1); }  /* Fortran Index */
 
/* Find the RHS Index for the bin */
  while (1) {
    if (xi[i] >= *x) {
      iz2=i;
      return (iz2+1); }  /* Fortran Index */
    else {
      i++;
    }
  }
}
