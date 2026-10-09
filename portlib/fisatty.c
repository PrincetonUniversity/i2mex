/* fisatty.c */
 
/* d.mccune, 18 Jun 1996 */
 
/* fortran callable `isatty' */
/* callable as fortran integer function FISATTY(IARG)  */
/*   IARG = 0:  inquiry about stdin                    */
/*        = 1:  inquiry about stdout                   */
/*        = 2:  inquiry about stderr                   */
/*                                                     */
/* return value:                                       */
/*        = 0:  unit attached to file                  */
/*        = 1:  unit attached to tty device            */
 
#include <stdio.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(fisatty)(iarg)
     int *iarg;   /* file descriptor arg, passed by reference */
{
  int i,j;
 
  int isatty();
 
  i = *iarg;
 
  j = isatty(i);
 
  return j;
}
