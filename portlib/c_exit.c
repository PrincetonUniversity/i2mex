/* c_exit.c */
/*                                     jim.conboy@jet.uk  11Dec2008 
   normal exit callable from fortran 
   avoids segv (after ftn stop ) in do_global_dtors_aux when linking 
   F90 & cpp mixed source  -
   ref http://www.archivum.info/gcc-bugs@gcc.gnu.org/2008-02/msg02037.html  (hint) */
 
#include <stdlib.h>
#include <stdio.h>

/* function names link differently depending on OS */
 
#include "f77name.h"
 
void F77NAME(c_exit )()
{
  printf(" ...portlib/c_exit.c called.\n");
  _Exit(0);
}
void F77NAME(c_exit_)()
{ F77NAME(c_exit )(); }

