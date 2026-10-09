/* cdelete.c */
/* fortran callable routine to delete a file */
/* **UNIX only** */
 
#include <unistd.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(cdelete)(fname)
     char* fname;    /* name of file to delete */
{
  int istat;
 
  istat = unlink(fname);
  if ( istat == 0 ) {
    return 0;
  } else {
    return 1;
  }
 
}
