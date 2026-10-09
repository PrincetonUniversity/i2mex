/* crename.c */
/* fortran callable routine to rename a file */
 
#include <stdio.h>
#include <errno.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(crename)(fromname,toname)
     char* fromname; /* name of file to rename */
     char* toname;   /* name to rename it to */
{
  int istat;
 
  istat = rename(fromname,toname);
  if ( istat == 0 ) {
    return 0;
  } else {
    perror("rename failed");
    return 1;
  }
 
}
