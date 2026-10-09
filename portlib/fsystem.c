#include <stdlib.h>
#include  <errno.h>
#include <stdio.h>
#include "f77name.h"
/* fortran access to C (unix) library "system" call */
/* fortran must pass a null terminated byte string  */

/* 02/09/2011 CLF: check errno                      */ 
int F77NAME(fsystem)(const char* cmd)
{
  int istat,errsv;
  float vm;
  float rm;

  /* DMC  get_proc_mem(&vm,&rm); */
  /* DMC  printf(" --> fsystem: vmem = %g rmem = %g\n",vm,rm); */
  /* DMC  printf("              prior to calling system()\n"); */
  errno = 0;
  istat = system(cmd);
  if ( istat != 0 ) {
     printf("fsystem: system(cmd) returns istat = %d\n", istat);
     printf("    cmd: %s\n\n", cmd);
     errsv = errno;
     if ( errsv != 0 ) {
       perror("fsystem->perror:");
       get_proc_mem(&vm,&rm);

       printf(" --> fsystem: vmem = %g rmem = %g\n",vm,rm);
       printf("%s\n",cmd);
       if ( errsv == ENOMEM){
	 printf("fsystem: Insufficient storage space is available\n");
       } else if ( errsv == EAGAIN ){
	 printf("fsystem: Resource unavailable\n");
       } else {
	 printf("fsystem: errno = %d\n", errno);
       }
       printf("  (fsystem failure)\n\n");
       return errsv;
     }
  }
  return istat;
}
