#include <stdlib.h>
#include  <errno.h>
#include <stdio.h>
#include "f77name.h"
/* fortran access to C (unix) library "system" call */
/* fortran must pass a null terminated byte string  */

/* 02/09/2011 CLF: check errno                      */ 
int F77NAME(fsystem_mem)()
{
  float vm;
  float rm;

  get_proc_mem(&vm,&rm);
  printf(" --> fsystem: vmem = %g rmem = %g\n",vm,rm);
  printf("              used memory by system()\n");

  return 0;
}
