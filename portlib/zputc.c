/* zputc.c */
/* putchar callable from fortran */
 
#include <stdio.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(zputc)(iarg)
     int *iarg;   /* integer form -- ascii character to write */
{
  int j,istat;
 
  j = putchar(*iarg);
  istat = fflush(NULL);
 
  return j;
}
