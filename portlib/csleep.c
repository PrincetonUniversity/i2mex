/* csleep.c */
/* sleep callable from fortran */
 
#include <unistd.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
void F77NAME(csleep )(iarg)
     int *iarg;   /* integer form -- ascii character to write */
{
  int j;

  j = sleep ( *iarg );
 
}
