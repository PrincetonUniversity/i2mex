/* cpclock.c */
/*   cpu clock -- not for all systems */
 
#include <stdio.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
int F77NAME(cpclock)()
{
  return clock();
}
