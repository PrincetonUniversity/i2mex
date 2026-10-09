/* Linux chdir routine */
 
#include "f77name.h"
 
/* function names link differently depending on OS */
 
int F77NAME(linux_chdir)(name)
char *name;
{
      return chdir(name);
}
 
 
int F77NAME(linux_chdir_)(name)
char *name;
{ 
	return( F77NAME(linux_chdir)(name) );
}
 
 
