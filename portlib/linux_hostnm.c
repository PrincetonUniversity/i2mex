/* Linux hostnm routine */
 
#include <unistd.h>
 
/* function names link differently depending on OS */
 
#include "f77name.h"
 
#if defined __LINUX
int F77NAME(hostnm)(name,len)
#else
int F77NAME(hostnm_dummy)(name,len)
#endif
char *name;
int len;
{
      return gethostname(name,len);
}


#if defined __LINUX
int F77NAME(hostnm_)(name,len)
#else
int F77NAME(hostnm_dummy_)(name,len)
#endif
char *name; int len;
{
	return( gethostname(name,len) );
}
