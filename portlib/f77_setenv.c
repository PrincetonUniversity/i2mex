#include "f77name.h"
#include <stdlib.h>

/* f77 callable setenv: DMC Jul 2008, portability not verified */
/*     linux manpages describe setenv as conforming to BSD 4.3 */
/*         if environment variable (var) is already defined */
/*            its value is overwritten! */
/* successful competion returns 0 */

int F77NAME(f77_setenv)(const char* var, const char* val)
{
  return setenv(var,val,1);
}
