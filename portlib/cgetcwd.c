/*  f77 callable routine to get cwd for Cray X1 only */

#include "f77name.h"
#include <unistd.h>
#include <errno.h>
 
void cgetcwd_dummy() {}

/* ------------------------------------------------------------------*/
/* C getcwd which returns an integer error code instead of a pointer */

int F77NAME(int_getcwd)(char* buf, int* size) {
    size_t sz ;
    char*  result ;
     
    if (buf==NULL) return 1 ;
    
    sz = (*size)>0 ? *size : 0;
    result = getcwd(buf,sz) ;
    return result==NULL ? errno : 0 ;
}
