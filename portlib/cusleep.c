/* cusleep.c */
/* usleep callable from fortran using an integer argument*/
/*   cusleep(int* imicro) */
/*      imicro = sleep for this number of microseconds */

#include <unistd.h>
#include "f77name.h"
 
void F77NAME(cusleep )(int* iarg) {
    int j;
    j = usleep ((unsigned long)(*iarg)); 
}
