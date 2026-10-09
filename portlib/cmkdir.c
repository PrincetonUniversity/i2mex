#include <errno.h>
#include <stdio.h> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include "f77name.h"

/*  attempt to make a directory; return "errno" status code  */

int F77NAME(cmkdir)(const char *path)
{ 
  int istat;

  errno = 0;

  istat = mkdir(path, 0775 );

#ifdef __DEBUG
  fprintf(stderr, " cmkdir %s => istat=%d, errno=%d \n", path, istat, errno);
#endif

  if ( errno == 17 ) errno = 0;  /* ignore: "File exists"... */
  if ( errno != 0 ) {
#ifdef __SUN
    /* fprint strerror(errno) does not work on Solaris */ 
     if ( errno == 30 ) {
         fprintf(stderr, " cmkdir  %s =>  Read Only Filesystem\n", path);
     }else{
         fprintf(stderr, " cmkdir %s => errno=%d \n", path, errno);
     }
#else
    fprintf(stderr, " cmkdir: strerror(%d): %s\n", errno, strerror(errno));
#endif
  }

  return ( errno );
}
