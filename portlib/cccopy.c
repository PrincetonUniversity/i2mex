/* Copy a file as a byte stream -- any file format -- fortran callable */

#include <unistd.h>
#include <stdio.h> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <sys/types.h>
#include <fcntl.h> 
#include "f77name.h"

#define BUFSIZE    256*1024 

/*  
    This C routine copies a file (any format) 
    
    NOTE: the output file permissions are copied from the input file 
    permissions; calling process must have read access of course.

*/
int F77NAME(cccopy)(cinput,coutput)

     const char *cinput;
     const char *coutput;
{ 
  int mystatus, done, numread; 
 
  int imode;
  struct stat fileStat;

  int infd, outfd, iwrite, istat;
  char buf[BUFSIZE]; 
    
  /*   Open the input file */ 
  /*   If input file open fails, exit */

  if(stat(cinput, &fileStat) < 0) {    
    fprintf( stderr, "cccopy: could not stat input file %s.\n", cinput);
    return -1;
  }

  imode = 0;

  if ( fileStat.st_mode & S_IRUSR ) imode = imode + 256;
  if ( fileStat.st_mode & S_IWUSR ) imode = imode + 128;
  if ( fileStat.st_mode & S_IXUSR ) imode = imode + 64;

#if ! defined(__MINGW32__)
  if ( fileStat.st_mode & S_IRGRP ) imode = imode + 32;
  if ( fileStat.st_mode & S_IWGRP ) imode = imode + 16;
  if ( fileStat.st_mode & S_IXGRP ) imode = imode + 8;

  if ( fileStat.st_mode & S_IROTH ) imode = imode + 4;
  if ( fileStat.st_mode & S_IWOTH ) imode = imode + 2;
  if ( fileStat.st_mode & S_IXOTH ) imode = imode + 1;
#endif
  
  if ( (infd = open( cinput, O_RDONLY ) ) == -1 ) {
    fprintf( stderr, "cccopy: input file %s does not exist.\n", cinput);
    mystatus = -1;
  }
  else {
    mystatus = 0;
  }

  if ( mystatus == 0 ) {
    if ( (outfd = open( coutput, O_CREAT|O_WRONLY|O_TRUNC, imode ) ) == -1 ) {
      fprintf( stderr, "cccopy: cannot open output file %s.\n", coutput);
      mystatus = -1;
    }
    else {
      mystatus = 0;
    }
  }

  if ( mystatus < 0 ) return ( -1 );

  /*  OK files opened successfully  */

  mystatus = 0;
  done = 0;

  while ( !done ) {
    numread = read( infd, buf, BUFSIZE );
    if ( numread > 0 ) {
      write( outfd, buf, numread );
    }
    else {
      close( infd );
      close( outfd );
      done = 1;
    }
  }
  return 0;
} 
