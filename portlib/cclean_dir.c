/* Remove file contents of a directory -- make a clean directory */
/* NOTE -- hidden files not removed; subdirectories not touched. */

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include "f77name.h"
 
#define FPSIZE 200

int F77NAME(cclean_dir)(cpath)

     const char *cpath;

{
  struct stat fileStat;
  DIR *dp;
  struct dirent *ep;     

  char fpath[FPSIZE];

  int cdelete();
  int iret;

  /* -- start of executable code -- */

  dp = opendir(cpath);

  if (dp == NULL) {
    fprintf(stderr, "cclean_dir: directory not found: %s\n", cpath);
    return 1;
  } else {

    iret = 0;

    while ( ep = readdir(dp) )
      {
	strcpy(fpath,cpath);
	strcat(fpath,ep->d_name);

	if ( ep->d_name[0] == '.' ) continue;

	if(stat(fpath,&fileStat) < 0) {

	  fprintf(stderr, "cclean_dir: stat failed on: %s\n", fpath);
	  iret = 1;

	} else {

	  if( (S_ISDIR(fileStat.st_mode)) ) {
	    fprintf(stderr, "cclean_dir: skipped subdirectory: %s\n", fpath);

	  } else {
            /* not a directory, not ahidden file: delete it */

	    if ( F77NAME(cdelete)(fpath) != 0 ) iret = 1;

	  }
	}
      }
  }
  return iret;
}
