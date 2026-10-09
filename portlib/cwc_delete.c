/* Remove file contents of a directory -- make a clean directory */
/* NOTE -- hidden files not removed; subdirectories not touched. */

#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include "f77name.h"
 
#define FPSIZE 200

int F77NAME(cwc_delete)(cpath,cprefix,csuffix)

     const char *cpath;     /* directory path */
     const char *cprefix;   /* filename prefix (could be empty) */
     const char *csuffix;   /* filename suffix (could be empty) */

{
  struct stat fileStat;
  DIR *dp;
  struct dirent *ep;     

  char fpath[FPSIZE],fnam[FPSIZE];

  int cdelete();
  int iret;

  int ii,iskip;
  size_t ilen_prefix,ilen_suffix,ilen_filename;

  /* -- start of executable code -- */

  ilen_prefix = strlen(cprefix);
  ilen_suffix = strlen(csuffix);

  dp = opendir(cpath);

  if (dp == NULL) {
    fprintf(stderr, "cwc_delete: directory not found: %s\n", cpath);
    return 1;
  } else {

    iret = 0;

    while ( ep = readdir(dp) )
      {
	strcpy(fpath,cpath);
	strcat(fpath,ep->d_name);
	strcpy(fnam,ep->d_name);

	if ( fnam[0] == '.' ) continue;

	if(stat(fpath,&fileStat) < 0) {

	  fprintf(stderr, "cwc_delete: stat failed on: %s\n", fpath);
	  iret = 1;

	} else {

	  if( (S_ISDIR(fileStat.st_mode)) ) {
	    fprintf(stderr, "cwc_delete: skipped subdirectory: %s\n", fpath);

	  } else {
            /* not a directory, not ahidden file: check for deletion */

	    ilen_filename = strlen(fnam);

	    iskip = 0;

	    ii = 0;
	    while ( ii < ilen_prefix )
	      {
		if ( cprefix[ii] != fnam[ii] ) iskip = 1;
		ii++;
	      }

	    ii = 0;
	    while ( ii < ilen_suffix )
	      {
		if ( csuffix[ii] != fnam[ilen_filename - ilen_suffix + ii] )
		  iskip = 1;
		ii++;
	      }

	    if( iskip == 0 ) {
	      /* prefix & suffix match, delete file */
	      if ( F77NAME(cdelete)(fpath) != 0 ) iret = 1;
	    }

	  }
	}
      }
  }
  closedir(dp);
  return iret;
}
