#include <errno.h>
#include <stdio.h>
#include <string.h>
 
#include "f77name.h"
 
void F77NAME(portlib_isascii)(fname,iansr,ier)
     const char *fname;   /* name of file */
     int *iansr;          /* return 1 if ascii, 0 if binary, -1 if error */
     int *ier;            /* returns 0 if OK, errno value if error on OPEN */
{
 
#define IO_ERRTEST errno != 0
 
#define TESTSIZE 256
 
  FILE *f;
  int  istat;
  int j,c,dtest[TESTSIZE];
 
  int ascii_ct,bin_ct,null_ct;
 
/* executable code */
 
  errno = 0;

  f = fopen(fname, "r");
 
  if ( f == NULL )
    {
      *iansr = -1;
      *ier = errno;
      printf("isascii: file open error: %s\n",fname);
      if ( errno != 0 )
	{
	  perror(fname);
	}
      else
	{
	  *ier=1;
	}
      return;
    }
 
  j = 0;
  while ( j < TESTSIZE )
    {
      dtest[j]=0;
      j++;
    }
 
  errno = 0;
  ascii_ct = 0;
  bin_ct = 0;
  null_ct = 0;
 
  j = 0;
  while ( j < TESTSIZE )
    {
      c = getc( f );
      if ( IO_ERRTEST )
	{
	  j = TESTSIZE;
	}
      else
	{
	  dtest[j] = c;
	  if ( c == 0 ) null_ct++;
	  if ( c >= 32 && c < 127 )
	    {
	      ascii_ct++;
	    }
	  else
	    {
	      if ( c > 0 && c < 7 ) bin_ct++;
	      if ( c > 13 && c < 32 ) bin_ct++;
	      if ( c >= 127 ) bin_ct++;
	    }
	}
      j++;
      /*      printf(" char %d: %d\n", j, c);  */
    }
 
  /*  printf(" %s:\n", fname); */
  /*  printf(" #nulls: %d, #printables: %d, #exceptionals: %d\n",null_ct,ascii_ct,bin_ct) */
  /*  printf("    total %d (< %d OK)\n", null_ct+ascii_ct+bin_c) */
 
  *ier = 0;
  *iansr = 0;
 
  if ( bin_ct == 0 )
    {
      if ( ascii_ct > 0 && null_ct < TESTSIZE/4 ) *iansr=1;
    }
 
  errno = 0;
  istat = fclose ( f );
  return;
}

void F77NAME(portlib_isascii_)(fname,iansr,ier)
     const char *fname;   /* name of file */
     int *iansr;          /* return 1 if ascii, 0 if binary, -1 if error */
     int *ier;            /* returns 0 if OK, errno value if error on OPEN */
{
	F77NAME(portlib_isascii)(fname,iansr,ier);
}
