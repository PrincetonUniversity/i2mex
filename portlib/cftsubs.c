/* cftsubs.c */
 
/* d.mccune, 22 Mar 1996 */
 
/* low level i/o routines callable from fortran code */
/*  support open/close, read/write of bare stream files in fixed */
/*  512 byte chunks. */
 
/*  use a fortran-like lun indexed file table, locally stored */
 
#include <errno.h>
#include <stdio.h>
#include <string.h>
 
#include "f77name.h"
 
static int lunflag[140] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };
 
static FILE *lunfile[140];
 
/* function names link differently depending on OS -- see "f77name.h" */
 
#define IO_ERRTEST errno != 0

#ifndef FILE_OPEN_WRITE
#define FILE_OPEN_WRITE fopen(fname, facod)
#endif
 
/* -------------------------------------------------------------- */
/* cftopn  ...  open a byte stream file for a fortran program     */
/*                                                                */
/*    fortray style lun [btw 0 and 100 -- range check by caller]  */
/*    filename and access code are passed                         */
 
void F77NAME(cftopn)(ilun,fname,facod,iquiet,ier)
     int *ilun;   /* logical unit number */
     const char *fname;  /* filename */
     char *facod;  /* access code, 'w' or 'r' */
     int *iquiet; /* =1 for quiet mode; suppress file open error message */
     int *ier;    /* completion code (errno), returned. */
 
{
  FILE *f;
  int i;
 
  i = *ilun;
 
  *ier = 0;
 
  if ( lunflag[i] != 0 )
    {
      fprintf(stderr, " cftopn:  lun %d already in use.\n", *ilun);
      *ier = -1;
      return;
    }
  if ( *fname == 0 )
    {
      fprintf(stderr, " cftopn:  null filename\n");
      *ier = -2;
      return;
    }
 
  errno = 0;
 
  if ( facod[0] == 'w' )
    {
    }
  else
    {
      f = fopen(fname, facod);
    }
 
  if ( f == NULL )
    {
      if ( *iquiet != 1 )
        fprintf(stderr,
                " cftopn(%d,\"%s\",\"%s\") failed:  %s\n",
                *ilun, fname, facod, strerror(errno));
      *ier = errno;
      return;
    }
 
/*  success  */
 
  lunflag[i] = 1;
 
  lunfile[i] = f;
 
  return;
}
 
 
/* -------------------------------------------------------------- */
/* cftcls  ...  close a byte stream file for a fortran program    */
/*                                                                */
/*    fortray style lun [btw 0 and 100 -- range check by caller]  */
/*    is passed.                                                  */
 
void F77NAME(cftcls)(ilun)
     int *ilun;   /* logical unit number */
 
{
  int i;
  int s;
 
  i = *ilun;
 
  if ( lunflag[i] == 0 )
    {
      /* file never opened -- return w/o complaining */
      return;
    }
 
  errno = 0;
  s = 0;
  s = fclose( lunfile[i] );
  if ( s == EOF ) perror(" cftcls:  Close failed.");
 
  lunflag[i] = 0;
 
  return;
}
 
/* -------------------------------------------------------------- */
/* cftrwn  ...  rewind a byte stream file for a fortran program   */
/*                                                                */
/*    fortray style lun [btw 0 and 100 -- range check by caller]  */
/*    is passed.                                                  */
 
void F77NAME(cftrwn)(ilun)
     int *ilun;   /* logical unit number */
 
{
  int i;
  int s;
 
  i = *ilun;
 
  if ( lunflag[i] == 0 )
    {
      fprintf(stderr," cftrwn:  file never opened on lun=%d\n", *ilun);
      return;
    }
 
  rewind( lunfile[i] );
 
  return;
}
 
 
/* -------------------------------------------------------------- */
/* cftcwr  ...  write a byte stream file for a fortran program   */
/*               512 bytes at a time                              */
/*                                                                */
/*    fortray style lun [btw 0 and 100 -- range check by caller]  */
/*    is passed.                                                  */
 
void F77NAME(cftcwr)(ilun,xstream,ier)
     int *ilun;    /* logical unit number */
     char xstream[512];    /* 512 bytes to write  */
     int *ier;     /* completion code */
{
  int i;
  int j;
  int s;
  int c;
 
  i = *ilun;
 
  *ier = 0;
  if ( lunflag[i] == 0 )
    {
      fprintf(stderr," cftcwr:  file never opened on lun=%d\n", *ilun);
      *ier = -1;
      return;
    }
  errno = 0;
  j = 0;
  while ( j < 512 )
    {
      c = xstream[j];
      s = putc(c, lunfile[i] );
      if ( IO_ERRTEST )
        {
          fprintf(stderr," cftcwr:  write error on lun=%d\n", *ilun);
	  fprintf(stderr," strerror:  %s\n", strerror(errno));
          *ier = 1;
          j = 512;
        }
      j++;
    }
 
  return;
}
 
/* -------------------------------------------------------------- */
/* cftcrd  ...  read a byte stream file for a fortran program    */
/*               512 bytes at a time                              */
/*                                                                */
/*    fortray style lun [btw 0 and 100 -- range check by caller]  */
/*    is passed.                                                  */
 
void F77NAME(cftcrd)(ilun,xstream,ier)
     int *ilun;    /* logical unit number */
     char xstream[512];    /* 512 bytes to write  */
     int *ier;     /* completion code */
{
  int i;
  int j;
  int s;
  int c;
 
  i = *ilun;
 
  *ier = 0;
  if ( lunflag[i] == 0 )
    {
      fprintf(stderr," cftcrd:  file never opened on lun=%d\n", *ilun);
      *ier = -1;
      return;
    }
 
  errno = 0;
  j = 0;
  while ( j < 512 )
    {
      c = getc(lunfile[i]);
      xstream[j] = c;
      if ( IO_ERRTEST )
        {
          fprintf(stderr," cftcrd:  read error on lun=%d\n", *ilun);
	  fprintf(stderr," strerror:  %s\n", strerror(errno));
          *ier = 1;
          j = 512;
        }
      j++;
    }
 
  return;
}
