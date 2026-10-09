#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* c_trscalar   -- wrapper for trprofil
                08/09/01 C.Ludescher
   $Id: c_trprofil.c,v 1.7 2009-01-30 16:19:58 Rob_Andre Exp $
*/
void c_trprofil(char *disk, char *dir, char *runid , char *dname,
                int *maxtimes, int *maxdata,
                char *label, char *units, int *itype , int *numx,
                int *numtimes, float *times, float *fdat0,
	        int *ierr)
{
#include "f77name.h"
  /* assume "disk" contains an environment e.g. $RESULTDIR
     needs to be translated and passed as "full path" to trscalar
  */
 
  int ldev, ldir, lr, lnam, ll, lu;
  char *p;
  char path[100];
  char dummy[2];
 
     ll=strlen(label);
     lu=strlen(units);
     ldev=strlen(disk);
     ldir=strlen(dir);
     lr=strlen(runid);
     lnam=strlen(dname);
 
     if (disk[0] == '$') {
       if (p = (char *)getenv(disk+1)) strcpy(path,p);
       strcat(path,dir);
       ldir=strlen(path);
       if (!strcmp(disk,"$ARCDIR")){
	 if (path[ldir-3] == '.') path[ldir-3] = '/';
       }
       dummy[0] = ' ';
       ldev = 1;
       F77NAME(trprofil)(dummy, path, runid, dname,
		 maxtimes, maxdata,
		 label, units, itype, numx,
		 numtimes, times, fdat0,
		 ierr, ldev, ldir, lr, lnam, ll, lu);
     }
     else {
       F77NAME(trprofil)(disk, dir, runid, dname,
		 maxtimes, maxdata,
		 label, units, itype, numx,
		 numtimes, times, fdat0,
		 ierr, ldev, ldir, lr, lnam, ll, lu);
     }
}


void c_trprofil_connect(char *disk, char *dir, char *runid , char *dname,
			int *idrun, char *label, char *units, int *ifcn, int *itype , int *numx,
			int *numtimes, int *imds, int *ierr)
{
#include "f77name.h"
  /* assume "disk" contains an environment e.g. $RESULTDIR
     needs to be translated and passed as "full path" to trscalar
  */
  
  int ldir ;
  char *p;
  char path[100];
  char dummy[2];
  /* fprintf(stderr,"inside c_trprofil_connect\n") ; */
  if (disk[0] == '$') {
    if (p = (char *)getenv(disk+1)) strcpy(path,p);
    strcat(path,dir);
    ldir=strlen(path);
    if (!strcmp(disk,"$ARCDIR")){
      if (path[ldir-3] == '.') path[ldir-3] = '/';
    }
    strncpy(dummy, " ",2) ;
    F77NAME(trprofil_connect)(dummy, path, runid, dname,
			      idrun, label, units, ifcn, itype, numx,
			      numtimes, imds, ierr);
  }
  else {
    F77NAME(trprofil_connect)(disk, dir, runid, dname,
			      idrun, label, units, ifcn, itype, numx,
			      numtimes, imds, ierr);
  }
}


void c_trprofil_fetch(int *idrun, char *dname, int *ifcn, int *imds,
		      int *maxtimes, int *maxdata, float *times, float *fdat0,
		      int *ierr)
{
#include "f77name.h"
  int lnam ;

  lnam=strlen(dname);

  if (0) {
    F77NAME(trprofil_fetch)(idrun, dname, ifcn, imds,
			    maxtimes, maxdata, times, fdat0, ierr, 
			    lnam);
  }
  else {
    F77NAME(trprofil_fetch_cstring)(idrun, dname, ifcn, imds,
				    maxtimes, maxdata, times, fdat0, ierr);
  }
}
