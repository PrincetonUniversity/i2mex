#include <string.h>
#include <stdlib.h>
/* c_trscalar   -- wrapper for trscalar
                08/09/01 C.Ludescher

   $Id: c_trscalar.c,v 1.7 2009-01-30 19:25:11 Rob_Andre Exp $
*/
void c_trscalar(char *disk, char *dir, char *runid, char *dname,
                int *maxtimes,
                char *label, char *units,
                int *ntimes, float *times, float *sdata,
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
     lr=strlen(runid);
     lnam=strlen(dname);
     ldev=strlen(disk);
     ldir=strlen(dir);
 
     if (disk[0] == '$') {
       if(p =  (char *)getenv(disk+1)) {
	 strcpy(path,p);
       }
       strcat(path,dir);
       ldir=strlen(path);
       if (!strcmp(disk,"$ARCDIR")){
	 if (path[ldir-3] == '.') path[ldir-3] = '/';
       }
       dummy[0] = ' ';
       ldev = 1;
       F77NAME(trscalar)(dummy, path, runid, dname,
		 maxtimes,
		 label, units,
		 ntimes, times, sdata,
		 ierr, ldev, ldir, lr, lnam, ll, lu);
     }
     else {
       F77NAME(trscalar)(disk, dir, runid, dname,
		 maxtimes,
		 label, units,
		 ntimes, times, sdata,
		 ierr, ldev, ldir, lr, lnam, ll, lu);
     }
}


void c_trscalar_connect(char *disk, char *dir, char *runid, char *dname,
			int *idrun, char *label, char *units, int *ifcn,
			int *numtimes, int *imds, int *ierr)
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
  lr=strlen(runid);
  lnam=strlen(dname);
  ldev=strlen(disk);
  ldir=strlen(dir);
 
  if (disk[0] == '$') {
    if(p =  (char *)getenv(disk+1)) {
      strcpy(path,p);
    }
    strcat(path,dir);
    ldir=strlen(path);
    if (!strcmp(disk,"$ARCDIR")){
      if (path[ldir-3] == '.') path[ldir-3] = '/';
    }
    dummy[0] = ' ';
    ldev = 1;
    F77NAME(trscalar_connect)(dummy, path, runid, dname,
			      idrun, label, units, ifcn,
			      numtimes, imds, ierr, ldev, ldir, lr, lnam, ll, lu);
  }
  else {
    F77NAME(trscalar_connect)(disk, dir, runid, dname,
			      idrun, label, units, ifcn,
			      numtimes, imds, ierr, ldev, ldir, lr, lnam, ll, lu);
  }
}



void c_trscalar_fetch(int *idrun, char *dname, int *ifcn, int *imds,
		      int *maxtimes, float *times, float *sdata,
		      int *ierr)
{
#include "f77name.h"
  int lnam ;
  
  lnam=strlen(dname);
  
  if (0) {
    F77NAME(trscalar_fetch)(idrun, dname, ifcn, imds,
			    maxtimes, times, sdata, ierr, 
			    lnam);
  }
  else {
    F77NAME(trscalar_fetch_cstring)(idrun, dname, ifcn, imds,
				    maxtimes, times, sdata, ierr);
  }
}
