#define _GNU_SOURCE                              /* avoid whinges about pointers  */
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include "f77name.h"

/*

   fortran access to shell command, returns stdout        
   fortran must pass a null terminated byte string

   Ref http://pubs.opengroup.org/onlinepubs/009695399/functions/popen.html

   Example ( from portlib/mkdir ) 
    call fpopen( 'echo '//cdir//char(0), cdir_ex )             ! expand string

   Author(?)  :  jim.conboy@ccfe.ac.uk 
   Cobbled together from various examples :  25May2011

   Status:
   25May2011     Written & tested on Fedora 10 (Cambridge) ONLY

   Modifications:
   <date>  <user>
           <action>

*/
int F77NAME(fpopen)(const char* cmnd, char* cbuf )
{
  int istat;

  FILE *fp;                                      /* Ptr to pipe                    */
  char line[256];			         /* line of data from unix command */
  const char  cbl = ' ';
  const char cnl = '\n';
  char  *c ;
 
  c = &line[0];
  cbuf[0] = cbl ; cbuf[1] = '\0';                /* Return single blank on error  */

  fp = popen( cmnd, "r");	                 /* Issue the command.	          */

  if( fp ) {
     while ( fgets( line, sizeof line, fp))      /* Read a line			  */
     {
       /* printf("%s",  line);                       Debug !                      */
       cbuf = stpcpy( cbuf, c ); 
     }
     cbuf -= 1 ; *cbuf = '\0' ;                  /* zap trailing newline          */
     pclose(fp);
  }
}

