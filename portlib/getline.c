/*
  $Id: getline.c,v 1.6 2001-08-14 14:15:18 xshare Exp $
 
  NOTE: OSF VERSION NORE _ IN GETLINE_
 
  This simple Fortran-callable interface to the GNU "readline"
  library was written by Joe Freeman (General Atomics, San Diego,
  California, USA) on 16 October 2000, originally for use in the
  General Atomics Fusion Group's "review" data analysis program.
 
  modified dmc 3 11 2000 -- for use in ureadsub, portlib ...
  C strings, not fortran strings, are passed.  see portlib/cstring.for
 
*/
 
#include "f77name.h"
 
#ifdef __GETLINE_EDITOR
 
#include <stdio.h>
#include <stdlib.h>
#include <readline/readline.h>
#include <readline/history.h>
 
int   nchars ;
 
static char *line_read = (char *) NULL ;
extern       HIST_ENTRY **history_list () ;
 
#endif
 
void F77NAME(getline) (char *prompt, char *line, int *linlen) /* strings */
 
{
#ifdef __GETLINE_EDITOR
  *linlen = 0;
  if (line_read)
  {
    free (line_read) ;
    line_read = (char *) NULL ;
  }
  line_read = readline (prompt) ;
  if ( line_read == (char *) 1 ) {
    *linlen = -99;
    return;  /* signals dummy readline routine */
  }
  if ( line_read == (char *) NULL ) {
    *linlen = -1;
    return;  /* signals EOF */
  }
  if (strcmp (line_read, "*history") == 0)
  {
    HIST_ENTRY **list ;
    register int k ;
    list = history_list () ;
    if (list)
    {
      for (k = 0 ; list[k] ; k++)
        printf ("%2d: %s\r\n", k, list[k] -> line) ;
    }
  }
  else if (line_read && *line_read)
  {
    add_history (line_read) ;
  }
  nchars = strlen (line_read) ;
  strcpy (line, line_read) ;
  *linlen = nchars;
#else
  *linlen = -99;  /* signals dummy readline routine, no command editing */
#endif
  return;
}
