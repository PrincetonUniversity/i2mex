#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>

#include "c_logmod.h"

/*
  Mods
     25Jan2010        jim.conboy@ccfe.ac.uk
                      preset c_logmsg_level  to -1, else if an infolog call preceeds openfile,
                      c_genlog is not called, no file open is attempted, & any value in $LOG_LEVEL
                      is ignored..
                      The df value of 1 is restored in c_openlog
 */

static FILE*  c_logmod_file = NULL ;   /* log file */ 
static char   c_logmod_name[256] ;     /* file name */
static int    c_logmod_level = -1 ;    /* log level, <0 for uninitialized*/
static int    c_logmod_isbuf = 0 ;     /* nonzero if buffered */
static pid_t  c_logmod_pid ;           /* process id, set in openlog */
static char   c_logmod_time[64] ;      /* formatted time */
static int    c_logmod_istime = 1 ;    /* true for time stamped messages */
static int    c_logmod_cpu=-1 ;        /* cpu index */
static char   c_logmod_file_tag[16] ;  /* cpu tag for file name */
static char   c_logmod_msg_tag[16] ;   /* cpu tag for message */
static size_t c_logmod_depth ;         /* number of enterLogs - exitLogs */

#define C_LOGMOD_MAXDEPTH 7   /* maximum depth of enterLog, exitLog pairs */
#define C_LOGMOD_OFFDEPTH 2   /* number of spaces to offset for each depth */

#define C_LOGMOD_BUFFDEPTH C_LOGMOD_MAXDEPTH*C_LOGMOD_OFFDEPTH+1 

static void c_logtime() {
  time_t         now ;
  time_t*        pnow ;
  struct timeval tday ;
  struct tm*     snow ;
  char           cbuf[64] ;
  int            frac ;

  if (gettimeofday(&tday,NULL)==0) {
    pnow = &(tday.tv_sec) ;
    frac = tday.tv_usec/10000 ;
  }
  else {
    pnow = &now ;
    time(pnow) ;
    frac = 0 ;
  }
  
  snow = localtime(pnow) ;

  strftime(cbuf, 63, "%d%b%Y %H:%M:%S.%%02d ", snow) ;
  snprintf(c_logmod_time, 63, cbuf, frac) ;
}

/*
  Set the c_logmod_level log level from the environment variable LOG_LEVEL.
     LOG_LEVEL = 0 or 'INFO'
               = 1 or 'WARN'  (default)
               = 2 or 'ERROR'
               = 3 or 'NOMSG'   no log messages
*/
void c_setup_loglevel_env() {
  char* cp ;    
  char* ep ;
  int ic ;  
  
  cp = getenv("LOG_LEVEL") ;

  if (cp!=NULL) {
    ic = strtol(cp, &ep, 0) ;                   /* try to convert to an integer log level */
    if (ep!=cp) {
      F77NAME(c_setloglevel)(&ic) ;
    }
    else {                                      /* try to identify log level by name */
      ep=cp ;
      while ((*ep)==0x20 && (ep-cp)<10) ep++ ;  /* remove leading spaces */
      
      if (strncasecmp(ep,"inf",3)==0) {
	ic=0 ;
      }
      else if (strncasecmp(ep,"err",3)==0) {
	ic=2 ;
      }
      else if (strncasecmp(ep,"nom",3)==0) {
	ic=3 ;
      }
      else {
	ic=1 ;
      }
      F77NAME(c_setloglevel)(&ic) ;
    }
  }
  else{
    ic = 1 ;
    F77NAME(c_setloglevel)(&ic) ;
  }
}

/*
  Open a log file

   file       = file name
   mode       = open mode  
                 '1' -> do nothing if already open otherwise like 'REPLACE'
                 'n' -> new file created, error if file already exists
                 'a' -> reopen an existing file and append new messages
                 'r' -> if old file exists, delete and open a new file
                 'm' -> if old file exists, rename old file as <file.name>~  and open a new file 
   iostat     = return code, nonzero on error
   comm_world = not currently used
   icpu       = cpu index starting from 0 or -1 to not use cpu number
   numcpu     = maximum number of cpus
   gnv        = nonzero to fetch log level from LOG_LEVEL environment variable
*/
void F77NAME(c_openlog)(char* file, char* mode, int* iostat, int* comm_world, int* icpu, int* numcpu, int* gnv) {
  char s ;              /* single character open mode */
  int ier, iexist ;     /* nonzero on error, nonzero if file exists */
  int ic, mc ;          /* cpu index starting from 0 and maximum number of cpus */
  int kc, kn ;          /* cpu index starting from 1, number of digits to write out */
  char cfmt[16] ;       /* formatting string */
  char* cp ;    
  char* ep ;
  char file_new[257] ;  /* allocate len(c_logmod_name)+1 */
  size_t ilen,itlen,ix,iy,it ; 

  int F77NAME(c_loglevel)();

  *iostat = 0 ;
  s = mode==NULL ? '1' : *mode ;

  c_logmod_depth = 0 ;

  /* --- set logging level --- */
  if (gnv==NULL || (*gnv)!=0) {
    c_setup_loglevel_env() ;
  }

  /* --- close existing log file --- */
  if (c_logmod_file!=NULL) {
    if (s=='1') return ;

    /* write out trailer message */
    if (c_logmod_istime) {
      c_logtime() ;
    }
    else {
      c_logmod_time[0]=0 ;
    }
    fprintf(c_logmod_file, "%s%s-closeLog: PID=%d info_level=%d\n", 
	    c_logmod_time, c_logmod_msg_tag, c_logmod_pid, F77NAME(c_loglevel)()) ;

    ier = fclose(c_logmod_file) ;
    c_logmod_file = NULL ;

    if (ier!=0) {
      fprintf(stderr, "?c_openlog: error closing open log file\n") ;
      *iostat=1 ; return ;
    }
  }
  if (s=='1') s='r' ;

  /* --- check for existing file --- */
  c_logmod_pid = getpid() ;

  if (file==NULL || strlen(file)==0) {
    sprintf(c_logmod_name, "portlib_%d.log", c_logmod_pid) ;  /* no file name so use "portlib_<pid>.log" */
  }
  else if (strlen(file)>255) {
    fprintf(stderr, "?c_openlog: file name is too long (max 255 chars): '%s'\n",file) ;
    *iostat=1 ; return ;    
  }
  else {
    strncpy(c_logmod_name, file, 256) ;
  }
  file = c_logmod_name ;

  ic = icpu!=NULL   ? *icpu   : -1 ;
  mc = numcpu!=NULL ? *numcpu : -1 ;

  if (ic>=0) {
    /* modify file name to accomodate cpu id */
    c_logmod_cpu=ic ;

    kc = ic+1 ;
    mc = kc>mc ? kc : mc ;
    mc = mc-1 ;            /* highest cpu index which will be written */
    kn = 0 ;              
    do {
      kn++ ; 
      mc=mc/10 ;
    } while(mc>0) ;

    snprintf(cfmt,16,"_%%0%dd",kn) ;
    snprintf(c_logmod_file_tag, 16, cfmt, ic) ;
    snprintf(cfmt,16,"[%%0%dd] ",kn) ;
    snprintf(c_logmod_msg_tag, 16, cfmt, ic) ;

    itlen = strlen(c_logmod_file_tag) ;
    ilen  = strlen(c_logmod_name) ;

    if ((itlen+ilen)>255) {
      fprintf(stderr, "?c_openlog: file name with cpu id is too long (max 255 chars): '%s'\n",file) ;
      *iostat=1 ; return ;    
    }

    ix = ilen ;                        /* put _nn tag at c_logmod_name[ix] */
    for (it=0 ; it<ilen ; it++) {
      if (c_logmod_name[ilen-(it+1)]=='.') {
	ix = ilen-(it+1) ;
	break ;
      } ;
    }
    for (iy=ilen-1 ; iy>=ix ; iy--) c_logmod_name[iy+itlen] = c_logmod_name[iy] ;      /* make room for _nn tag */
    for (iy=0 ; iy<itlen ; iy++)    c_logmod_name[iy+ix]    = c_logmod_file_tag[iy] ;  /* put in _nn tag */
    c_logmod_name[itlen+ilen]=0 ;
  }
  else {
    c_logmod_cpu         = -1 ;
    c_logmod_file_tag[0] =  0 ;
    c_logmod_msg_tag[0]  =  0 ;
  }
  
  iexist = access(file, F_OK)==0 ;

  /* --- open the log --- */
  if (s=='r' && iexist) {
    ier = unlink(file) ;    /* remove existing file */
    if (ier!=0) {
      fprintf(stderr, "?c_openlog: unable to delete existing file '%s'\n",file) ;
      *iostat=1 ; return ;
    }
    iexist = 0 ;
  }

  if (s=='m' && iexist) {
    ilen  = strlen(file) ;
    for (it=0 ; it<ilen ; it++) {
      file_new[it] = file[it] ;
    }
    file_new[ilen]   = '~' ;
    file_new[ilen+1] = 0 ;
    ier = rename(file,file_new) ;    /* rename existing file to file~ */
    if (ier!=0) {
      ier = unlink(file_new) ;    /* remove existing file_new first */
      if (ier!=0) {
	fprintf(stderr, "?c_openlog: unable to remove existing file '%s'\n",file_new) ;
	*iostat=1 ; return ;
      }
      ier = rename(file,file_new) ;    /* rename existing file to file~ */
      if (ier!=0) {
	fprintf(stderr, "?c_openlog: unable to rename existing file '%s'\n",file) ;
	*iostat=1 ; return ;
      }
    }
    iexist = 0 ;
  }
  
  if (s=='n'||s=='r'||s=='m') {
    if (iexist) {
      fprintf(stderr, "?c_openlog: can not open file, file alreay exists: '%s'\n",file) ;
      *iostat=1 ; return ;
    }
  
    c_logmod_file = fopen(file,"w") ;

    if (c_logmod_file==NULL) {
      fprintf(stderr, "?c_openlog: error opening new log file: '%s'\n",file) ;
      *iostat=1 ; return ;
    }
  }
  else if (s=='a') {
    if (!iexist) {
      fprintf(stderr, "?c_openlog: can not append to nonexistant file: '%s'\n",file) ;
      *iostat=1 ; return ;
    }

    c_logmod_file = fopen(file,"a") ;

    if (c_logmod_file==NULL) {
      fprintf(stderr, "?c_openlog: error opening log file for append: '%s'\n",file) ;
      *iostat=1 ; return ;
    }
  }
  else {
    fprintf(stderr, "?c_openlog: unknown file mode option = %c\n",s) ;
    *iostat=1 ; return ;
  }

  /* write out header message */
  if (c_logmod_istime) {
    c_logtime() ;
  }
  else {
    c_logmod_time[0]=0 ;
  }
  fprintf(c_logmod_file, "%s%s-openLog: PID=%d info_level=%d\n", 
	  c_logmod_time, c_logmod_msg_tag, c_logmod_pid, F77NAME(c_loglevel)()) ;
} 

void F77NAME(c_closelog)(int* iostat) {
  int ier ;              /* error flag */

  *iostat = 0 ;

  c_logmod_depth = 0 ;

  int F77NAME(c_loglevel)();

  if (c_logmod_file==NULL) {
    fprintf(stderr, "?c_closelog: a log file is not open\n") ;
    *iostat=1 ; return ;
  }

  /* write out trailer message */
  if (c_logmod_istime) {
    c_logtime() ;
  }
  else {
    c_logmod_time[0]=0 ;
  }
  fprintf(c_logmod_file, "%s%s-closeLog: PID=%d info_level=%d\n", 
	  c_logmod_time, c_logmod_msg_tag, c_logmod_pid, F77NAME(c_loglevel)()) ;

  /* close log */
  ier = fclose(c_logmod_file) ;
  c_logmod_file = NULL ;

  if (ier!=0) {
    fprintf(stderr, "?c_closelog: error closing the log file\n") ;
    *iostat=1 ; return ;
  }
}

/*
  Return nonzero if the log file is opened.
*/
int F77NAME(c_isopenlog)() {
  return c_logmod_file==NULL ? 0 : 1 ;
}

void F77NAME(c_getlogname)(int* csize, char* name, int* iret) {
  size_t clast ;   /* position of last character in name */

  if (c_logmod_file==NULL) {
    *iret=-1 ; return ;
  }

  *iret = strlen(c_logmod_name) ;
  if (*iret >= *csize) {                 /* need room for final null */
    *iret=-2 ; return ;
  }
  
  clast = *csize ;
  strncpy(name, c_logmod_name, clast) ;
  name[clast]=0 ;                        /* make sure this is null terminated */
}

void F77NAME(c_setloglevel)(int* ilevel) {
  int i,j ;

  i = *ilevel ;
  j = c_logmod_level ;
  if (c_logmod_file!=NULL) {
    if ( j == 0 ) {
      F77NAME(c_infologi)("c_setloglevel: info_level set to:", ilevel) ;
    } else {
      F77NAME(c_warnlogi)("c_setloglevel: info_level set to:", ilevel) ;
    }
  }
  c_logmod_level = i>3 ? 3 : (i<0 ? 0 : i) ;
} ;

int F77NAME(c_loglevel)() {
  if (c_logmod_level<0) c_setup_loglevel_env() ;
  return c_logmod_level ;
} ;

void F77NAME(c_setlogbuffered)(int* ibuff) {
  c_logmod_isbuf = (*ibuff)!=0 ;
  if (!c_logmod_isbuf && c_logmod_file!=NULL) fflush(c_logmod_file) ;
} ;

void F77NAME(c_setlogtime)(int* itime)  {
  c_logmod_istime = (*itime)!=0 ;
} ;

/*
  Write out a log message at any acceptable log level with optional integer and double values.
  Open a log file if not already open.
     jlevel = level of the message, must be >=c_logmod_level to be written to the log file
              Mod: level evaluated modulo 10, suppress leading identifier when jlevel>=10
     msg    = message
     iv     = optional integer value to write after message
     dv     = optional double  value to write after message
*/
static void F77NAME(c_genlog)(int jlevel, char* msg, int* iv, double* dv) {
  int    ier ;       /* open error flag */
  char   c[2] ;      /* character identifier of level */
  int    ilevel ;    /* actual log level of message */
  size_t ideep ;     /* depth of enterLog,exitLog pairs */

  char   dbuf[C_LOGMOD_BUFFDEPTH] ;
  char*  p ;
  size_t i ;
  int    k ;

  ilevel = jlevel%10 ;

  if (c_logmod_level<0) c_setup_loglevel_env() ; /* log level not set so set from LOG_LEVEL environment variable */

  if (ilevel<c_logmod_level) return ; 

  if (c_logmod_file==NULL) {
    k = 0 ;   /* don't set log level from environment because already set somehow */

    F77NAME(c_openlog)(NULL, NULL, &ier, NULL, NULL, NULL, &k) ;   /* try to open portlib_<pid>.log */
    if (ier!=0) {
      fprintf(stderr, "?c_genlog: log file open error from message logging\n") ;   /* might get a lot of these messages */
      return ;
    } ;
  }

  c[1]=0 ;

  if (jlevel>=10) {
    c[0]=0 ;
  }
  else if (ilevel==0) {
    c[0]='%' ;
  }
  else if (ilevel==1) {
    c[0]='!' ;
  }
  else {
    c[0]='?' ;
  }

  if (c_logmod_istime) {
    c_logtime() ;
  }
  else {
    c_logmod_time[0]=0 ;
  }

  ideep = C_LOGMOD_OFFDEPTH*(c_logmod_depth<=C_LOGMOD_MAXDEPTH ? c_logmod_depth : C_LOGMOD_MAXDEPTH) ;

  p = dbuf ;
  for (i=0 ; i<ideep ; i++) *(p++)=' ' ;
  *p=0 ;

  if (iv!=NULL) {
    fprintf(c_logmod_file, "%s%s%s%s%s %d\n", c_logmod_time, c_logmod_msg_tag, dbuf, c, msg, *iv) ;
  }
  else if (dv!=NULL) {
    fprintf(c_logmod_file, "%s%s%s%s%s %e\n", c_logmod_time, c_logmod_msg_tag, dbuf, c, msg, *dv) ;
  }
  else {
    fprintf(c_logmod_file, "%s%s%s%s%s\n", c_logmod_time, c_logmod_msg_tag, dbuf, c, msg) ;
  }

  if (!c_logmod_isbuf) fflush(c_logmod_file) ;
};

void F77NAME(c_infolog)(char* msg) {
  F77NAME(c_genlog)(0,msg,NULL,NULL) ;
}

void F77NAME(c_infologi)(char* msg, int* ivalue) {
  F77NAME(c_genlog)(0,msg,ivalue,NULL) ;
}

void F77NAME(c_infologd)(char* msg, double* dvalue) {
  F77NAME(c_genlog)(0,msg,NULL,dvalue) ;
}

void F77NAME(c_warnlog)(char* msg) {
  F77NAME(c_genlog)(1,msg,NULL,NULL) ;
}

void F77NAME(c_warnlogi)(char* msg, int* ivalue) {
  F77NAME(c_genlog)(1,msg,ivalue,NULL) ;
}

void F77NAME(c_warnlogd)(char* msg, double* dvalue) {
  F77NAME(c_genlog)(1,msg,NULL,dvalue) ;
}

void F77NAME(c_errorlog)(char* msg) {
  F77NAME(c_genlog)(2,msg,NULL,NULL) ;
}

void F77NAME(c_errorlogi)(char* msg, int* ivalue) {
  F77NAME(c_genlog)(2,msg,ivalue,NULL) ;
}

void F77NAME(c_errorlogd)(char* msg, double* dvalue) {
  F77NAME(c_genlog)(2,msg,NULL,dvalue) ;
}

/*
  writes out a warning level log of the form
  ... \Enter\ sub

  Should be paired with a c_exitlog(sub)

  sub = subroutine name
*/
void F77NAME(c_enterlog)(char* sub) {
  char   cbuf[64] ;
  size_t lsub, lbuf ;

  if (c_logmod_file==NULL) return ;   /* do not call before log file has been opened */

  lsub = strlen(sub) ;
  if (lsub>50) lsub=50 ;
  
  strcpy(cbuf,"\\Enter\\ ") ;
  lbuf = strlen(cbuf);

  strncpy(cbuf+lbuf,sub,lsub) ;

  cbuf[lsub+lbuf]=0 ;
  
  F77NAME(c_genlog)(11,cbuf,NULL,NULL) ;

  if (c_logmod_depth<C_LOGMOD_MAXDEPTH) c_logmod_depth++ ;
}

/*
  writes out a warning level log of the form
  ... /Exit/ sub

  Should be paired with a c_enterlog(sub)

  sub = subroutine name
*/
void F77NAME(c_exitlog)(char* sub) {
  char   cbuf[64] ;
  size_t lsub, lbuf ;

  if (c_logmod_file==NULL) return ;   /* do not call before log file has been opened */

  lsub = strlen(sub) ;
  if (lsub>50) lsub=50 ;
  
  strcpy(cbuf,"/Exit / ") ;
  lbuf = strlen(cbuf);

  strncpy(cbuf+lbuf,sub,lsub) ;

  cbuf[lsub+lbuf]=0 ;

  if (c_logmod_depth>0) c_logmod_depth-- ;

  F77NAME(c_genlog)(11,cbuf,NULL,NULL) ;
}

