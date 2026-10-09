#include "f77name.h"
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <sys/wait.h>
#include <fcntl.h>

extern char **environ ;

#define MAX_ARGS  16
#define MAX_HEAP 512
#define SLEEP_OPEN_COMMAND 5         /* number of seconds between retries of opening the command fifo */
#define OPEN_RETRY 5                 /* maximum number of open retries */
#define MAX_CWD_PATH 181             /* buffer size for holding the cwd sent to the shell server, based on mpi_env_mod.maxlen_val */

int execsystem_extern = 0 ;          /* set nonzero to use these statically defined variables */

char*  execsystem_execute ;          /* point to executable */
char*  execsystem_pargs[MAX_ARGS] ;  /* globally accessible argument pointers */
char   execsystem_heap[MAX_HEAP] ;   /* place to store argument data */
size_t execsystem_nheap ;            /* amount of execsystem_heap used */

pid_t execsystem_pid ;         /* pid returned by fork */

char  execsystem_sent_cwd[MAX_CWD_PATH] ;  /* hold current working directory as sent to shell server */

FILE* exec_fcommand = NULL ;   /* stream to fifo for commands */ 
int   exec_dcommand = -1 ;     /* possible file descriptor for fcommand */
FILE* exec_fresult  = NULL ;   /* stream to fifo for results */

extern FILE* fdopen(int, const char*) ;

/* ============================= c_execsystem_fork ============================= */

/*
  Execute a command with a list of arguments by forking a separate process.

    iexec    = option for running
                0 -> run with arguments in heap memory, menvs>=0 allowed.  It is assumed that the heap memory of the
                     parent can be accessed by the forked child.  Fork is called with execvp() so the path will be 
                     searched for the executable.
                1 -> store arguments in static memory, a full path to the executable is required, menvs<0 required and
                     there will be no verbose messages from the child.  The forked child gets the arguments in its own
                     memory space.  Fork is called with execv() so absolute path must be used.
    kverbose = nonzero for verbose messages to stdout>0 or stderr<0
    execute  = name of executable which will be searched on the path
    margs    = number of arguments in args string
    args     = a concatenated list of arguments each terminated by 0 in the string
    menvs    = number of environment variables or -1 to use the default environment
    envs     = a concatenated list of environment variable definitions each terminated by 0 in the string

  The status code of the child is returned.  Special status codes (might also come
  from the child process),
    40 = fork failed, stderr should contain the reason for the failure
    41 = execve failed in the child, stderr should contain the reason for the failure
    42 = logic error in this subroutine
    43 = the wrong child finished while waiting for forked child
    44 = child abnormally exited for unknown reason (missing test in this subroutine?)
    45 = error changing the environment in the child

    50 = failure using static memory, too many arguments
    51 = failure using static memory, arguments are too long
    52 = failure using static memory, executable path is too long
    53 = failure using static memory, can not modify environment variables
*/
int c_execsystem_fork(const int* iexec, const int* kverbose, const char* execute, 
		      const int* margs, const char* args, 
		      const int* menvs, const char* envs) {

  int   sts, csts, kts ;  /* status for parent and child */
  pid_t child_pid ;       /* child pid */

  size_t i,k ;            /* loop variables */
  int    ii ;
  FILE* vstd ;     /* verbose file stream */
  FILE* cstd ;

  char** pargs ;   /* list of null terminated strings for the arguments */
  char** penvs ;   /* list of null terminated strings for the environment */

  const char* sptr ;     /* modifiable */
  char**      xptr ;
  char*       zptr ;
  char*       zz   ;

  int    nverbose ;      /* int version */
  size_t nargs, nenvs ;  /* int versions */
  int    kenvs ;         /* nonzero for custom envrionment */

  /* ---- init ---- */

  execsystem_extern = (*iexec)!=0 ;

  nverbose = *kverbose ;
  nargs    = (*margs)<0 ? 0 : *margs ;   /* non-negative */
  nenvs    = (*menvs)<0 ? 0 : *menvs ;
  
  kenvs = (*menvs)>=0 ? 1 : 0 ;

  if (nverbose>0) {
    vstd = stdout ;
  }
  else if (nverbose<0) {
    vstd = stderr ;
  }
  else {
    vstd = NULL ;
  }

  /* ---- allocate and convert ---- */
  if (execsystem_extern) {
    if (nargs>=MAX_ARGS) {
      fprintf(stderr, "?c_execsystem: too many arguments = %d, max = %d\n", (int)nargs, MAX_ARGS-1) ;
      return 50 ;
    }

    if (kenvs) {
      fprintf(stderr, "?c_execsystem: a custom environment is not currently supported when using static memory\n") ;
      return 53 ;
    }

    execsystem_nheap = 0 ;

    pargs = execsystem_pargs ;

    k = strlen(execute) ;

    if ((execsystem_nheap+k+1)>MAX_HEAP) {
      fprintf(stderr, "?c_execsystem: executable is too long, max char = %d\n", MAX_HEAP) ;
      return 52 ;
    }
    execsystem_execute =  execsystem_heap+execsystem_nheap ;
    execsystem_nheap  += k+1 ;

    strncpy(execsystem_execute,execute,k+1) ;
  }
  else {
    pargs = (char**)malloc((nargs+1)*sizeof(char*)) ;
  }
  
  pargs[nargs] = NULL ;

  sptr = args ;
  for (i=0 ; i<nargs ; i++) {
    k = strlen(sptr) ;
    
    if (execsystem_extern) {
      if ((execsystem_nheap+k+1)>MAX_HEAP) {
	fprintf(stderr, "?c_execsystem: executable+arguments is too long, max char = %d\n", MAX_HEAP) ;
	return 51 ;
      }
      pargs[i] = execsystem_heap+execsystem_nheap ;
      execsystem_nheap += k+1 ;
    }
    else {
      pargs[i] = (char*)malloc((k+1)*sizeof(char)) ;
    }
    strncpy(pargs[i],sptr,k+1) ;

    sptr += k+1 ;
  }

  if (kenvs) {
    penvs = (char**)malloc((nenvs+1)*sizeof(char*)) ;
    
    penvs[nenvs] = NULL ;

    sptr = envs ;
    for (i=0 ; i<nenvs ; i++) {
      k        = strlen(sptr) ;
      penvs[i] = (char*)malloc((k+1)*sizeof(char)) ;

      strncpy(penvs[i],sptr,k+1) ;

      sptr += k+1 ;
    }
  }
  else {
    penvs = NULL ;    
  }
  
  /* ------ fork ------- */
  if (nverbose!=0) fprintf(vstd,"%%c_execsystem:forking...\n") ;

  execsystem_pid = fork() ;

  switch(execsystem_pid) {
  case -1:
    /* -------- fork error -------- */
    perror("?c_execsystem fork failed") ;  /* write out errno result to stderr */

    sts = 40 ;  /* no child */

    break ;

  case 0:
    /* --------- child ---------- */
    if (execsystem_extern==0) {
      if (nverbose>0) {    /* reset verbose stream */
	cstd = stdout ;
      }
      else if (nverbose<0) {
	cstd = stderr ;
      }
      else {
	cstd = NULL ;

      }

      if (nverbose!=0) {
	fprintf(cstd, ">c_execsystem: child command = \"%s\"\n",execute) ;
	for (i=0 ; i<nargs ; i++) {
	  fprintf(cstd, ">c_execsystem: arg[%d] = \"%s\"\n",(int)i,pargs[i]) ;
	}
      }

      if (kenvs) {
	if (nverbose!=0) fprintf(cstd, ">c_execsystem: changing the environment\n") ;
	
        kts = clearenv() ;
	if (kts) {
	  if (nverbose!=0) fprintf(cstd, ">c_execsystem: error calling clearenv()\n") ;
	  exit(45) ;
	}
	if (nenvs>0 && penvs!=NULL) {
	  xptr = penvs ;
	  while((*xptr)!=NULL) {
	    k    = strlen(*xptr) ;
	    zptr = (char*)malloc((k+1)*sizeof(char)) ;
	    strncpy(zptr,*xptr,k+1) ;
	    zz = strchr(zptr,'=') ;
	    ii = 1 ;
	    if (zz!=NULL && zz!=zptr) {     /* = found and key is not length 0 */
	      *(zz++) = '\0' ;              /* terminate key with 0, point zz to value */
	      if ((*zz)!='\0') {            /* value is not length 0 */
	        kts = setenv(zptr,zz,1) ;
	        if (kts) {
	          if (nverbose!=0) fprintf(cstd, ">c_execsystem: error calling setenv(%s)\n",*xptr) ;
	          exit(45) ;
	        }
	        else {
	          ii=0 ;   /* ok */
	        }
	      }
	    }
	    if (ii!=0 && nverbose!=0) fprintf(cstd, ">c_execsystem: rejected environment '%s'\n",*xptr) ;
	    free(zptr) ;
	    xptr++ ;
	  }
	}
      }

      execvp(execute, pargs) ;                       /* arguments on heap */
    }
    else {
      execv(execsystem_execute, execsystem_pargs) ;  /* arguments in static memory */
    }

    perror("?>c_execsystem child exec failure") ;

    exit(41) ;   /* kill the child */

  default:
    if (nverbose!=0) fprintf(vstd,"%%c_execsystem: I am the parent, started child pid = %d\n",execsystem_pid) ;

    sts = 0 ;
  }

  /* ---------- parent --------- */
  if (execsystem_pid==0) {
    fprintf(stderr, "?>c_execsystem: child unexpectedly in wrong place\n") ;
    exit(42) ;   /* kill wayward child */
  }

  if (sts==0) {
    /* --- wait for child --- */
    child_pid = wait(&csts) ;
    if (child_pid==-1) {
       if (nverbose!=0) fprintf(vstd, "%%c_execsystem: child finished with pid=-1, sleep(1) and retry\n") ;
       sleep(1) ;
       child_pid = wait(&csts) ;
    }
    if (nverbose!=0) fprintf(vstd, "%%c_execsystem: child finished with pid = %d\n",child_pid) ;
  }

  /* --- free -- */
  if (execsystem_extern==0) {
    for (i=0 ; i<nargs ; i++) free(pargs[i]) ;
    free(pargs) ;
  }

  if (kenvs) {
    for (i=0 ; i<nenvs ; i++) free(penvs[i]) ;
    free(penvs) ;
  }

  if (sts!=0) return sts ;

  /* --- handle child results --- */
  if (child_pid == -1) {
    perror("?c_execsystem child returned with pid = -1") ;
    return 43 ;  
  }
  else if (child_pid != execsystem_pid) {
    fprintf (stderr,"?c_execsystem: finished waiting but the wrong child pid returned\n") ;
    return 43 ;
  }

  if (WIFEXITED(csts)) {
    sts = WEXITSTATUS(csts) ;
    if (nverbose!=0) fprintf(vstd,"%%c_execsystem: child exited normally, status = %d\n",sts) ;
  }
  else if (WIFSIGNALED(csts)) {
    sts = WTERMSIG(csts) ;
    if (vstd!=NULL) fprintf(vstd, "?c_execsystem: child exited by uncaught signal, number = %d\n",sts) ;
  }
  else if (WIFSTOPPED(csts)) {
    sts = WSTOPSIG(csts) ;
    if (vstd!=NULL) fprintf(vstd,"?c_execsystem: child exited by stop signal, number = %d\n",sts) ;
  }
  else {
    printf("%%c_execsystem: unknown exit of child\n") ;
    sts = 44 ;
  }

  return sts ;
}



/* ================================= c_connect_shell_server ============================= */

/*
  Connect to a running shell_server through two fifos.

     ioptarg = select option
                 -1 -> close any open fifo streams and return
                  0 -> close any open fifo streams, open new ones and test (usual first call)
                  1 -> check fifo streams are open and test
                  2 -> check fifo streams are open and return
     command = name of the fifo where commands will be sent
     result  = name of the fifo where the results of the command will be fetched
     msg     = string to write out in test or null for default, must not contain ' otherwise
               you will get the default

  This subroutine connects and then sends one echo command to test the shell server.
  Returns 0 on success or nonzero on failure.  exec_fcommand and exec_fresult 
  will be set to the opened fifo streams.

  Notes:
  The goal is to not allow this subroutine to hang since it will be called from rank 0 in TRANSP.  The shell server
  might hang so there should be a mechanism to kill it.

  Usually shell_server will be started first but this should work if this subroutine is called first.  shell_server
  will open command_fifo for reading and block.  This subroutine will repeatedly test (about 25 seconds) opening
  command_fifo for writing.  If it fails, shell_server will be stuck in the block but this subroutine will return
  as a failure.  After command_fifo is connected, this subroutine will try sending a command.  The result_fifo
  was deleted at the beginning and won't be created in shell_server until the test was executed and a return value
  is ready to be sent.  This subroutine will test for the existence of result_fifo and fail if not found.  If found,
  it will block waiting for the connection.  It should be unlikely for a failure in shell_server between creating
  the result_fifo and returning a result.  If this subroutine fails to connect to result_fifo, shell_server will 
  again be blocking. 
*/
int F77NAME(c_connect_shell_server)(int* ioptarg, const char* command, const char* result, const char* msg) {
  int ier ;     /* error flag */
  int istat ;   /* status */
  int ic ;      /* count open retries */
  int iopt ;    /* actual option */

  char* p ;
  const char* pc ;

  char* mycom ; /* trimmed command fifo name */
  char* myres ; /* trimmed result fifo name */
  char  buf[64] ;
  const char* mydef = "hello, fifo" ;

  size_t k ;

  ier  = 0 ;
  iopt = *ioptarg ;

  if (iopt<=0) {
    *execsystem_sent_cwd = 0 ;

    if (exec_fcommand!=NULL) {
      istat = fclose(exec_fcommand) ;
      exec_fcommand = NULL ;
      if (istat!=0) {
	perror("?c_connect_shell_server: error closing open fifo command stream") ;
	ier = 1 ;
      }
    }
    exec_dcommand = -1 ;
    
    if (exec_fresult!=NULL) {
      istat = fclose(exec_fresult) ;
      exec_fresult = NULL ;
      if (istat!=0) {
	perror("?c_connect_shell_server: error closing open fifo result stream") ;
	ier = 1 ;
      }
    }

    if (ier!=0) {
      fprintf(stderr,"?c_connect_shell_server: error closing open fifos\n") ;
      return ier ;
    }
  }

  if (iopt<0) {
    return ier ;
  }
  else if (iopt>0) {
    /* -- streams should already be opened -- */
    if (exec_fcommand==NULL || exec_fresult==NULL) {
      fprintf(stderr,"?c_connect_shell_server: one of the fifo streams has not been opened\n") ;
      ier = 1 ;
    }
    if (ier!=0) *execsystem_sent_cwd = 0 ;

    if (ier!=0 || iopt>=2) return ier ;
  }

  /* --- want to open fifos ---*/
  pc = command ;
  while ((*pc)!=0 && isblank(*pc)) pc++ ;  /* eat up leading whitespace */
  k = strlen(pc) ;
  
  mycom = (char*)malloc((k+1)*sizeof(char)) ;
  strncpy(mycom,pc,k+1) ;

  p = mycom ;
  while((*p)!=0 && isblank(*p)==0) p++ ; *p=0 ;  /* trim */

  pc = result ;
  while ((*pc)!=0 && isblank(*pc)) pc++ ;
  k = strlen(pc) ;
  
  myres = (char*)malloc((k+1)*sizeof(char)) ;
  strncpy(myres,pc,k+1) ;

  p = myres ;
  while((*p)!=0 && isblank(*p)==0) p++ ; *p=0 ;

  /* --- remove myres so as not to block on unused fifo --- */
  if (iopt==0 && ier==0) {
    istat = access(myres, F_OK) ;    /* is there something here */

    if (istat==0) {
      istat = unlink(myres) ;
      if (istat<0) {
	fprintf(stderr,"?c_connect_shell_server: failure to delete fifo or file '%s'\n",myres) ;
	perror("?c_connect_shell_server") ;
	ier=1 ;
      }
    }
  }

  /* --- open command fifo --- */
  if (iopt==0 && ier==0) {
    if (1) {
      /* --- non blocking open, with fixed number of retries --- */
      ic = 0 ;
      while(exec_dcommand<0 && ic<OPEN_RETRY) {
	if (ic>0) sleep(SLEEP_OPEN_COMMAND) ;

	exec_dcommand = open(mycom,O_WRONLY|O_NONBLOCK) ;
	ic++ ;
	if (exec_dcommand<0) {
	  fprintf(stderr,"!c_connect_shell_server: failed to open fifo '%s', sleeping %d seconds, %d more retries\n",
		  mycom,SLEEP_OPEN_COMMAND,(OPEN_RETRY-ic)) ;
	}
      }
      if (exec_dcommand<0) {
	fprintf(stderr, "?c_connect_shell_server: error connecting to command fifo '%s'\n",mycom) ;
	perror("?c_connect_shell_server") ;
	ier=1 ;
      }
      else {
	exec_fcommand = fdopen(exec_dcommand,"w") ;     /* blocking open */
      }
    }
    else {
      exec_fcommand = fopen(mycom,"w") ;
    }

    if (ier!=0 || exec_fcommand==NULL) {
      fprintf(stderr, "?c_connect_shell_server: error connecting to command fifo '%s'\n",mycom) ;
      if (ier==0) perror("?c_connect_shell_server") ;   /* fopen or fdopen failed */
      if (exec_dcommand>0) {
	close(exec_dcommand) ;
	exec_dcommand=-1 ;
      }
      ier=1 ;
    }
  }

  /* --- write test to command and open result fifo --- */
  if (ier==0) {
    pc = (msg!=NULL && strchr(msg,'\'')==NULL) ? msg : mydef ;  /* a ' will confuse shell_server */

    fprintf(exec_fcommand,"&echo echo '%s'\n",pc) ;
    fflush(exec_fcommand) ;

    if (iopt==0) {
      /* -- connect to result fifo -- */
      sleep(2) ;  
      ic = 0 ;
      while (access(myres, F_OK)<0 && ic<OPEN_RETRY) {
	ic++ ;
	fprintf(stderr,"!c_connect_shell_server: can not find fifo '%s', sleeping %d seconds, %d more retries\n",
		myres,SLEEP_OPEN_COMMAND,(OPEN_RETRY-ic)) ;
	sleep(SLEEP_OPEN_COMMAND) ;
      }

      exec_fresult = fopen(myres,"r") ;
      if (exec_fresult==NULL) {
	fprintf(stderr, "?c_connect_shell_server: error connecting to result fifo '%s'\n",myres) ;
	perror("?c_connect_shell_server") ;
	ier=1 ;
      }
    }
  }

  /*  --- read back result of test --- */
  if (ier==0) { 
    p = fgets(buf,sizeof(buf),exec_fresult) ;

    if (p==NULL) {
      fprintf(stderr, "?c_connect_shell_server: error reading test result from fifo '%s'\n",myres) ;

      if (feof(exec_fresult)!=0) {
	fprintf(stderr, "?c_connect_shell_server: EOF flag set on fifo\n") ;
      }
      else if (ferror(exec_fresult)!=0) {
	fprintf(stderr, "?c_connect_shell_server: Error flag set on fifo\n") ;
      }
      ier=1 ;
    }

    if (sscanf(p,"%d",&istat)!=1) {
      fprintf(stderr, "?c_connect_shell_server: error reading fifo result integer\n") ;
      ier=1 ;
    }
    else if (istat!=0) {
      fprintf(stderr, "?c_connect_shell_server: bad result from fifo = %d\n",istat) ;
      ier=1 ;
    }      
  }

  free(mycom) ;
  free(myres) ;

  if (ier!=0) *execsystem_sent_cwd = 0 ;

  return ier ;
}


/* ============================== c_shell_server_send ============================ */

/*
  Send a command to the shell server.  A connection to the server must have already
  been made.
       com = command string

  Returns the status of the command.  If the command expects input, I would expect this to hang.
  Could probably be solved by changing stdin in the forked process to /dev/null.
*/

int F77NAME(c_shell_server_send)(const char* com) {
  int   istat ;
  int   ier ;
  int   idum  = 0 ;
  char* cdum1 = "dummy" ;
  char* cdum2 = "dummy" ;

  size_t k,j ;
  const char* pc ;
  char*  mycom ;
  char*  p ;
  char   buf[64] ;

  ier = 0 ;

  idum  = 2 ;
  istat = F77NAME(c_connect_shell_server)(&idum, cdum1, cdum2, cdum2) ;

  if (istat!=0) {
    fprintf(stderr, "?c_shell_server_send: fifos are not open -- call c_connect_shell_server() first\n") ;
    return 1 ;
  }

  pc = com  ;
  while((*pc)!=0 && isblank(*pc)) pc++ ;  /* remove leading whitespace */
  k = strlen(pc) ;

  if (k==0) return 0 ;  /* empty command */

  mycom = (char*)malloc((k+2)*sizeof(char)) ;  /* room for \n */
  strncpy(mycom,pc,k+1) ;

  p = strchr(mycom,'\n') ;    /* remove existing line feed and terminate string */

  if (p!=NULL) {
    (*p) = 0 ;
    k = strlen(mycom) ;
    if (k==0) {
      free(mycom) ; return 0 ;
    }
  }

  p = mycom+k-1 ;    /* end of string */

  while((p!=mycom) && isblank(*p)) p-- ;   /* find last non blank char */
  *(++p)='\n' ;
  *(++p)=0 ;         /* trim string and terminate with \n */

  /*  -- send command --- */
  k = 1 ;
  j = strlen(mycom) ; ;
  istat = fwrite(mycom, j, k, exec_fcommand) - k ;  /* expect a return value of number of items written */

  if (istat==0) {
    istat = fflush(exec_fcommand) ;
    if (istat!=0) {
      perror("?c_shell_server_send: bad flush") ;
    }
  }

  if (istat!=0) {
    fprintf(stderr, "?c_shell_server_send: error sending command '%s'",mycom) ;
    ier=1 ;
  }
  else {
    /*  --- read back result of test --- */
    p = fgets(buf,sizeof(buf),exec_fresult) ;
    
    if (p==NULL) {
      fprintf(stderr, "?c_shell_server_send: error reading result from fifo\n") ;

      if (feof(exec_fresult)!=0) {
	fprintf(stderr, "?c_shell_server_send: EOF flag set on fifo\n") ;
      }
      else if (ferror(exec_fresult)!=0) {
	fprintf(stderr, "?c_shell_server_send: Error flag set on fifo\n") ;
      }
      ier = 1 ;
    }
    else {
      if (sscanf(p,"%d",&ier)!=1) {
	fprintf(stderr, "?c_shell_server_send: error reading fifo result integer\n") ;
	ier=1 ;
      }
    }
  }

  free(mycom) ;

  return ier ;  
}

/* ============================= c_execsystem ============================= */

/*
  Execute a command with a list of arguments by forking a separate process.

    iexec    = option for running
                0 -> run with arguments in heap memory, menvs>=0 allowed
                1 -> store arguments in static memory, a full path to the executable is required, menvs<0 required and
                     there will be no verbose messages from the child
                2 -> run command to already opened and connected (c_connect_shell_server) shell_server. menvs<0 is
                     required.  If there are three arguments with args[0]="sh" and args[1]="-c", the command will be
                     run as a shell command otherwise the arguments will be packaged in a string for exec. Arguments 
                     with both ' and " not allowed when run as an exec.  The working directory on the server will
                     be changed to the cwd if necessary.
    kverbose = nonzero for verbose messages to stdout>0 or stderr<0
    execute  = name of executable which will be searched on the path
    margs    = number of arguments in args string
    args     = a concatenated list of arguments each terminated by 0 in the string
    menvs    = number of environment variables or -1 to use the default environment
    envs     = a concatenated list of environment variable definitions each terminated by 0 in the string

  The status code of the child is returned.  Special status codes (might also come
  from the child process),
    40 = fork failed, stderr should contain the reason for the failure
    41 = execve failed in the child, stderr should contain the reason for the failure
    42 = logic error in this subroutine
    43 = the wrong child finished while waiting for forked child
    44 = child abnormally exited for unknown reason (missing test in this subroutine?)
    45 = error changing the environment in the child
    46 = cwd buffer is too small, increase MAX_CWD_PATH
    47 = cwd command to shell_server failed

    50 = failure using static memory, too many arguments
    51 = failure using static memory, arguments are too long
    52 = failure using static memory, executable path is too long
    53 = failure using static memory, can not modify environment variables
*/
int F77NAME(c_execsystem)(const int* iexec, const int* kverbose, const char* execute, 
			  const int* margs, const char* args, 
			  const int* menvs, const char* envs) {
  int istat ;    /* status of command */
  int nargs ;    /* number of arguments */

  size_t ktot ;  /* number of chars needed for command */
  size_t k ;     /* length of string */
  size_t i,j ;

  char* com ;    /* command */
  char* p ;
  char* q ;
  char  c ;
  char  buf[MAX_CWD_PATH] ;  /* cwd */

  const char** pargs ; /* pointers to strings in args */
  const char*  pc ;
  FILE* vstd ;         /* verbose file stream */

  istat = 0 ;

  if (*iexec!=2) {
    /* --- fork within process  --- */
    istat = c_execsystem_fork(iexec, kverbose, execute, margs, args, menvs, envs) ;
  }
  else {
    /* --- shell server --- */
    if (*kverbose>0) {
      vstd = stdout ;
    }
    else if (*kverbose<0) {
      vstd = stderr ;
    }
    else {
      vstd = NULL ;
    }
    
    /* -- update cwd on server if necessary -- */
    i = MAX_CWD_PATH-1 ;

    q = getcwd(buf,i) ;
 
    if (q==NULL) return 46 ;    /* cwd is too long for buffer */

    buf[MAX_CWD_PATH-1] = 0 ;
    execsystem_sent_cwd[MAX_CWD_PATH-1] = 0 ;

    ktot = strlen(execsystem_sent_cwd) ;
    k    = strlen(buf) ;

    if (k==ktot) {
      if (strncmp(execsystem_sent_cwd,buf,MAX_CWD_PATH)!=0) *execsystem_sent_cwd = 0 ;
    }
    else {
      *execsystem_sent_cwd = 0 ;
    }

    if (*execsystem_sent_cwd == 0) {
      /* -- update cwd path on server -- */
      q = "&CD '" ;                   /* special command for shell_server */
      j = strlen(q) ;

      com = (char*)malloc((k+j+5)*sizeof(char)) ;  /* len(buf)+j+1 plus slop */

      strcpy(com, q) ;                    
      strncpy(com+j, buf, k) ;
      strcpy(com+j+k, "'") ; 

      if (vstd!=NULL) fprintf(vstd,"%%c_execsystem:%s\n",com) ;

      istat = F77NAME(c_shell_server_send)(com) ;

      free(com) ;

      if (istat!=0) return 47 ;

      strncpy(execsystem_sent_cwd,buf,MAX_CWD_PATH) ;  /* cd was successful */
    }

    /* -- handle command -- */
    if (*menvs>=0) {
      fprintf(stderr, "?c_exec_system: a custom environment is not available when using a shell_server\n") ;
      istat=1 ;
    }
      
    nargs = *margs ;

    if (nargs>0) {
      pargs = (const char**)malloc(nargs*sizeof(char*)) ;

      pc = args ;
      for (i=0 ; i<nargs ; i++) {
	pargs[i] = pc ;
        pc += strlen(pc) + 1 ;
      }
    }

    if (istat==0 && nargs==3 && strcmp(pargs[0],"sh")==0 && strcmp(pargs[1],"-c")==0) {
      /* --- shell command --- */
      if (strchr(pargs[2],'\n')!=NULL) {
	fprintf(stderr, "?c_exec_system: for the shell server, an argument may not contain \\n\n") ;
	istat = 1 ;
      }
      else {
	if (vstd!=NULL) fprintf(vstd,"%%c_execsystem:%s\n",pargs[2]) ;

	istat = F77NAME(c_shell_server_send)(pargs[2]) ;
      }
    }
    else if (istat==0) {
      /* --- exec --- */
      ktot = strlen(execute)+1 ;  /* don't need to count 0 but include anyway just in case */

      for (i=0 ; i<nargs ; i++) {
	if (strchr(pargs[i],'"')!=NULL && strchr(pargs[i],'\'')!=NULL) {
	  fprintf(stderr, "?c_exec_system: for an exec command with the shell server, an argument may not contain both \" and '\n") ;
	  istat = 1 ;
	}
    
	if (strchr(pargs[i],'\n')!=NULL) {
	  fprintf(stderr, "?c_exec_system: for the shell server, an argument may not contain \\n\n") ;
	  istat = 1 ;
	}
 	  
	ktot += strlen(pargs[i])+1 ;
      }
      ktot += 1 + 3*nargs + 2 + 10 ;  /* &, \s"" for each argument, \n,0 and a little bit extra */
      
      if (istat==0) {
	com = (char*)malloc(ktot*sizeof(char)) ;
	p   = com ;
	
	*p++= '&' ;                      /* mark as an exec command */
	k   = strlen(execute) ;
	strncpy(p,execute,k) ; p+=k ;    /* do not copy 0 */
	
	for (i=0 ; i<nargs ; i++) {
	  c = strchr(pargs[i],'"')!=NULL ? '\'' : '"' ;
	  
	  *p++ = ' ' ;                   /* space separated */
	  *p++ = c ;                     /* quoted in command */
	  
	  k = strlen(pargs[i]) ;
	  strncpy(p,pargs[i],k) ; p+=k ;
	  *p++ = c ;                     /* close quote */
	}
	*p++ = '\n' ;
	*p++ = 0 ;
	
	if (vstd!=NULL) fprintf(vstd,"%%c_execsystem:%s\n",com) ;

	istat = F77NAME(c_shell_server_send)(com) ;

	free(com) ;
      }
    }
    if (nargs>0) free(pargs) ;
  }

  return istat ;
}
