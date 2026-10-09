#ifndef C_LOGMOD_H
#define C_LOGMOD_H

#include "f77name.h"


/*
 Open a log file. A log file will be automatically opened to portlib_<processID>.log
 if a message is written but a log file is not currently opened.  If a cpu number is
 given in icpu, then the file name will be modified with '_nn' added just before the
 last period in the file argument.  As an example if file='foo.log' and icpu=3 and
 numcpu=16 then the file name will be 'foo_03.log'.

    file   = file name (input)
    mode   = select the method of opening the log file (input)
             '1' -> do nothing if already open otherwise like 'r' (default)
             'n' -> new file created, error if file already exists
             'a' -> reopen an existing file
             'r' -> if old file exists, delete and open a new file
    iostat = 0->ok, nonzero->error (output)
    icpu   = cpu number if nonnegative
    numcpu = maximum number of cpus (for formatting)
    gnv    = nonzero to set the logging level from the LOG_LEVEL environment variable (default is nonzero)
*/
extern void F77NAME(c_openlog)(char* file, char* mode, int* iostat, int* comm_world,  int* icpu, int* numcpu,int* gnv) ;

#define openLog F77NAME(c_openlog)


/*
 Close an open log file.
    iostat = 0->ok, nonzero->error (output)
*/
extern void F77NAME(c_closelog)(int* iostat) ;

#define closeLog F77NAME(c_closelog)


/*
  Return the name of the open log file.
      csize = number of available characters in name including terminating null (input)
      name  = log file name is returned (output)
      iret  = number of characters in log file name not including null or <0 if there was an error (output)
              -1 -> no open log file
              -2 -> name[csize] is too small to hold log file name and terminating null
*/
extern void F77NAME(c_getlogname)(int* csize, char* name, int* iret) ;

#define getLogName F77NAME(c_getlogname)


/*
 Set the level of logging.  The levels should be used as,
   ilevel = log level (input)
     0 -> info level, detailed information not to be used for production runs
     1 -> warn level (default), typically used to mark entrance and exit
          of major subroutines or other significant messages
     2 -> error level, used to write an error before aborting the execution
*/
extern void F77NAME(c_setloglevel)(int* ilevel) ;

#define setLogLevel F77NAME(c_setloglevel)


/*
  Return the logging level
*/
extern int F77NAME(c_loglevel)() ;

#define logLevel F77NAME(c_loglevel)


/*
  Change whether the log file is being buffered.  The default is unbuffered.
     ibuff = nonzero to buffer the log, 0 to flush and unbuffer the log file.
*/
extern void F77NAME(c_setlogbuffered)(int* ibuff) ;

#define setLogBuffered F77NAME(c_setlogbuffered)


/*
  Change whether the log messages are time stamped.  The default is time stamped.
     itime = nonzero to time stamp the log messages, 0 for no time stamp
*/
extern void F77NAME(c_setlogtime)(int* itime) ;

#define setLogTime F77NAME(c_setlogtime)


/*
  Write out an info level message to the log file.
    msg = message
*/
extern void F77NAME(c_infolog)(char* msg) ;

#define infoLog F77NAME(c_infolog)


/*
  Write out an info level message and an integer to the log file.
    msg    = message
    ivalue = integer value
*/
extern void F77NAME(c_infologi)(char* msg, int* ivalue) ;

#define infoLogi F77NAME(c_infologi)


/*
  Write out an info level message and a double to the log file.
    msg    = message
    dvalue = double value
*/
extern void F77NAME(c_infologd)(char* msg, double* dvalue) ;

#define infoLogd F77NAME(c_infologd)


/*
  Write out a warn level message to the log file.
    msg = message
*/
extern void F77NAME(c_warnlog)(char* msg) ;

#define warnLog F77NAME(c_warnlog)


/*
  Write out a warn level message and an integer to the log file.
    msg    = message
    ivalue = integer value
*/
extern void F77NAME(c_warnlogi)(char* msg, int* ivalue) ;

#define warnLogi F77NAME(c_warnlogi)


/*
  Write out a warn level message and a double to the log file.
    msg    = message
    dvalue = double value
*/
extern void F77NAME(c_warnlogd)(char* msg, double* dvalue) ;

#define warnLogd F77NAME(c_warnlogd)


/*
  Write out an error level message to the log file.
    msg = message
*/
extern void F77NAME(c_errorlog)(char* msg) ;

#define errorLog F77NAME(c_errorlog)


/*
  Write out an error level message and an integer to the log file.
    msg    = message
    ivalue = integer value
*/
extern void F77NAME(c_errorlogi)(char* msg, int* ivalue) ;

#define errorLogi F77NAME(c_errorlogi)


/*
  Write out an error level message and a double to the log file.
    msg    = message
    dvalue = double value
*/
extern void F77NAME(c_errorlogd)(char* msg, double* dvalue) ;

#define errorLogd F77NAME(c_errorlogd)


#endif
