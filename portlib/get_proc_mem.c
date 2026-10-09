/**************************************************************************
get_proc_mem(float *v_mem, float *r_mem)

05/15/2009 CLF: modified get_vmem_()
                to return vmem & rss in KB as "floating point"


S. Ethier (ethier@pppl.gov)

06/18/2003     FUNCTION get_vmem_()
              ======================
  This function returns the size of the virtual memory used by the current
  process. It can be called from a Fortran or C code and it basically
  reads the /proc/$PID/stat file to get the information, just like the
  "ps" and "top" commands.

  This file can be compiled as a stand-alone test program with:
         cc -DSTANDALONE get_mem.c -o read_stat

07/19/2004  S. Ethier  Update:

  The function reads many other quantities from /proc/$PID/stat. You just
  need to change the return value (ex. "return rss;" to return the
  resident memory size). See the STANDALONE section for examples...

 
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
/* #include <asm/page.h>  */    /* for PAGE_SHIFT and PAGE_SIZE (in bytes) */
#include <sys/param.h>          /* for HZ (number of clock cycles per sec)*/

static pid_t pid;
static char  stat_file[128];
static int   pidnotknown=1;  /* pidnotknown=1 is true, used
                                to get the PID only once */

void get_proc_mem_(float *v_mem, float *r_mem)
{
   char state[8], comm[80];
   unsigned long i;
   FILE *fd;
   int            pid2, ppid, pgrp, session, tty, tpgid;
   unsigned long  flags, minflt, cminflt, majflt, cmajflt;
   int            utime, stime, cutime, cstime, counter, priority;
   unsigned long  timeout, itrealvalue;
   int            starttime;
   unsigned long  rss, rlim;
   float          vsize;
   int            res_mem;

   long sz = sysconf(_SC_PAGESIZE);

   if (pidnotknown) {
      pid = getpid();
      sprintf(stat_file, "/proc/%d/stat", pid);
      pidnotknown = 0;
   }

   fd = fopen(stat_file, "r");
   fscanf(fd, "%d%s%s%d%d%d%d%d%u%u%u%u%u%d%d%d%d%d%d%u%u%u%g%u%u",
         &pid2, comm, state, &ppid, &pgrp, &session, &tty, &tpgid,
         &flags, &minflt, &cminflt, &majflt, &cmajflt,
         &utime, &stime, &cutime, &cstime, &counter, &priority,
         &timeout, &itrealvalue,
         &starttime,
         &vsize, &rss, &rlim);
   fclose(fd);
 
#ifdef STANDALONE
   printf("pid = %d     stat_file = %s\n", pid, stat_file);
   printf("pid2 = %d  comm = %s  state = %s\n", pid2, comm, state);
   printf("ppid = %d  pgrp = %d  session = %d  tty = %d  tpgid = %d\n",
          ppid, pgrp, session, tty, tpgid);
   printf("flags=%u  minflt=%u  cminflt=%u  majflt=%u  cmajflt=%u\n",
          flags, minflt, cminflt, majflt, cmajflt);
   printf("utime=%d  stime=%d  cutime=%d  cstime=%d  counter=%d  priority=%d\n",
          utime, stime, cutime, cstime, counter, priority);
   printf("timeout=%u  itrealvalue=%u  starttime=%d\n",
          timeout, itrealvalue, starttime);
   printf("\nvsize=%g  rss=%u  rlim=%u\n",
          vsize, rss, rlim);
   /* res_mem = rss << (PAGE_SHIFT -10); */
   res_mem = (int)(rss*sz/1024);
   *v_mem = (vsize/1024.); /* kB */ 
   printf("  resident memory = %d kB  virtual memory size = %g kB \n\n\n",
          res_mem, *v_mem);
   /* sleep(10); */

#endif
   /*   res_mem = rss << (PAGE_SHIFT -10); *//* kB */ 
   res_mem=(int)rss;
   *r_mem= (float)(res_mem*sz)/1024.;
   *v_mem = (vsize/1024.); /* kB */ 

   /*   printf("  resident memory = %g kB  %d %d\n\n\n",
          *r_mem, res_mem, sz);
  */
   return;
}


#ifdef STANDALONE
int main()
{
   float vm;
   float rm;
   long i;
   float  *yy;

   get_proc_mem_(&vm,&rm);
   printf(" >> main v = %g r = %g\n",vm,rm);
   yy = calloc(262144, sizeof(float));
   get_proc_mem_(&vm,&rm);
   printf(" >> main v = %g r = %g\n",vm,rm);
   for (i = 0; i < 262144; i++) {
      yy[i] = 1.;
   }
   get_proc_mem_(&vm,&rm);
   printf(" >> main v = %g r = %g\n",vm,rm);
   free(yy);
   get_proc_mem_(&vm,&rm);
   printf(" >> main v = %g r = %g\n",vm,rm);

   return;
}
#endif


void get_proc_mem(float *v_mem, float *res_mem)
{
   get_proc_mem_(v_mem, res_mem);
}

 
/***********************************************************************
  stat   Status information about the process.  This is used by ps(1).

  The fields, in order, with their proper scanf(3) format specifiers, are:

     pid %d The process id.

     comm %s
             The filename of the executable, in parentheses.  This is visible
             whether or not the executable is swapped out.

     state %c
             One character from the string "RSDZT" where R is running, S is
             sleeping in an interruptible wait, D is sleeping in an
             uninterruptible wait or swapping, Z is zombie, and T is traced
             or stopped (on a signal).

     ppid %d
             The PID of the parent.

     pgrp %d
             The process group ID of the process.

     session %d
             The session ID of the process.

      tty %d The tty the process uses.

     tpgid %d
             The process group ID of the process which currently owns the tty
             that the process is connected to.

     flags %u
             The flags of the process.  Currently, every flag has the
             math bit set, because crt0.s checks for math emulation,
             so this is not included in the output.  This is probably
             a bug, as not every process is a compiled C program.  The
             math bit should be a decimal 4, and the traced bit is
             decimal 10.

     minflt %u
             The number of minor faults the process has made, those
             which have not required loading a memory page from disk.

     cminflt %u
             The number of minor faults that the process and its
             children have made.

     majflt %u
             The number of major faults the process has made, those
             which have required loading a memory page from disk.

     cmajflt %u
             The number of major faults that the process and its
             children have made.

     utime %d
             The number of jiffies that this process has been
             scheduled in user mode.

     stime %d
             The number of jiffies that this process has been
             scheduled in kernel mode.

     cutime %d
             The number of jiffies that this process and its children
             have been scheduled in user mode.

     cstime %d
             The number of jiffies that this process and its children
             have been scheduled in kernel mode.

     counter %d
             The current maximum size in jiffies of the process's next
             timeslice, or what is currently left of its current
             timeslice, if it is the currently running process.

     priority %d
             The standard nice value, plus fifteen.  The value is
             never negative in the kernel.

     timeout %u
             The time in jiffies of the process's next timeout.

     itrealvalue %u
             The time (in jiffies) before the next SIGALRM is sent to
             the process due to an interval timer.

     starttime %d
             Time the process started in jiffies after system boot.

     vsize %u
             Virtual memory size

      rss %u Resident Set Size: number of pages the process has in
             real memory, minus 3 for administrative purposes. This is
             just the pages which count towards text, data, or stack
             space.  This does not include pages which have not been
             demand-loaded in, or which are swapped out.

     rlim %u
             Current limit in bytes on the rss of the process (usually
             2,147,483,647).

     startcode %u
             The address above which program text can run.

     endcode %u
             The address below which program text can run.

     startstack %u
             The address of the start of the stack.

     kstkesp %u
             The current value of esp (32-bit stack pointer), as found
             in the kernel stack page for the process.

     kstkeip %u
             The current EIP (32-bit instruction pointer).

     signal %d
             The bitmap of pending signals (usually 0).

     blocked %d
             The bitmap of blocked signals (usually 0, 2 for shells).

     sigignore %d
             The bitmap of ignored signals.

     sigcatch %d
             The bitmap of catched signals.

     wchan %u
             This is the "channel" in which the process is waiting.
             This is the address of a system call, and can be looked
             up in a namelist if you need a textual name.  (If you
             have an up-to-date /etc/psdatabase, then try ps -l to see
             the WCHAN field in action)
***************************************************************************/

