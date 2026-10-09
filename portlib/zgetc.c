#include "f77name.h"
 
#include <stdio.h>

#ifdef __LINUX
#undef SYS5
#endif
 
#ifdef SYS5
#include <sys/termio.h>
#else		/* POSIX */
#include <termios.h>
#endif
 
#ifdef SYS5
struct termio init_tty, new_tty;
#else
struct termios init_tty, new_tty;
#endif
 
char F77NAME(zgetc)()
{
  int istat;
  char c;
 
#ifdef SYS5
  ioctl(0, TCGETA, &init_tty);
  ioctl(0, TCGETA, &new_tty);
#else
  tcgetattr(0, &init_tty);
  tcgetattr(0, &new_tty);
/*  ioctl(0, TCGETP, &init_tty); */
/*  ioctl(0, TCGETP, &new_tty); */
#endif
  new_tty.c_lflag = new_tty.c_lflag & ~(ECHO); /* turn off echo */
  new_tty.c_lflag = new_tty.c_lflag & ~(ICANON); /* turn off canonical */
  new_tty.c_cc[VMIN] = 1;		/* one character input */
  new_tty.c_cc[VTIME] = 0;	/* timer off */
#ifdef SYS5
  ioctl(0, TCSETA, &new_tty);
#else
  tcsetattr(0, TCSANOW, &new_tty);
/*  ioctl(0, TCSANOW, &new_tty); */
#endif

  istat = fflush(NULL);
  c = getchar();

/* canon(); reset to canonical mode */
#ifdef SYS5
  ioctl(0, TCSETA, &init_tty);
#else
  tcsetattr(0, TCSANOW, &init_tty);
/*  ioctl(0, TCSANOW, &init_tty); */
#endif

  return c;

}
