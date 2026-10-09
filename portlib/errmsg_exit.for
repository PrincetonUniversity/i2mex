      subroutine errmsg_exit(msg)
      character*(*) msg
C
C  write message and exit
C
      write(0,*) ' '
      write(0,*) msg
      write(0,*) ' '
C
      call bad_exit
      return
      end
