      subroutine fdelete(filename,ier)
      use iso_c_binding, only: c_char
C
C  dmc 18 June 1996 -- delete file
C
      character*(*) filename  ! file to delete...
C
C  ** delete a file **
C
C  return error code
C    ier=0 -- success
C    ier=1 -- error
C
C
      character(kind=c_char) fnbuf(500)
      integer cdelete
C----------------------------------
C
      ier=0
C
      ilf=len_trim(filename)
      if(ilf.eq.0) then
         write(6,*) ' ?fdelete:  filename argument is blank.'
         ier=1
         return
      endif
C
      call cstring(filename(1:ilf),fnbuf,'2C')
      ier=cdelete(fnbuf)
C
      return
      end
