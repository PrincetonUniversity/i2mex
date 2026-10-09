      subroutine fcopy(oldname,newname,ier)
C
C  dmc 15 June 2011 -- copy file w/o system (shell) call out
C
      use iso_c_binding, only: c_char
      implicit NONE
C
      character*(*), intent(in) :: oldname   ! old filename
      character*(*), intent(in) :: newname   ! new filename
      integer, intent(out) :: ier  ! completion status, 0=OK
C
C  ** copy a file **
C
C  return error code
C    ier=0 -- success
C    ier=1 -- error
C
C  cccopy.c writes a message on stdout if there is an error opening
C    oldname or newname files
C
      integer ilo,iln,cccopy
C
C
      character(kind=c_char) fnbuf1(500),fnbuf2(500)
C----------------------------------
C
      ier=0
C
      ilo=len_trim(oldname)
      iln=len_trim(newname)
      if(min(ilo,iln).le.0) then
         write(6,*) ' ?fcopy:  filename argument is blank.'
         ier=1
         return
      endif
C
      call cstring(oldname(1:ilo),fnbuf1,'2C')
      call cstring(newname(1:iln),fnbuf2,'2C')

      ier=cccopy(fnbuf1,fnbuf2)
      if(ier.ne.0) ier=1
C
      return
      end
