      subroutine fclean_dir(path,ier)
C
C  dmc 16 June 2011 -- delete all files in a directory (CAREFUL!)
C     "hidden" files and subdirectories are not touched.
C
      use iso_c_binding, only: c_char
      implicit NONE
C
      character*(*), intent(in) :: path  ! directory to be cleaned...
      integer, intent(out) :: ier        ! completion status, 0=normal
C
C  if (path) is blank, clean out all files in $cwd !!
C
C  return error code
C    ier=0 -- success
C    ier=1 -- error
C
      integer ilp
C
      character(kind=c_char) cpath(500)
      integer cclean_dir
C
C----------------------------------
C
      ier=0
C
      ilp=len_trim(path)
      if(ilp.eq.0) then
         write(6,*) ' ?fclean_dir:  directory path argument is blank.'
         ier=1
         return
      endif
C
      if(path.eq.' ') then
         call cstring('./',cpath,'2C')
      else if(path(ilp:ilp).ne.'/') then
         call cstring(path(1:ilp)//'/',cpath,'2C')
      else
         call cstring(path(1:ilp),cpath,'2C')
      endif

      ier=cclean_dir(cpath)
C
      return
      end
