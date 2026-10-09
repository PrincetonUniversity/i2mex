      subroutine sset_cwd(dir,ierr)
      use iso_c_binding, only: c_char
      character*(*) dir                 ! directory to cd to...
      integer ierr                      ! error status code, 0 = normal
 
      character(kind=c_char) dirbuf(256)
 
      integer chdir
      integer ild
C
      ild=max(1,len_trim(dir))

      if(ild.le.255) then
         call cstring(dir(1:ild),dirbuf,'2C')
         ierr = linux_chdir(dirbuf)
      else
         write(6,*) ' ?sset_ccwd:  directory path too long:'
         write(6,*) dir(1:ild)
         ierr=1
      endif

      return
      end
