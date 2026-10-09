      subroutine cstring(fstr,cstr,cmd)
C
C  Fortran-77 <-> C string converter
C
C  usage:  call cstring('this is a string',cstr,'2C')
C          ...fills the byte array cstr with a null terminated string
C             suitable for being passed to a C routine
C
C               **caution**
C               in  "call cstring(fstr,cstr,'2C')" trailing blanks in
C               fstr are also copied; if this is not desired use something
C               like
C                     ilstr=max(1,len_trim(fstr))
C                     call cstring(fstr(1:ilstr),cstr,'2C')
C
C
C          call cstring(fstr,cstr,'2F') copies the contents of the
C                byte array cstr into fstr, stopping when a terminating
C                null is encountered or the end of fstr is reached,
C                whichever happens first.  fstr is padded with blanks.
C
C==  NB 
C     call cstring( trim(fs), cs, '2C' )
C                                       IS INVALID & WILL FAIL if compiled in debug mode-
C     "There is an assignment to dummy argument 1, which is associated with a constant, 
C     "an expression that is not a variable..." 
C
C     Use cstr_2C( trim(fs), cs, '2C' )  instead..
C
      use iso_c_binding, only: c_char, c_null_char
      implicit none
C
 
C
C  fortran string:
C
      character*(*) fstr
C
C  C string (null terminated byte array)
C
      character(kind=c_char) cstr(*)
C
C  cmd string:  "2F" for cstr->fstr copy, "2C" for fstr->cstr copy.
C
      character*2 cmd
C
C--------------------------
C  local:
C
      integer i,ilen
C------------------------------------------------
C
C
      if(cmd.eq.'2F')                             then
         fstr=' '
         do i=1,len(fstr)
            if(cstr(i).eq.c_null_char) go to 10             ! stop at terminating null
            fstr(i:i)=cstr(i)
         enddo
 10      continue
C
      elseif(cmd.eq.'2C')                         then
         call cstr_2C( fstr, cstr, cmd )
                                                  else
         call errmsg_exit(' ??cstring:  invalid cmd value:  '//cmd)
         call bad_exit
      endif
C
      return
      end
      subroutine cstr_2C( fstr, cstr, cmd )
C.. cstr_2C                             Legally convert trim(fstr) to cstr
C                                       codesys/source/portlib/cstring
C
C
C     Author                            Jim.Conboy@ccfe.ac.uk
C     Version                           1.00,  22Dec2008
C     Modifications
C
C     1.00  22Dec2008                   Jim.Conboy@ccfe.ac.uk
C                                       Extracted from cstring
C----_^---------------------------------=========================================|
C
      use iso_c_binding, only: c_char, c_null_char
      implicit  none
C
C
C..Arguments
C
      character(len=*), intent(in)   :: fstr           ! the fortran string
      character(kind=c_char)         :: cstr(*)        ! the 0 terminated c string
      character(len=2), intent(in)   :: cmd            ! bkwd compatible..
C
C..local
      integer                        :: i              ! loop index 
     +                                 ,ilen           ! len string
C
C     character(len=*),parameter     :: $cr='cstr_2C'
C----_^------------------------------::================!.........................|
C
C     if( present(cmd))                then
C--      Check cmd if present, in case of editing error..
         if( cmd .ne. '2C')        then
         call errmsg_exit(' ??cstr_2C:  invalid cmd value:  '//cmd)
         call bad_exit
                                   endif
C                                      endif
C
         ilen=len(fstr)
         do i=1,ilen
            cstr(i)=fstr(i:i)
         enddo
         cstr(ilen+1)=c_null_char
         return
      end subroutine cstr_2C
