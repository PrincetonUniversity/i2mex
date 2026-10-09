      subroutine str_pad(str)
C
C  replace nulls in str with blanks
C
      character*(*) str
C
      ilen=len(str)
      do il=1,ilen
         if(ichar(str(il:il)).eq.0) str(il:il)=' '
      enddo
C
      return
      end
