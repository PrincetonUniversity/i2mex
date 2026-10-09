      subroutine rpbufsiz(zroutine,ihave,ineed)
C
C  subroutine to write an error message -- passed buffer too small
C
C  all input:
      character*(*) zroutine            ! name of routine with error
      integer ihave                     ! buffer size passed
      integer ineed                     ! buffer size needed
C
C----------------------
C
      write(lunzer(0),9001) zroutine,ineed,ihave
 9001 format(' ?',a,' -- buffer size too small; need = ',i6,
     >   ' passed = ',i6)
C
      return
      end
