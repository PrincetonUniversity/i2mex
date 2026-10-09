      SUBROUTINE SHOWDEFL(DEVNAM,LDEV,DIRNAM,LDIR)
C
      implicit none
C
C  GET THE CURRENT DEFAULT DISK AND DIRECTORY BY PARSING A DUMMY
C  FILENAME WITH THE RMS SERVICE CALL
C  DMC:  FOLLOWING HINTS BY M. THOMPSON
C
 
      CHARACTER*(*) DEVNAM  ! OUTPUT DEVICE NAME
      CHARACTER*(*) DIRNAM  ! OUTPUT DIRECTORY NAME
 
      INTEGER LDEV	      ! OUTPUT LENGTH OF DEVICE NAME
      INTEGER LDIR	      ! OUTPUT LENGTH OF DIRECTORY NAME
 
      integer istat
      integer getcwd
 
      ldev=0
      devnam=' '
C
      ldir=0
      dirnam=' '
      istat=getcwd(dirnam)
C
      call str_pad(dirnam)
C
      ldir=len_trim(dirnam)
C
      RETURN
      END
