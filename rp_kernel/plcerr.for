      subroutine plcerr(ier,zinput,ilast)
C
C  after error message
C  report position of a calculator error -- from parser or evaluator
C
      character*(*) zinput
C
      CHARACTER*1  HYPHEN
      CHARACTER*1  CARROT
C
      data hyphen/'-'/
      data carrot/'^'/
C
C ====================
C
      IF ( ILAST .LE. 70)  THEN
         WRITE(lunzer(0),8126) ZINPUT(1:ILAST)
 8126    FORMAT(/' ?PLCFXT: POSITION OF',
     1      ' FIRST ERROR IN INPUT LINE:'
     1      /1X, A)
         WRITE(lunzer(0),8127) (HYPHEN, I=3,IER), CARROT
 8127    FORMAT(1X, 79A1)
      ELSE                              ! INPUT IS MORE THAN ONE LINE
         WRITE(lunzer(0),8126) ZINPUT(1:70)
         IF (IER .LE. 72)  THEN
C		    .ERROR IS IN FIRST PART OF INPUT.
      	    WRITE(lunzer(0),8127) (HYPHEN, I=3,IER), CARROT
            WRITE(lunzer(0),8128) ZINPUT(71:ILAST)
         ELSE
C		    .ERROR IS IN SECOND PART OF INPUT.
      	    WRITE(lunzer(0),8128) ZINPUT(71:ILAST)
      	    WRITE(lunzer(0),8127) (HYPHEN, I=73,IER), CARROT
 8128       FORMAT(1X, A)
         END IF                         ! IER
 
      END IF                            ! ILAST
      WRITE(lunzer(0),8290)
 8290 FORMAT(/'?PLCFXT: Input line ignored. $ is unchanged'/)
 
      return
      end
