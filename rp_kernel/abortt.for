C----------------------------------------------------------------------
C  ABORT An RPLOT RUN AFTER FATAL ERROR CONDITION IS DETECTED
C    FORCE AN ARITHMETIC ERROR WHICH THE JOB CONTROL CODE SHOULD BE
C    ABLE TO DETECT.  Created 1/21/94 from ABORTR by T.B. Terpstra
C
      subroutine abortt
 
      call bad_exit
 
      END
