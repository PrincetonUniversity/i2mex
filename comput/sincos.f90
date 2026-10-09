!******************** START FILE SINCOS.FOR ; GROUP TKODE2 ******************
!.............................................................

      subroutine SINCOS(ZTHETA,JANG,SNTHTK,CSTHTK)
use iso_c_binding, only: fp => c_double

!
!	THIS subroutine CALCULATES:
!		SIN(ZTHETA)
!		COS(ZTHETA)
!		SIN(2*ZTHETA)
!		COS(2*ZTHETA)
!	ETC., UP TO JANG*ZTHETA
!
!
!		ZTHETA IS INPUT ANGLE
!		JANG IS HIGHEST N*ZTHETA TO CALCULATE
!		SNTHTK IS THE ARRAY CONTAINING SIN(N*ZTHETA)
!		CSTHTK IS THE ARRAY CONTAINING COS(N*ZTHETA)
!
      DIMENSION SNTHTK(JANG),CSTHTK(JANG)
!
!  LOCAL MEMORY (DMC 6 JUL 1994)
!
      DATA ZTHETAP/0.0/
      DATA ZSINP/0.0/
      DATA ZCOSP/1.0/
!
      SAVE ZTHETAP,ZSINP,ZCOSP
!
!--------------------------------------------------------------------
!
!  DMC -- USE LOCAL MEMORY FOR SPEED
!
      if(ZTHETA.NE.ZTHETAP) THEN
!
!  EVALUATE SIN,COS
!
        SNTHTK(1)=SIN(ZTHETA)
        CSTHTK(1)=COS(ZTHETA)
        ZTHETAP=ZTHETA
        ZSINP=SNTHTK(1)
        ZCOSP=CSTHTK(1)
      else
!
!  REUSE PREVIOUS RESULTS
!
        SNTHTK(1)=ZSINP
        CSTHTK(1)=ZCOSP
      end if
!
      do 100 I=2,JANG
        SNTHTK(I)=SNTHTK(I-1)*CSTHTK(1)+CSTHTK(I-1)*SNTHTK(1)
        CSTHTK(I)=CSTHTK(I-1)*CSTHTK(1)-SNTHTK(I-1)*SNTHTK(1)
100   continue
!
      return
      end
!
!******************** end FILE SINCOS.FOR ; GROUP TKODE2 ******************
