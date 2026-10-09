      subroutine plftmck(zsubr,lunt,ztime,zdata,inpts,iordr,zamax)
C
C  dmc -- break shared code out of various plftmk*.for routines in
C    rplot_sub
C
C    find max abs. value of data.
C    check time axis ordering of data
C
C  input:
      character*(*) zsubr               ! name of caller (for message)
      integer lunt                      ! i/o unit number for message
      integer inpts                     ! no. of data pts
      real ztime(inpts)                 ! timebase
      real zdata(inpts)                 ! the f(t) data
C
C  output:
C
      real zamax                        ! max(abs(zdata(i)),i=1 to INPTS)
      integer iordr                     ! no. of timebase pts out of order
C
C-------------------
C
C
C------------------------------------
C
      ZAMAX=ABS(ZDATA(1))
      IORDR=0
C
      IF(INPTS.GT.1) THEN
        DO IT=2,INPTS
          ZAMAX=AMAX1(ZAMAX,ABS(ZDATA(IT)))
          ITM1=IT-1
          IF(ZTIME(ITM1).GE.ZTIME(IT)) THEN
            IORDR=IORDR+1
            ils=max(1,len_trim(zsubr))
            WRITE(LUNT,9001) ITM1,ZTIME(ITM1),IT,ZTIME(IT),zsubr(1:ils)
 9001   FORMAT('  PT. NO. ',I6,' TIME VALUE = ',1PE11.4/
     >           '  PT. NO. ',I6,' TIME VALUE = ',1PE11.4/
     >		 ' %',A,' - TIME SERIES OUT OF ORDER')
          ENDIF
        ENDDO
      ENDIF
C
      return
      end
 
