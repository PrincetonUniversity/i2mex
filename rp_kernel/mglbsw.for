      subroutine mglbsw(iint,zunits,isw)
C
      character*(*) zunits
C
C  RPLOT units transformation routine
C
C input:
C  iint -- operator code (cf trread/rpfixuns for translation)
C in/out:
C  zunits -- units to be transformed
C out:
C  isw -- 0 if units transform could not be done
C         1 if units transform was done
C
C-----------------------------------------
C  local...
C
      CHARACTER*10 LBLSW(2,10,5)
C
      INTEGER ILSWA(5)
C
C  STANDARD LABEL FOLLOWED BY LABEL TO BE USED WHEN INTEGRATED
C  FUNCTION IS GRAPHED
C
      DATA LBLSW/'N/CM3     ','N         ',
     >             'N/CM**3   ','N         ',
     >             'WATTS/CM3 ','WATTS     ',
     >             'N/CM3/SEC ','N/SEC     ',
     >             'JLES/CM3  ','JOULES    ',
     >             'GRAMS/CM3 ','GRAMS     ',
     >             'Nt-M/CM3  ','Nt-M      ',
     >             'NtM-S/CM3 ','Nt-M-SEC  ',
     >             'NtMS2/CM3 ','Nt-M-SEC2 ',
     1             'EV        ','EV*CM**3  ',       ! Tbt 9/22/92
     >             'N/CM3/SEC ','N/CM2/SEC ',
     >             'WATTS/CM3 ','WATTS/CM2 ',
     >             'Nt-M/CM3  ','Nt-M/CM2  ',14*'          ',
     >             'AMPS/CM2  ','AMPS      ',
     >             'AMPS/CM2  ','AMPS      ',16*'          ',
     >             'N/CM3     ','N/CM**4   ',
     >             'N/CM**3   ','N/CM**4   ',
     >             'TESLA     ','TESLA/CM  ',
     >             '          ','CM**-1    ',
     >             'AMPS/CM2  ','AMPS/CM3  ',
     >             'EV        ','EV/CM     ',
     >             'PASCALS   ','PASCLS/CM ',
     >             'JLES/CM3  ','JLES/CM4  ',
     >             'VOLTS     ','VOLTS/CM  ',
     >             'SEC**-1   ','1/SEC/CM  ',
     >             'N/CM3     ','N/CM3/SEC ',   ! dmc Feb 1996
     >             'N/CM**3   ','N/CM3/SEC ',
     >             'TESLA     ','TESLA/SEC ',
     >             '          ','SEC**-1   ',
     >             'JLES/CM3  ','WATTS/CM3 ',
     >             'GRAMS/CM3 ','G/CM3/SEC ',
     >             'PASCALS   ','PAS/SEC   ',
     >             'EV        ','EV/SEC    ',
     >             'SEC**-1   ','SEC**-2   ',
     >             '          ','          '/
 
C
C  NUMBER OF ACTIVE LABEL SWITCHES BY OPERATION: ***
      DATA ILSWA/10,3,2,10,9/
C
C
C-----------------------------------------
C
      IF(IINT.EQ.5) THEN
        ZUNITS='CM**-1 '
        RETURN
      ELSE IF(IINT.EQ.6) THEN
        ZUNITS='CM '
        RETURN
      ENDIF
C
      IJINT=IINT
      if(IJINT.eq.9) IJINT=5
C
      CALL GLBSW(LBLSW(1,1,IJINT),ILSWA(IJINT),ZUNITS,ISW)
C
      return
      end
