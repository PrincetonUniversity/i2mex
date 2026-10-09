C******************** START FILE TKBRYS.FOR ; GROUP TKBLOAT ************
C-- END FILE #TKBRY2# ************************************
CFORMFEEDC-- START FILE #TKBRYS# ************************************
C-----------------------------------------------------------------
C  DMC - UTILITY ROUTINE FOR MOMENT SURFACES EXTRAPOLATION.
C
C  GIVEN A SURFACE AT THE SPECIFIED MOMENTS INDEX LOCATION THIS
C  ROUTINE RETURNS THE (R,Y) AS A FCN OF THETA, N PTS. FROM 0 TO
C  N-1/N * TWOPI; AND GIVES A ROUGH ESTIMATE OF THE CENTROID
C  (R0,Y0) DEFINED FROM AN APPROXIMATE PT TO PT LINE INTEGRATION.
C
      SUBROUTINE TKBRYSRP(ISURF,ZR,ZY,ZR0,ZY0,N)
C
      REAL ZR(N),ZY(N)
C
      DATA TWOPI/6.283185/
C
C------------------------------------
C
      CALL TKBMRYRP(ISURF,0.0,ZR(1),ZY(1))
C
      ZR0=0.0
      ZY0=0.0
      ZL=0.0
C
      DO 100 ITH=2,N
        ITHM1=ITH-1
C
        ZTH=TWOPI*FLOAT(ITH-1)/FLOAT(N)
        CALL TKBMRYRP(ISURF,ZTH,ZR(ITH),ZY(ITH))
C
        ZDL=SQRT((ZR(ITH)-ZR(ITHM1))**2+(ZY(ITH)-ZY(ITHM1))**2)
C
        ZL=ZL+ZDL
        ZR0=ZR0+ZDL*0.5*(ZR(ITH)+ZR(ITHM1))
        ZY0=ZY0+ZDL*0.5*(ZY(ITH)+ZY(ITHM1))
C
 100  CONTINUE
C
C  COMPLETE THE LOOP
C
        ZDL=SQRT((ZR(1)-ZR(N))**2+(ZY(1)-ZY(N))**2)
C
        ZL=ZL+ZDL
        ZR0=ZR0+ZDL*0.5*(ZR(1)+ZR(N))
        ZY0=ZY0+ZDL*0.5*(ZY(1)+ZY(N))
C
        ZR0=ZR0/ZL
        ZY0=ZY0/ZL
C
      RETURN
      END
C******************** END FILE TKBRYS.FOR ; GROUP TKBLOAT **************
