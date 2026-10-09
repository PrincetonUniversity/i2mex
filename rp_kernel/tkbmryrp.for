C-----------------------------------------------------------------
C  DMC - UTILITY ROUTINE FOR MOMENT SURFACES EXTRAPOLATION.
C   RPLOT VERSION - ADAPTED FROM TRANSP - DMC DEC 1988
C
C  RETURN THE R,Y PT. FOR GIVEN SURFACE, PASSED THETA VALUE
C
      SUBROUTINE TKBMRYRP(ISURF,ZTHETA,ZROUT,ZYOUT)
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
C--------------------------------------------
C
      CALL SINCOS ( ZTHETA, IMOM, SNTHTK, CSTHTK )
C
      if(immgeo.eq.0) then
C
C  SYMMETRIC GEOMETRY
C
        ZROUT = ZRMM0A(ISURF)
        ZYOUT = 0.0
C
        DO 10560 IMN = 1,IMOM
      	ZRMB2X = ZRMOMA(IMN, ISURF)
      	ZYMB2X = ZYMOMA(IMN, ISURF)
              ZROUT = ZROUT + ZRMB2X * CSTHTK(IMN)
              ZYOUT = ZYOUT + ZYMB2X * SNTHTK(IMN)
10560     CONTINUE
C
      else
C
C  ASYMMETRIC GEOMETRY
C
        ZROUT=ZRMCA(0,1,ISURF)
        ZYOUT=ZYMCA(0,1,ISURF)
        DO IMN=1,IMOM
           ZRMCL=ZRMCA(IMN,1,ISURF)
           ZRMSL=ZRMCA(IMN,2,ISURF)
           ZYMCL=ZYMCA(IMN,1,ISURF)
           ZYMSL=ZYMCA(IMN,2,ISURF)
           ZROUT=ZROUT+ZRMCL*CSTHTK(IMN)+ZRMSL*SNTHTK(IMN)
           ZYOUT=ZYOUT+ZYMCL*CSTHTK(IMN)+ZYMSL*SNTHTK(IMN)
        ENDDO
C
      endif
C
      RETURN
      END
C******************** END FILE TKBMRY.FOR ; GROUP TKBLOAT **************
