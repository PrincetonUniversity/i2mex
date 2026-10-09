C******************** START FILE DMGFXT.FOR ; GROUP PLDMGR ******************
C-------------------------------------------------------------
C  DMGFXT
C
C  READ A FCN OF TIME+ADDL COORDINATE INTO DATA AREA
C  OR UPDATE ACCESS CODE IF DATA IS ALREADY IN DATA AREA
C
C History:
C 11/19/97 DMC : support multiple runs
C 09/24/97 CAL : Handle netCDF files
C-------------------------------------------------------------
      SUBROUTINE DMGFXT(JFCN,IND)
C
      use datmgr_mod
      use cplotr_mod
      implicit none
C
      character*50 zdlbl
C
      logical inlcdf, inlmds, itransp
C
      real,allocatable :: zxbuf(:)
      
      integer JFCN, iidcdf, ISIZE, JTYP, INX, IREC, IT, IND, IPT, ier

      allocate(zxbuf(max(16384,NR0)))
C------------------------------------
C
C  JFCN-- FCN NUMBER (input)
C
C  FIRST CHECK IF DATA IS ALREADY AVAILABLE
C
      if(lrun_x.eq.0) then
         zdlbl=abr(jfcn)
         inlmds=nlmds
         inlcdf=nlcdf
         iidcdf=idcdf
      else
         zdlbl=rlbl(lrun_x)(1:lrlbl(lrun_x))//'!'//abr_x(jfcn,lrun_x)
         inlmds=nlmds_x(lrun_x)
         inlcdf=nlcdf_x(lrun_x)
         iidcdf=idcdf_x(lrun_x)
      endif
C
      itransp = transp_imbed .and. (lrun_x.eq.0)
C
      CALL DMDLOC(zdlbl,IND,ISIZE,IPT)
      IF(IND.GT.0) goto 999
C
C  DATA NEEDS TO BE READ IN-- ALLOCATE SPACE
C
      if(lrun_x.eq.0) then
         JTYP=ITYPR(JFCN)
         INX=NZONEX(JTYP)
         ISIZE=NTR*INX
      else
         JTYP=ITYPR_X(JFCN,lrun_x)
         INX=NZONEX_X(JTYP,lrun_x)
         ISIZE=NTR_X(lrun_x)*INX
      endif
C
      CALL RP_DMGALO(ISIZE,IND,5)
C
C  LABEL SPACE
      DMGLBL(IND)=zdlbl
C
C  READ IN FCN
C
      ier=0
      if (inlcdf) then
         call cdfmfrd(iidcdf, JFCN, ind, ier)
      else if (inlmds) then
         call mdsmfrd(JFCN, ind, ier)
         if (ier .ne. 0) then
            ntr=-1                      ! cludge to return error
         endif
      else if (itransp) then
         ipt=locd(ind)
         datbuf(ipt:ipt+isize-1)=0      ! TRANSP will provide the data later.
         nwds(ind)=isize
      else
         if(lrun_x.eq.0) then
            IREC=NTR+NTCORR
         else
            IREC=NTR_X(lrun_x)
         endif
         DO 100 IT=1,IREC
            CALL READMF(IT,JFCN,ZXBUF,INX)
            IF(LTWRIT(IT).or.(lrun_x.gt.0))
     1	 CALL DATADD(IND,ZXBUF,INX)
 100     CONTINUE
      end if
 999  continue
      deallocate(zxbuf)
      RETURN
      END
C******************** END FILE DMGFXT.FOR ; GROUP PLDMGR ******************
