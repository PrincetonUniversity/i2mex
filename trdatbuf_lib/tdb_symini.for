      subroutine tdb_symini(d,nzones)
      !
      !  set up workspace for presymmetrization & mapping of profiles
      !  size of workspace depends on target application's number of zones
      !  and the number of zones in the input data
      !

      use trdatbuf_obj
      implicit NONE
      type (trdatbuf) :: d
      integer, intent(in) :: nzones   ! no. of radial zones in caller's grid

      !----------------------------
      integer :: iwrk,i,j
      !----------------------------

      IWRK=2*(NZONES+2)
C
C=>TRDATGEN+
!
!    ******************************************
!    * TRDATGEN GENERATED CODE -- DO NOT EDIT *
!    ******************************************
!
!    code generation ends at C=>TRDATGEN- line
!
!
C
      IF((D%LFBOL.GT.0).AND.(D%NSYBOL.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXBOL + 1))
C
      IF((D%LFBPA.GT.0).AND.(D%NSYBPA.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXBPA + 1))
C
      IF((D%LFBPB.GT.0).AND.(D%NSYBPB.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXBPB + 1))
C
      IF((D%LFD2F.GT.0).AND.(D%NSYD2F.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXD2F + 1))
C
      IF((D%LFDE2.GT.0).AND.(D%NSYDE2.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDE2 + 1))
C
      IF((D%LFDF3.GT.0).AND.(D%NSYDF3.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDF3 + 1))
C
      IF((D%LFDF4.GT.0).AND.(D%NSYDF4.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDF4 + 1))
C
      IF((D%LFDF6.GT.0).AND.(D%NSYDF6.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDF6 + 1))
C
      IF((D%LFDFD.GT.0).AND.(D%NSYDFD.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDFD + 1))
C
      IF((D%LFDFH.GT.0).AND.(D%NSYDFH.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDFH + 1))
C
      IF((D%LFDFT.GT.0).AND.(D%NSYDFT.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXDFT + 1))
C
      IF((D%LFECF.GT.0).AND.(D%NSYECF.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXECF + 1))
C
      IF((D%LFGRB.GT.0).AND.(D%NSYGRB.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXGRB + 1))
C
      IF((D%LFLF3.GT.0).AND.(D%NSYLF3.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLF3 + 1))
C
      IF((D%LFLF4.GT.0).AND.(D%NSYLF4.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLF4 + 1))
C
      IF((D%LFLF6.GT.0).AND.(D%NSYLF6.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLF6 + 1))
C
      IF((D%LFLFD.GT.0).AND.(D%NSYLFD.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLFD + 1))
C
      IF((D%LFLFH.GT.0).AND.(D%NSYLFH.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLFH + 1))
C
      IF((D%LFLFT.GT.0).AND.(D%NSYLFT.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXLFT + 1))
C
      IF((D%LFNER.GT.0).AND.(D%NSYNER.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNER + 1))
C
      IF((D%LFNI3.GT.0).AND.(D%NSYNI3.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNI3 + 1))
C
      IF((D%LFNI4.GT.0).AND.(D%NSYNI4.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNI4 + 1))
C
      IF((D%LFNI6.GT.0).AND.(D%NSYNI6.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNI6 + 1))
C
      IF((D%LFNID.GT.0).AND.(D%NSYNID.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNID + 1))
C
      IF((D%LFNIH.GT.0).AND.(D%NSYNIH.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNIH + 1))
C
      IF((D%LFNIM.GT.0).AND.(D%NSYNIM.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNIM + 1))
C
      IF((D%LFNIT.GT.0).AND.(D%NSYNIT.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNIT + 1))
C
      IF((D%LFNMR.GT.0).AND.(D%NSYNMR.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXNMR + 1))
C
      IF((D%LFOMG.GT.0).AND.(D%NSYOMG.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXOMG + 1))
C
      IF((D%LFPRS.GT.0).AND.(D%NSYPRS.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXPRS + 1))
C
      IF((D%LFQPR.GT.0).AND.(D%NSYQPR.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXQPR + 1))
C
      IF((D%LFSBI.GT.0).AND.(D%NSYSBI.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXSBI + 1))
C
      IF((D%LFTER.GT.0).AND.(D%NSYTER.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXTER + 1))
C
      IF((D%LFTI2.GT.0).AND.(D%NSYTI2.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXTI2 + 1))
C
      IF((D%LFTI3.GT.0).AND.(D%NSYTI3.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXTI3 + 1))
C
      IF((D%LFTQI.GT.0).AND.(D%NSYTQI.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXTQI + 1))
C
      IF((D%LFV2F.GT.0).AND.(D%NSYV2F.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXV2F + 1))
C
      IF((D%LFVB2.GT.0).AND.(D%NSYVB2.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVB2 + 1))
C
      IF((D%LFVC3.GT.0).AND.(D%NSYVC3.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVC3 + 1))
C
      IF((D%LFVC4.GT.0).AND.(D%NSYVC4.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVC4 + 1))
C
      IF((D%LFVC6.GT.0).AND.(D%NSYVC6.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVC6 + 1))
C
      IF((D%LFVCD.GT.0).AND.(D%NSYVCD.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVCD + 1))
C
      IF((D%LFVCH.GT.0).AND.(D%NSYVCH.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVCH + 1))
C
      IF((D%LFVCT.GT.0).AND.(D%NSYVCT.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVCT + 1))
C
      IF((D%LFVEE.GT.0).AND.(D%NSYVEE.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVEE + 1))
C
      IF((D%LFVIE.GT.0).AND.(D%NSYVIE.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVIE + 1))
C
      IF((D%LFVMO.GT.0).AND.(D%NSYVMO.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVMO + 1))
C
      IF((D%LFVP2.GT.0).AND.(D%NSYVP2.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVP2 + 1))
C
      IF((D%LFVPO.GT.0).AND.(D%NSYVPO.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVPO + 1))
C
      IF((D%LFVPR.GT.0).AND.(D%NSYVPR.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVPR + 1))
C
      IF((D%LFVTR.GT.0).AND.(D%NSYVTR.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXVTR + 1))
C
      IF((D%LFZF2.GT.0).AND.(D%NSYZF2.EQ.1))
     >    IWRK = MAX(IWRK,(D%NXZF2 + 1))
C
      DO I = 1, D%NNSIM
        J = D%NISSIM(I)
        IF((D%LFSIM(I).GT.0).AND.(D%NSYSIM(J).EQ.1))
     >      IWRK = MAX(IWRK,(D%NXSIM(J) + 1))
      END DO
C
C
C=>TRDATGEN-
C
      DO 5 I=1,8
        d%NBX(I)=IWRK
        CALL TDB_WORKALO(d,IWRK,d%LBX(I))
 5    CONTINUE
C
      return
      end
