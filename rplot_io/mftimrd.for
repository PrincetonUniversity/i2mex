      subroutine mftimrd(mfluni,mfhdr,iblksz,time3,ntime3,ntr)
C
C  read "new format" MF.PLN time vector
C
      implicit NONE
C
      integer mfluni                    ! LUN of open MF.PLN file
      integer :: iblksz
      integer mfhdr(iblksz)             ! header record (already read in)
      integer :: ntime3
      real time3(ntime3)                ! time vector
      integer ntr                       ! actual no. of time pts used
C
C-------------------------
      integer :: itimes,itb,iblokt,irec,iwd1,iwd2,ii
      integer :: lunzer
C-------------------------
C
      itimes=mfhdr(2)
      if(itimes.gt.ntime3) then
         write(lunzer(0),8801) ntime3
 8801    format(' %mftimrd:  array limit = ',i5,' exceeded'/
     >          '  MF.PLN file timebase is truncated.')
         itimes=ntime3
      endif
C
      ntr=itimes
C
      IBLOKT=1+(NTR-1)/IBLKSZ
      DO ITB=1,IBLOKT
         IREC=MFHDR(3)+ITB-1
         IWD1=1+(ITB-1)*IBLKSZ
         IWD2=IWD1+IBLKSZ-1
         READ(MFLUNI,REC=IREC) (TIME3(II),II=IWD1,IWD2)
      ENDDO
C
      return
      end
 
