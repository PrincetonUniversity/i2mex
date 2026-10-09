      subroutine mfhdrd(mfluni,nfxt,ier)
      use mfblok_mod,only : iblksz,MFHDR,MFHDR_x
      use cplotr_mod, only : naxxtra,lrun_x
C
C  read the "new format" MF.PLN header
C
      integer MFLUNI                    ! lun of open MF.PLN file
      integer nfxt                      ! function count
      integer ier                       ! completion code returned.
      integer :: MFHDR1(iblksz),i_size
      integer, allocatable, dimension(:,:) :: itmp
C
C-----------------
C
      READ(MFLUNI,REC=1) (MFHDR1(II),II=1,IBLKSZ)
      IF(MFHDR1(1).NE.NFXT) THEN
         WRITE(lunzer(0),9003) MFHDR1(1),NFXT
 9003    FORMAT(' ?RPLOT DATA ERROR:  NO. OF FCNS IN MF FILE = ',I5/
     $        '  NUMBER INDICATED IN TF FILE = ',I5)
         IER=2
         return
      ELSE IF(MFHDR1(2).EQ.0) THEN
         WRITE(lunzer(0),9004)
 9004    FORMAT(' ?RPLOT MF DATA CONTAINS ZERO TIME PTS')
         IER=2
         return
      endif
      ISIZH=MFHDR1(1)+3+IBLKSZ/8
      IBLOKH=1+(ISIZH-1)/IBLKSZ
      ISIZH=IBLOKH*IBLKSZ
      if(lrun_x.eq.0) then
         if(.not.allocated(MFHDR)) then
            allocate(MFHDR(ISIZH))
            MFHDR=0
         else
            if(size(MFHDR).ne.ISIZH) then
               deallocate(MFHDR)
               allocate(MFHDR(ISIZH))
            endif
         endif
         MFHDR(1:IBLKSZ)=MFHDR1(1:IBLKSZ)
      else
         if(.not.allocated(MFHDR_x)) then
            allocate(MFHDR_x(ISIZH,naxxtra))
            MFHDR_x=0
         else
            i_size=size(MFHDR_x,dim=1)
            if(i_size.lt.ISIZH) then
               allocate(itmp(i_size,naxxtra))
               itmp=0
               itmp(1:i_size,1:naxxtra)=MFHDR_x(1:i_size,1:naxxtra)
               deallocate(MFHDR_x)
               allocate(MFHDR_x(ISIZH,naxxtra))
               MFHDR_x=0
               MFHDR_x(1:i_size,1:naxxtra)=itmp(1:i_size,1:naxxtra)
               deallocate(itmp)
            endif
         endif
         MFHDR_X(1:IBLKSZ,lrun_x)=MFHDR1(1:IBLKSZ)
      endif
C     
C     READ REST OF HEADER - DATA ADDRESSES
C     
      IF(IBLOKH.GT.1) THEN
         IPREV=IBLKSZ
         DO IREC=2,IBLOKH
            IWD1=IPREV+1
            IWD2=IPREV+IBLKSZ
            if(lrun_x.eq.0) then
               READ(MFLUNI,REC=IREC) (MFHDR(II),II=IWD1,IWD2)
            else
               READ(MFLUNI,REC=IREC) (MFHDR_x(II,lrun_x),II=IWD1,IWD2)
            endif
            IPREV=IWD2
         ENDDO
      ENDIF
C     
      return
      end
      
