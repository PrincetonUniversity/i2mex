      subroutine trintrp(inxnew,intimes,ztime,zfdat,zxdat,ixref,inx,
     >   zdata)
c
      use datmgr_mod
      use cplotr_mod
c
c  given zfdat(inxnew,intimes) and zxdat(inxnew,intimes) interpolate
c  the data onto the current run's time-space grid.  argument `ixref'
c  indicates which space grid to use
c
      integer, intent(in) :: inxnew     ! no. of x index pts in input data
      integer, intent(in) :: intimes    ! no. of time pts in input data
      real, intent(in) ::  ztime(intimes)        ! input data timebase
      real, intent(in) ::  zfdat(inxnew,intimes) ! input data
      real, intent(in) ::  zxdat(inxnew,intimes) ! x data (time dependent)
c
      integer, intent(in) :: ixref      ! corresponding x data, current run
      integer, intent(in) :: inx        ! no. of x pts, current run
c
      real, intent(out) :: zdata(inx,ntr)  ! interpolated data (out)
c
      EXTERNAL XIDENT                   ! FUNCTION PASSED TO XINTER
c
c  local:
c
      real :: ztmin,ztmax
      real, dimension(:,:), allocatable :: zfdatt,zxdatt
c-----------------------------------------------
c
      allocate(zfdatt(inxnew,ntr),zxdatt(inxnew,ntr))
C
C  TIME INTERPOLATE THE DATA & x axis
C
      ztmin = minval(ztime(1:intimes))
      ztmax = maxval(ztime(1:intimes))
      DO 1100 IT=1,NTR
         if(nlxtrap0.and.
     1        ((time3(it).lt.ztmin).or.(time3(it).gt.ztmax))) then
            zfdatt(1:inxnew,it)=0.0
            if(time3(it).lt.ztmin) then
               zxdatt(1:inxnew,it)=zxdat(1:inxnew,1)
            else
               zxdatt(1:inxnew,it)=zxdat(1:inxnew,intimes)
            endif
         else
            CALL XINTER(XIDENT,TIME3(IT),ZTIME,intimes,
     1           IT0,IT0P1,ZTI,ZTIC,IEX)
C
            DO 190 IX=1,INXNEW
               zfdatt(ix,it)=ztic*zfdat(ix,it0)+zti*zfdat(ix,it0p1)
               zxdatt(ix,it)=ztic*zxdat(ix,it0)+zti*zxdat(ix,it0p1)
 190        CONTINUE
C
         endif
 1100 CONTINUE
C
C  SPACE INTERPOLATE THE DATA
C
      CALL DMGFXT(ixref,IND)
      IPX=LOCD(IND)
C
      do it=1,ntr
         do ix=1,inx
            zxval=datbuf(ipx-1+inx*(it-1)+ix)
            call xinter(XIDENT,zxval,zxdatt(1,it),inxnew,
     >         ix0,ix0p1,zxi,zxic,iex)
 
            zdata(ix,it)=
     >         zxic*zfdatt(ix0,it)+zxi*zfdatt(ix0p1,it)
C
         enddo
      enddo
C
      deallocate(zfdatt,zxdatt)
      return
      end
