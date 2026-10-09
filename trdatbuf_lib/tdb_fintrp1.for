C------------------------------------------------------------
C  interpolation for "standard" trdatbuf f(t) items
C
C---------------
C  Background information:  a code generator is used to create
C  data channels for TRANSP/trdatbuf-- this makes it easy to add 
C  new data channels e.g. as new measurements are invented on
C  tokamak experiments.
C
C  The input file for the code generator, trdatgen, is a file 
C  called trdatgen.spec.  This file defines profile channels,
C  scalar channels, and "special handling" channels.  An f(t)
C  item is a scalar channel.
C
C  Each channel has a 3 character name or "trigraph".  For 
C  example, the name for plasma current is "CUR".  The header
C  entry for CUR in trdatgen.spec looks like:
C
C      *CUR(TIM)        { 'plasma current',     'amps' }
C
C  Another example:
C      *VSF(TIM)        { 'surface voltage',    'volts' }
C
C  The generated code uses the trigraphs CUR and VSF to generate various
C  named component connected with the lookup and time interpolation of 
C  plasma current and surface voltage, respectively.
C
C  Any header entry in trdatgen.spec that starts with
C
C      *<trigraph>(TIM)
C
C  specifies a scalar channel, the data for which can be interpolated
C  using the software described here, provided that the data exists
C  for a given trdatbuf object.
C
C  For any given instance of a trdatbuf object, all scalar channel data
C  is defined over a single shared time base, i.e. a 1d strictly ascending
C  finite sequence of time values, in seconds.  The time base and the 
C  associated scalar channel data are not necessarily evenly spaced in
C  time.
C
C---------------
C  What to do in your code (what follows assumes you have instantiated
C  a trdatbuf object "d" and have read its contents in from a trdatbuf
C  file.  Your code would then have access to the declarations
C      use trdatbuf_module
C      type (trdatbuf) :: d
C  but of course a different name could be chosen for the trdatbuf 
C  object.  Then:
C
C      logical :: cur_exist
C      integer :: cur_addr
C
C  (do this maybe once at the start or restart of a run-- but not every 
C  timestep, since the results won't change and the lookup could be slow):
C      cur_exist = TDB_PRESENT1(d,"CUR",cur_addr)
C
C  and later when an interpolation to a particular time is needed:
C
C      real*8 :: ztime     ! desired time (in) (seconds).
C      integer :: it1,it2  ! time bin indices
C      real*8 :: zf1,zf2   ! linear interpolation factors
C      real*8 :: cur_now   ! the CUR data at the desired time (out)
C      real*8 :: vsf_now   ! the VSF data at the desired time (out)
C
C  call TDB_LOOKUP1(d,ztime,it1,it2,zf1,zf2)         ! time bin lookup
C  cur_now = TDB_FINTRP1(d,cur_addr,it1,it2,zf1,zf2) ! get CUR at ztime
C  vsf_now = TDB_FINTRP1(d,vsf_addr,it1,it2,zf1,zf2) ! get VSF at ztime.
C
C  Note that the TDB_LOOKUP1 call is done once; it sets it1,it2,zf1,zf2
C  which are then reused in multiple subsequent linear interpolations,
C  avoiding repeated calculations of the lookup in the scalar channels'
C  time base inside the trdatbuf object.  This takes advantage of the 
C  fact that all scalar channels share a common time base, for efficiency.
C
C  The interfaces to the logical function TDB_PRESENT1, the subroutine
C  TDB_LOOKUP1, and the real*8 function TDB_FINTRP1, are all defined by
C  use association via the module, i.e. "use trdatbuf_module".
C
C------------------------------------------------------------
C  TDB_LOOKUP1
C
C  PREPARE TO INTERPOLATE ON A FUNCTION OF TIME
C
C  ZTIME-- TIME TO INTERPOLATE TO
C  ILOC-- LOCATION OF FCN IN BUFFER
c
c  interpolation formula for data stored at location ILOC
c  will be ...
c
c    result = d%datbuf(iloc+it1)*zf1 + d%datbuf(iloc+it2)*zf2
C
c  ** see REAL*8 function TDB_FINTRP1, below...
c
c------------
C
      SUBROUTINE TDB_LOOKUP1(d,ztime,it1,it2,zf1,zf2)
C
      use trdatbuf_obj
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C  no use of SAVEd variables; repetitive calls at the same time will
C  cause look up calculation to be repeated.
C
      type (trdatbuf) :: d
      real*8, intent(in) :: ztime
      integer, intent(out) :: it1,it2  ! time offsets btw which ZTIME lies
      real*8, intent(out) :: zf1,zf2   ! linear interpolation factors
C
      integer :: int,ilt,ilt1,ilt2
      real*8 :: ztinit,zftime,zfac
      real*8, parameter :: ONE = 1.0d0
C
C-------------------------------------
C
      ilt=d%ltime1
      int=d%ntime1
      ztinit = d%datbuf(ilt)
      zftime = d%datbuf(ilt+int-1)
C
C check bounds
C
      if(ztime.le.ztinit) then
         it1=0
         it2=1
         zf1=1
         zf2=0
         return
      endif
      if(ztime.ge.zftime) then
         it1=int-2
         it2=int-1
         zf1=0
         zf2=1
         return
      endif
C
C initial lookup estimate
C
      zfac=(ztime-ztinit)/(zftime-ztinit)
      it1=zfac*int
      it1=max(0,min(int-2,it1))
      ilt1=ilt+it1

      do
         ilt2=ilt1+1
         if(d%datbuf(ilt1).gt.ztime) then
            ilt1=ilt1-1
            cycle
         else if(d%datbuf(ilt2).lt.ztime) then
            ilt1=ilt1+1
            cycle
         else
            exit
         endif
      enddo
C
C time bin found:
      it1=ilt1-ilt
      it2=it1+1
      ztinit=d%datbuf(ilt1)
      zftime=d%datbuf(ilt2)
C
C displacement within the bin:
      zf2=(ztime-ztinit)/(zftime-ztinit)
      zf1=ONE-zf2

      return
      end
C-------------------
      real*8 function TDB_FINTRP1(d,iloc,it1,it2,zf1,zf2)
C
      use trdatbuf_obj
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
      type (trdatbuf) :: d
      integer, intent(in) :: iloc     ! data location
      integer, intent(in) :: it1,it2  ! time offsets btw which ZTIME lies
      real*8, intent(in) :: zf1,zf2   ! linear interpolation factors
C
      integer lunmsg_tdb
C
      if(iloc.le.0) then
         write(lunmsg_tdb(0),*) ' ?TDB_FINTRP1: invalid address: ',iloc
         TDB_FINTRP1 = 0
      else
         TDB_FINTRP1 = d%datbuf(iloc+it1)*zf1 + d%datbuf(iloc+it2)*zf2
      endif
C
      return
      end
