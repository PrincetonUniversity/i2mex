      subroutine rpdims(istype,irank,idims,zxabb,ier)
C
C  return dimensioning information for given function subtype code
C
C  subroutine "rplabel" returns subtype codes given the function name.
C
      use cplotr_mod

      integer istype                    ! input:  function subtype code
C
C         istype = -1   -- scalar f(t) function or multigraph
C         istype = 1,2,... -- various f(x,t) profile functions / multigraphs
C
      integer irank                     ! output:  rank or dimensionality
      integer idims(*)                  ! output:  sizes of dimensions
C
C         idims(1) = 1st (contiguous storage) dimension
C         idims(2) = 2nd dimension -- if irank.ge.2
C         ...etc...
C
C         idims(j) is not referenced if j.gt.irank.
C         nothing above irank=2 exists in TRANSP databases as of Aug 1999
C         (dmc).
C
      character*(*) zxabb(*)            ! for profiles:  X axis function(s)
C
C         the non-temporal coordinate is represented by the named function
C         this is a default specification; other X axis functions can be
C         chosen.
C
      integer ier                       ! exit code, 0=OK, 1= invalid istype.
C
      ier=1
      zxabb(1)=' '
C
      if(istype.eq.-1) then
         ier=0
         irank=1
         idims(1)=ntt
      else if((istype.gt.0).and.(istype.le.nxr)) then
         ier=0
         irank=2
         idims(1)=nzonex(istype)
         idims(2)=ntr
         zxabb(1)=xndabb(istype)
      endif
C
      if(ier.ne.0) then
         write(lunzer(0),1001) istype
 1001    format(' ?rpdims:  invalid rplot data type code:  ',i6)
      endif
C
      return
      end
