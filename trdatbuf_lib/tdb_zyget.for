      real*8 function tdb_zyget(d,zid,i)
 
c   for 3d profile set indicated, get i'th y axis value
c
c   e.g. for multiple impurity model (sim) helper function
c   get stored charge state information...
 
      use trdatbuf_obj
      implicit NONE
      type (trdatbuf) :: d
      character*3, intent(in) :: zid  ! 3d profile set identifier
      integer,intent(in) :: i         ! get i'th element of y axis
 
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
C------------------
      IF(ZID.EQ.'SIM') THEN
        tdb_zyget = d%datbuf(d%lySIM+i-1)
      ENDIF
C
C=>TRDATGEN-
 
      return
      end
