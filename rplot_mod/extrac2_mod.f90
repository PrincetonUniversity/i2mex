!*************** START FILE EXTRAC.BLK ; GROUP EXTRAC *************
!
! **EXTRAC2.BLK** -- same as EXTRAC.BLK but with smaller parameters
!    used to give EXTRAC like capabilities to RPLOT but not blow out
!    its memory requirements.  DMC 23 Jul 1996
!
!  COMMON BLOCK AND PARAMETERS FOR UFILES EXTRAC UTILITY
!
module extrac2_mod
  implicit none

  integer, parameter :: MAXIND = 42000
  integer, parameter :: MAXTOT = 420000
  integer, parameter :: MAXSC = 500
  !
  !  DATA BUFFERS AND LABELS
  REAL X(MAXIND),Y(MAXIND),ZWK(MAXIND),ZWK1(MAXIND), &
       F(MAXTOT),SCV(MAXSC)
  REAL XB(MAXIND),YB(MAXIND),FB(MAXTOT)
  LOGICAL LX(MAXIND),LY(MAXIND)
  REAL XTARG(MAXIND),YTARG(MAXIND),XTRAP(3,2),YTRAP(3,2)
  INTEGER IORDSC(MAXSC)
  integer :: ILUNI,ILUNO,ILUNT,ILUNC,IPROC
  integer :: ISHOT,NUMSC,NXACT,NYACT,ICOMPR
  integer :: ILUNX,ILUNY,NXTARG,NYTARG
  integer :: LXTRAP,LYTRAP
  !
  CHARACTER*10 ZSDATE,XLAB(3),YLAB(3),FLAB(3),ZTITLE(3)
  CHARACTER*4 IDEV
  CHARACTER*10 LYABL(5),SCLAB(3,MAXSC)
  CHARACTER*3 CMPSTA(2)

end module extrac2_mod
