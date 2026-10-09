      subroutine idchek(zid,idnum,iclass,iuser,iadj)
C
      use cplotr_mod
C
C  determine what "zid" points at:
C
C  input:
      character*(*) zid                 ! putative function/mg id, input
C
C  output:
      integer idnum                     ! id# within class, or zero if undef.
C
      integer iclass                    ! id class
C               0:  undefined
C               1:  scalar function
C               2:  profile function
C               3:  multigraph
C
C              99:  ** illegal name **
C
      integer iuser                     ! user accessibility flag
C               0:  not user modifiable
C               1:  is user modifiable
C
      integer, intent(in) :: iadj       ! =0 --> "$" allowed,
C                                       ! =1 --> "$" NOT allowed, in zid
C
C----------------------------------------------------
C----------------------------------------------------
C
      idnum=0
      iclass=0
      iuser=0
C
      call idchek0(zid,iclass,iadj)
C
      if(iclass.eq.99) return
C
      idnum=ifind_ordr(abt,iordrt,nft,zid)
      if(idnum.gt.0) then
         iclass=1
         if(idnum.gt.(nft0+2)) iuser=1
         go to 1000
      endif
C
      idnum=ifind_ordr(abr,iordrr,nfxt,zid)
      if(idnum.gt.0) then
         iclass=2
         if(idnum.gt.nfxt0) iuser=1
         go to 1000
      endif
C
      idnum=ifind_ordr(abb,iordrb,nbal,zid)
      if(idnum.gt.0) then
         iclass=3
         go to 1000
      endif
C
 1000 continue
      return
      end
c--------------------------------------------------------------
      subroutine idchek0(zid,iclass,iadj)
c
c  check one id for syntax
c
      character*(*) zid                 ! identifier
      integer iclass                    ! set to 99 if error occurs
c                                       ! set to 0 if OK
      integer, intent(in) :: iadj       ! drop last <iadj> chars from zlegal
c
c     if iadj.gt.0 also enforce length limit=10
c
c---------------
C  local:
c
      integer ic,indx,ili,ilz
c
      character*38 zlegal
C
      DATA ZLEGAL/'1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ_$'/
c
c---------------------------
c
      ilz=len_trim(zlegal)-iadj
c
      iclass=0
c
      ili=len_trim(zid)
c
      if(iadj.gt.0) then
         if(ili.gt.10) then
            call zermsg(
     >         ' ?idchek:  identifier too long:  '//zid)
            iclass=99
            return
         endif
      endif
c
      do ic=1,ili
         indx=index(zlegal(1:ilz),zid(ic:ic))
         if(indx.le.0) then
            call zermsg(
     >         ' ?idchek:  illegal character in identifier:  '//zid)
            iclass=99
         else if((ic.eq.1).and.(indx.le.10)) then
            call zermsg(
     >         ' ?idchek:  1st character illegal:  a digit:  '//zid)
            iclass=99
         endif
      enddo
c
      return
      end

      subroutine idchek_setadj(nam,iadj)

      !  if name is of form <name>$<runid> set iadj=0
      !  otherwise set iadj=1

      implicit NONE

      character*(*), intent(in) :: nam
      integer, intent(out) :: iadj

      !-----------------------------------------
      integer :: idollr,ilen,ic
      !-----------------------------------------

      iadj=1
      idollr=0

      ilen = len(trim(nam))

      if((nam(1:1).ne.'$').and.(nam(ilen:ilen).ne.'$')) then
         do ic=1,ilen
            if(nam(ic:ic).eq.'$') idollr=idollr+1
         enddo
      endif

      if(idollr.eq.1) iadj=0       ! allow <name>$<runid> identifer on LHS

      return
      end
