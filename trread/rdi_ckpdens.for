c-------------------------------------------------------------------------
c
      subroutine rdi_ckpdens(infuns,zfuns,ias,izs,ierr)
c
c  check contents of TRANSP run's PDENS multigraph against known names...
c
      implicit NONE
c
c  input:
      integer infuns                 ! number of functions
      character*(*) zfuns(infuns)    ! PDENS multigraph contents
c
c  output:
      integer ias(infuns)            ! A of species (or 0 to ignore)
      integer izs(infuns)            ! Z of species (or 0 to ignore)
c
      integer ierr                   ! 0=OK; else index to unrecognized item
c
c---------------------------------
c
      integer, parameter :: nignore=6
      integer, parameter :: naccept=6
c
      character*10 ztest,zignore(nignore),zaccept(naccept)
      integer iaa(naccept),izz(naccept)
c
      integer i,j,k
c
      data zignore/'NE','NIMP','BDENS','NFI','NALPHA','NMINI'/
      data zaccept/'NH','ND','NT','NHE3','NHE4','NLITH'/
      data iaa/1,2,3,3,4,6/
      data izz/1,1,1,2,2,3/
c
c---------------------------------
c
      ias=0
      izs=0
      ierr=0
c
      do i=1,infuns
         ztest=zfuns(i)
         call uupper(ztest)
         k=0
         do j=1,naccept
            if(ztest.eq.zaccept(j)) then
               k=j
               ias(i)=iaa(j)
               izs(i)=izz(j)
               exit
            endif
         enddo
         if(k.eq.0) then
            do j=1,nignore
               if(ztest.eq.zignore(j)) then
                  k=j
                  exit
               endif
            enddo
         endif
         if(k.eq.0) then
            ierr=k
            exit
         endif
      enddo
c
      return
      end
