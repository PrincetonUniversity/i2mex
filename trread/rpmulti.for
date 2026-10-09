      subroutine rpmulti(zname,istype,zlbl,zuns,infuns,isigns,zmembs,
     >   ier)
C
C  fetch information on a named multigraph
C
      use cplotr_mod

      character*(*) zname               ! multigraph name, input
C
      integer istype                    ! x axis type code, output
      character*(*) zlbl                ! multigraph label (C*32), output
      character*(*) zuns                ! multigraph units (C*16), output
      integer infuns                    ! number of function members, output
      integer isigns(*)                 ! sign code for each member, output
      character*(*) zmembs(*)           ! name of each member, output
C
      integer ier                       ! completion code, output, 0=OK
C
C  abnormal completion means:  multigraph name invalid.
C
      character*10 znami
C
C--------------------------------
C
      znami=zname
      call trcaps(znami)
C
      ilmemb=len(zmembs(1))
      ila=len(abt(1))
      ilact=0
C
      call rplabel(znami,zlbl,zuns,imulti,istype)
      if(imulti.ne.1) then
         call zermsg(' ?rpmulti:  not a multigraph name:  '//znami)
         ier=1                          ! not a multigraph
         infuns=0
         istype=0
         zlbl=' '
         zuns=' '
      else
         ier=0
         ia=ifind_ordr(abb,iordrb,nbal,znami)
         infuns=infb(ia)
         do i=1,infuns
            ifa=ifunb(i,ia)
            if(ifa.lt.0) then
               isigns(i)=-1
               ifa=-ifa
            else
               isigns(i)=1
            endif
            if(iintb(ia).eq.1) then
               zmembs(i)=abt(ifa)(1:min(ila,ilmemb))
               ilact=max(ilact,len_trim(abt(ifa)))
            else
               zmembs(i)=abr(ifa)(1:min(ila,ilmemb))
               ilact=max(ilact,len_trim(abr(ifa)))
            endif
         enddo
      endif
C
      if(ilact.gt.ilmemb) then
         call zermsg(
     >     ' ?rplist:  passed character array element width too small.')
         write(lunzer(0),*) '  passed width = ',ilmemb,
     >      ' current need = ',ilact,'  abt/abr width = ',ila
         ier=2
      endif
C
      return
      end
