      subroutine rplist(ztype,istype,alist,namax,ngot,ier)
C
C  return sorted list of RPLOT names of requested type.
C
      use cplotr_mod

      character*(*) ztype               ! input:  list type
C
C  choose ztype = "scalar" or "profile" or "multigraph"
C
      integer istype                    ! input:  list subtype control
C
C  if ztype = "scalar", istype is ignored.  for "multi" or "profile":
C
C  set istype = 0 to not restrict the list by subtype.
C  set istype = -1 for scalars only
C  set istype = +N for type N profiles only
C
      integer namax                     ! input:  size of output list
C
      character*(*) alist(namax)        ! output:  the list
C
      integer ngot                      ! actual number of items in list
C
      integer ier                       ! error code, 0=OK
C
C  ztype can be 'profile', 'scalar', or 'multigraph'
C   actually only the first 5 characters are tested with a case blind
C   compare.
C
C  ztype can contains a substring specification, using the syntax
C    <list-type>:<substring>
C  if this is done, only items which contain the substring are counted
C  as part of the list.  For example, 'profile:NIMP_' refers to profiles
C  whose names contain the substring "NIMP_".
C
C  note that because of RPLOT calculator action, these lists can grow
C  dynamically and may need to be refetched...
C
C  see subroutine rpnlist ... to get current list size
C  see subroutine rptype ... to get subtype data associated with names.
C
C  error codes (no message written):
C
C    ier=1  invalid "ztype" input
C    ier=2  width of "alist" elements less than what is required
C           (n.b. character*10 is standard)
C    ier=3  length of "alist" does not accomodate actual list.
C
C  local:
C
      character*5 ztest
      character*21 substr
      logical ck_rplist
C
C--------------------------------------
C
      ngot=0
C
      ilpass=len(alist(1))
      ila=len(abt(1))
      ilb=len(abb(1))
      ilact=0
C
      ier=0
      ilz=len(ztype)
      ztest=ztype(1:min(ilz,5))
      call trcaps(ztest)
C
      substr=' '
      icolon=index(ztype,':')
      if(icolon.gt.0) substr=ztype(icolon+1:min(icolon+10,ilz))
      call trcaps(substr)
C
      isize=0
      if(ztest.eq.'SCALA') then
         do i=1,nft
            j=iordrt(i)
            if(ck_rplist(abt(j),substr)) then
               isize=isize+1
               if(isize.gt.namax) go to 90
               alist(isize)=abt(j)(1:min(ilpass,ila))
               ilact=max(ilact,len_trim(abt(j)))
            endif
         enddo
      else if(ztest.eq.'PROFI') then
         do i=1,nfxt
            j=iordrr(i)
            if(ck_rplist(abr(j),substr)) then
               if(istype.eq.0) then
                  isize=isize+1
                  if(isize.gt.namax) go to 90
                  alist(isize)=abr(j)(1:min(ilpass,ila))
                  ilact=max(ilact,len_trim(abr(j)))
               else
                  if(itypr(j).eq.istype) then
                     isize=isize+1
                     if(isize.gt.namax) go to 90
                     alist(isize)=abr(j)(1:min(ilpass,ila))
                     ilact=max(ilact,len_trim(abr(j)))
                  endif
               endif
            endif
         enddo
      else if(ztest.eq.'MULTI') then
         do i=1,nbal
            j=iordrb(i)
            if(ck_rplist(abb(j),substr)) then
               if(istype.eq.0) then
                  isize=isize+1
                  if(isize.gt.namax) go to 90
                  alist(isize)=abb(j)(1:min(ilpass,ilb))
                  ilact=max(ilact,len_trim(abb(j)))
               else if(istype.eq.-1) then
                  if(iintb(j).eq.1) then
                     isize=isize+1
                     if(isize.gt.namax) go to 90
                     alist(isize)=abb(j)(1:min(ilpass,ilb))
                     ilact=max(ilact,len_trim(abb(j)))
                  endif
               else
                  if(iintb(j).eq.0) then
                     if1=iabs(ifunb(1,j))
                     if(itypr(if1).eq.istype) then
                        isize=isize+1
                        if(isize.gt.namax) go to 90
                        alist(isize)=abb(j)(1:min(ilpass,ilb))
                        ilact=max(ilact,len_trim(abb(j)))
                     endif
                  endif
               endif
            endif
         enddo
      else
         call zermsg(' ?rplist:  invalid list type:  '//ztype)
         ier=1                          ! unrecognized list type
      endif
C
      if(ilact.gt.ilpass) then
         call zermsg(
     >     ' ?rplist:  passed character array element width too small.')
         write(lunzer(0),*) '  passed width = ',ilpass,
     >      ' current need = ',ilact,'  abt/abr width = ',ila
         ier=2
      endif
C
      ngot=isize                        ! OK...
      go to 100
C
C  error -- list too long
C
 90   continue
      ier=3
C
 100  continue
      return
      end
