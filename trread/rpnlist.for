C-----------------
C  rpnlist -- return size of an RPLOT/TRANSP database symbol list
C
      subroutine rpnlist(ztype,istype,isize)
C
      use cplotr_mod
C
      character*(*) ztype               ! input:  list type
C
C  choose ztype = "scalar" or "profile" or "multigraph"
C
      integer istype                    ! input:  list subtype
C
C  if ztype = "scalar", istype is ignored.  for "multi" or "profile":
C
C  set istype = 0 to not restrict the list by subtype.
C  set istype = -1 for scalars only
C  set istype = +N for type N profiles only
C
      integer isize                     ! output:  list size
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
C  see subroutine rplist ... to get the actual names.
C  see subroutine rptype ... to get subtype data associated with names.
C
C  if ztype is invalid, isize=0 is returned without further explanation.
C
C  local:
C
      character*5 ztest
      character*10 substr
      logical ck_rplist
C
C---------------------------
C
      isize=0
      ilz=len(ztype)
      ztest=ztype(1:min(5,ilz))
      call trcaps(ztest)
C
      substr=' '
      icolon=index(ztype,':')
      if(icolon.gt.0) substr=ztype(icolon+1:min(icolon+10,ilz))
      call trcaps(substr)
C
      if(ztest.eq.'SCALA') then
         do i=1,nft
            if(ck_rplist(abt(i),substr)) isize=isize+1
         enddo
      else if(ztest.eq.'PROFI') then
         do i=1,nfxt
            if(ck_rplist(abr(i),substr)) then
               if(istype.eq.0) then
                  isize=isize+1
               else
                  if(itypr(i).eq.istype) isize=isize+1
               endif
            endif
         enddo
      else if(ztest.eq.'MULTI') then
         do i=1,nbal
            if(ck_rplist(abb(i),substr)) then
               if(istype.eq.0) then
                  isize=isize+1
               else if(istype.eq.-1) then
                  if(iintb(i).eq.1) isize=isize+1
               else
                  if(iintb(i).ne.1) then
                     if1=iabs(ifunb(1,i))
                     if(itypr(if1).eq.istype) isize=isize+1
                  endif
               endif
            endif
         enddo
      else
         call zermsg(' ?rplist:  invalid list type:  '//ztype)
      endif
C
      return
      end
