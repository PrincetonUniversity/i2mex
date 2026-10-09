      subroutine splitn_lhs(zlhs,imatch,indx,imaxrank,ierr)
C
C  parse LHS of namelist line
C  Find variable name, find indexing.
C
C   <namelist-name>[(indx1,indx2,...)] = RHS...
C
      use splitn_module
      implicit NONE
C
      character*(*),intent(in) :: zlhs  ! input line's LHS
C                         **NO leading or trailing blanks**
C
      integer, intent(out) :: imatch    ! output variable match index
C                ....0 if no match
C
      integer, intent(in) :: imaxrank    ! max rank of namelist variables
      integer, intent(out) :: indx(imaxrank) ! output start index (arrays only)
      integer, intent(out) :: ierr      ! status code: 0=OK
C
C--------------------
      integer ilen,i,iw1,iw2,ilparen,irparen,ilb
      integer ibefore,jj,irank,inumdim,ic,iprev,id,ir
      integer idcod_splitn
C
      character*50 ztest
      character*20 berr
C
      integer, parameter :: maxlen = 32  ! change to 32 to allow 32 char names
C
C--------------------
C
C  1.  find var. name; 1st non-blank...
C
      imatch=0
      do i=1,imaxrank
         indx(i)=0
      enddo
C
      ierr=0
C
      ilen=len(zlhs)
C
      iw1=1
C
      ilparen=index(zlhs,'(')

      if(ilparen.eq.0) then
         iw2=ilen
         irparen=0
      else
         irparen=ilen
         do i=ilparen-1,1,-1
            if((zlhs(i:i).ne.' ').and.(zlhs(i:i).ne.char(9))) then
               iw2=i
               exit
            endif
         enddo
         if(zlhs(irparen:irparen).ne.')') then
            write(6,*) ' ?splitn_lhs: closing parenthesis missing or ',
     >         'misplaced.'
            ierr=1
         endif
      endif
C
      ztest=zlhs(iw1:iw2)
      call uupper(ztest)
      if((iw2-iw1+1).gt.maxlen) then
         write(6,*) ' ?splitn_lhs: name too long or invalid: "',
     >      zlhs(iw1:iw2),'"'
         imatch=0
      else
         call iorder(ztest,ibefore,imatch)
      endif
      if(imatch.eq.0) then
         if(.NOT.is_deleted(ztest)) then
            write(6,*) ' ?splitn_lhs: unrecognized name in namelist: "',
     >           zlhs(iw1:iw2),'"'
            jj=var_order(ibefore)
            ierr=1
         else
            write(6,*) ' %splitn_lhs: commenting out deleted item: "',
     >           zlhs(iw1:iw2),'"'
         endif
      endif
C
      if(ierr.ne.0) return
      if(imatch.eq.0) return
C
      imatch=ibefore                    ! this is the item order#
C
      jj=var_order(imatch)
      irank=varlist(jj)%rank
C
C  steerability check
C
      if(kupdate.gt.0) then
         if(varlist(jj)%steerable.le.1) then
            write(6,*) 
     >           ' ?splitn_lhs: variable in update block not '//
     >           'updatable: "',zlhs(iw1:iw2),'"'
            write(6,'(a,1pe13.6)') '  ~UPDATE_TIME (sec) is: ',
     >           tup(kupdate)
            write(6,*) '    update block index: ',kupdate
            ierr=1
         endif
      endif
      if(ierr.ne.0) return
C
      do i=1,irank
         indx(i)=varlist(jj)%dims(1,i)
      enddo
C
C  OK have match & initial address.  Look for = sign or subscripts
C
      if(ilparen.eq.0) return  ! initial address OK
C
C  subscripts inside zlhs(ilparen:irparen)
C
      inumdim=0
      iprev=ilparen
      ic=ilparen
C
 60   continue
      ic=ic+1
      if(zlhs(ic:ic).eq.',') go to 70
      if(zlhs(ic:ic).eq.')') go to 70
      go to 60
C
 70   continue
      inumdim=inumdim+1
      if(inumdim.gt.irank) go to 98
      if(iprev+1.gt.ic-1) go to 99
      id=idcod_splitn(zlhs(iprev+1:ic-1),ierr)
      if(ierr.ne.0) go to 99
      indx(inumdim)=id
      if(zlhs(ic:ic).eq.')') then
         if(inumdim.lt.irank) go to 98
      else
         iprev=ic
         go to 60                       ! next dim
      endif
C
C  check subscript bounds
C
      do ir=1,irank
         berr='out of bounds low'
         if(indx(ir).lt.varlist(jj)%dims(1,ir)) go to 97
         berr='out of bounds high'
         if(indx(ir).gt.varlist(jj)%dims(2,ir)) go to 97
      enddo
C
C  successful decode of address.
C
 80   continue
C
      return
C
C  errors
C
 97   continue
      ilb=len_trim(berr)
      write(6,9701) ir,berr(1:ilb)
 9701 format(' ?splitn_lhs -- subscript bounds check error:'/
     >'  subscript #',i1,' is ',a,'.')
      ierr=1
      return
C
 98   continue
      write(6,'('' ?splitn_lhs -- wrong no. of subscripts.'')')
      ierr=1
      return
C
 99   continue
      write(6,'('' ?splitn_lhs -- decode error in subscript.'')')
      ierr=1
      return
C
      end
