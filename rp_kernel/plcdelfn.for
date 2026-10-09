      subroutine plcdelfn(zinput,ier)
C
      use datmgr_mod
      use cplotr_mod

      character*(*) zinput              ! name of function to delete
      integer ier                       ! status code on completion
C
C  ier=0 -- function was deleted successfully
C  ier=-1 -- function does not exist so need not be deleted
C  ier=1 -- function exists but deletion is not allowed.
C
      logical idlock
C
C	   .See if name is a user defined profile function.
C
      ier=0
      inamt=0
      inamr=0
 
      ! see if it is in memory... only user defined functions can
      ! be deleted and these only exist in memory

      Call DmdLoc(Zinput, Ind, Isize, Ipt)
 
      If (Ind .NE. 0 ) Then
         Ifxt0p = NFXT0+1
         infxt=nfxt
C               .See if name is user defined.
         Do 300 I=Infxt,Ifxt0p, -1
            If (Zinput .Eq. ABR(I)) Then
               Call ZerMsg(  '%PLCDEL:  Deleted user '
     1            // 'defined function "' // Zinput
     1            // '" from memory')
               call mgunref(zinput) ! remove mg references to fcn
               call plcdel0(i,zinput)  ! nfxt decremented
               inamr=i
               call mg_fixref(0,inamr) ! update MG prof refs to indices > inamr
            End If
 300     Continue
C
         if(inamr.gt.0) then
C
C	     .Delete from memory.
 
            idlock = no_delete
            no_delete=.FALSE.
 
            Call DMIDEL(Ind)
 
            no_delete = idlock

         endif
 
      End If                            ! Ind
 
      if(inamr.gt.0) return
C
C  profile function was not deleted; look for a scalar function
C
      call dmgfotx(2,ipt,ier)
C
      ifr=ifind_ordr(abr,iordrr,nfxt,zinput)
      ift=ifind_ordr(abt,iordrt,nft,zinput)
C
      if((ifr.eq.0).and.(ift.eq.0)) then
C
C  completely bogus function name
C
         ier=-1
         call zermsg(
     >      ' %plcdelfn:  non-existent function, cannot delete:  '//
     >      zinput)
         return
      endif
C
C  delete scalar function if eligible
C
      if(ift.gt.(nft0+2)) then
         inamt=ift
         call mgunref(abt(ift))         ! remove mg references to fcn
         call aordr_del(abt,iordrt,naxfot,nft,zinput) ! nft decremented
         call zermsg(' %plcdelfn:  scalar function deleted:  '//zinput)
         do j=ift,nft
            labelt(j)=labelt(j+1)
            unitst(j)=unitst(j+1)
            ttagt(j)=ttagt(j+1)
            iadold=ipt+j*ntt-1
            iadnew=iadold-ntt
            do it=1,ntt
               datbuf(iadnew+it)=datbuf(iadold+it)
            enddo
         enddo
         ttagt(nft+1)=0.0
         call mg_fixref(1,inamt) ! update MG scalar refs to indices > inamr
      endif
C
      if((inamr.eq.0).and.(inamt.eq.0)) then
C
C  function exists but is not user defined and cannot be deleted.
C
         ier=1
         ilz=len_trim(zinput)
         call zermsg(
     >      ' %plcdelfn:  '//zinput(1:ilz)//
     >      ' is a run defined function and'//
     >      ' cannot be deleted.')
      endif
C
      return
      end
C
C  actually delete the indicated profile function from memory -- labels, etc.
C
      subroutine plcdel0(i,zinput)

      use cplotr_mod
      character*(*) zinput
      integer i
C
      call aordr_del(abr,iordrr,naxfxt,nfxt,zinput) !nfxt decremented
C
      if(i.le.nfxt) then
         do j=i,nfxt
            labelr(j)=labelr(j+1)
            unitsr(j)=unitsr(j+1)
            itypr(j)=itypr(j+1)
            ttagr(j)=ttagr(j+1)
         enddo
         ttagr(nfxt+1)=0.0
      endif
C
      return
      end
