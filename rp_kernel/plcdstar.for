      subroutine plcdstar
C
C  delete all user profiles and scalars
C
      use datmgr_mod
      use cplotr_mod
C
      Character*21 Zinput
C
      logical idlock
C
C--------------------------------------------------------
C  check scalars
C
      call dmgfotx(2,ipt,ier)
C
C  delete the profiles...
C
      Ifxt0p = NFXT0+1
      If (Ifxt0p .LE. NFXT) Then
C               .See if name is user defined.
         idlock = no_delete
         no_delete = .FALSE.            ! enable DATMGR table entry deletions
         Do 300 I=Nfxt,Ifxt0p,-1
            Zinput = ABR(I)
C	             .See if name is a user defined function.
 
            Call DmdLoc(Zinput, Ind, Isize, Ipt)
 
            If (Ind .NE. 0 ) Then
               call mgunref(zinput)
               Call Zermsg(' %PLCDSTAR: Deleting user defined'
     1            // ' profile:' // Zinput)
               Call DMIDEL(Ind)
               call plcdel0(i,zinput)
            Else
C   		       .Error name not found
               Call Zermsg(' ?DlcDelAl: Code error, User fcn ' //
     1            Zinput // ' not found in memory')
               Call AbortT
            End If                      ! Ind
 300     Continue
         no_delete=idlock
 
      End If                            ! Ifxtp
C
C  now delete all user scalars
C
      inftp=nft
      nft=nft0+2
      if(nft.lt.inftp) then
         do i=nft+1,inftp
            call mgunref(abt(i))
         enddo
         call zermsg(' %PLCDSTAR:  all user defined scalars deleted.')
         call aordr(iordrt,abt,nft)     ! recompute name ordering
      endif
C
      return
      end
