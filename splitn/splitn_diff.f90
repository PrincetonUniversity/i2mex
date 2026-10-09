subroutine splitn_diff(ilun,idiff,ierr)

  ! compare contents of current and previous namelist
  ! write warning messages on unit ILUN (if it is non-zero)
  ! number of variables with differences found returned in IDIFF
  ! number of changed "non-TRDAT-steerable" variables retuned in IERR
  ! IERR = -1 returned if there is no previous namelist

  ! NOTE: this compares the main namelist values only.  Update namelists
  ! are ignored: all variables in update namelists are TRDAT-steerable.

  !  Mods
  !    06Apr2009    jim.conboy@jet.uk    BUGFIX
  !                 idiff is inout
  !                 skip 'lookbehind' if no differences ; label do/enddo 
  !                 ToDo :
  !                        Another arg to control output level
  !                        rank 2 variables - subscripts ??
  !
  !    29Mar2009    jim.conboy@jet.uk [ build date 'Jan 24 04:30:36', JET vn 44
  !                 New - print all differences, if idiff /= 0 on entry
  !                 ( See splitn_comp )
  !
  !    29Mar2009    jim.conboy@jet.uk [ build date 'Jan 24 04:30:36', JET vn 43]
  !                 BUGfix - test only 1st char of type, for character variable
  !

  use splitn_module
  implicit NONE

  integer, intent(in)    :: ILUN  ! i/o channel for warnings, or 0 for "quiet"
  integer, intent(inout) :: IDIFF ! # of changed variables found
  integer, intent(out)   :: IERR  ! # of non-steerable changed variables or...
  !                              -1 if there is no previous namelist to compare

  !-----------------------------------------------
  ! local
  integer, parameter            :: nszmax = 20& ! max diffs /variable to print
                                  ,ncbuf  = 20

  integer                       :: ivar,jj,isize,i,ii,iadd,ilog

  integer                       :: isame        ! last element with same difference
!
  logical                       :: lprt_all     ! .T to print all differences
  !
  character(len=2)              :: csteer       ! flag 'non steerable' parameter
  character(len=8)              :: cix          ! prints range of var index
  !
  character(len=ncbuf),&
           dimension(2,nszmax)  ::  cdbuf       !  list of diffs
  integer, dimension(nszmax)    ::  ixszdf &    !  index of different elements
                                   ,ixszdf2     !  .. ptr for consec diff elems


  !-----------------------------------------------
  ! verify presence of prior namelist for comparison...

  lprt_all = (idiff.ne.0)

  idiff = -1
  ierr = -1
  if(.not.previous) then
     if(ilun.gt.0) then
        write(ilun,*) ' %splitn_diff -- no prior namelist => ierr = -1'
     endif
     return
  endif

  !-----------------------------------------------
  ! prior namelist exists

  idiff = 0
  ierr = 0

  lvar: do ivar = 1,nvars
     jj = var_order(ivar)

     isize=1
     do i=1,varlist(jj)%rank
        isize=isize*(varlist(jj)%dims(2,i)-varlist(jj)%dims(1,i)+1)
     enddo
     ii = varlist(jj)%addr

     iadd = 0  ! count differences, this variable...

     lelem: do i=ii,ii+isize-1

        if(varlist(jj)%type.eq.'I')         then
           if(intbuf(i).ne.intbuf_p(i))  then
              iadd = iadd + 1
              ixszdf(iadd) = 1+i-ii
              write(cdbuf(1,iadd),'(i20)')  intbuf(i)
              write(cdbuf(2,iadd),'(i20)')  intbuf_p(i)
           endif

        else if(varlist(jj)%type.eq.'L')     then
           ilog = 0
           if(logbuf(i)) ilog = ilog + 1
           if(logbuf_p(i)) ilog = ilog + 1
           if(ilog.eq.1)                 then  ! 0 means both F, 2 means both T
              iadd = iadd + 1
              ixszdf(iadd) = 1+i-ii
              write(cdbuf(1,iadd),'(19x,l1)')  logbuf(i)
              write(cdbuf(2,iadd),'(19x,l1)')  logbuf_p(i)
           endif

        else if(varlist(jj)%type.eq.'R')    then
           if(rbuf(i).ne.rbuf_p(i))      then
              iadd = iadd + 1
              ixszdf(iadd) = 1+i-ii
              write(cdbuf(1,iadd),'(f20.6)')  rbuf(i)
              write(cdbuf(2,iadd),'(f20.6)')  rbuf_p(i)
           endif

        else if(varlist(jj)%type.eq.'D')    then
           if(dbuf(i).ne.dbuf_p(i))      then
              iadd = iadd + 1
              ixszdf(iadd) = 1+i-ii
              write(cdbuf(1,iadd),'(1PG20.6)')  dbuf(i)
              write(cdbuf(2,iadd),'(1PG20.6)')  dbuf_p(i)
           endif

        else if(varlist(jj)%type(1:1).eq.'C')       then
!d-        print *, trim(varlist(jj)%name), &
!d-             trim(chbuf(i)), trim(chbuf_p(i))
           if(chbuf(i).ne.chbuf_p(i))     then
              iadd = iadd + 1
              ixszdf(iadd)  = 1+i-ii
              cdbuf(1,iadd)(:ncbuf) = chbuf(i)
              cdbuf(2,iadd)(:ncbuf) = chbuf_p(i)
           endif
!
        endif
        if ( iadd == 0 )                                cycle lelem
!       lookbehind to avoid printing 'same' differences

        if(iadd.gt.0) then
           ixszdf2(iadd) = ixszdf(iadd)
           if( iadd == 1 )                                 then
              isame = 1
           else 
              if(      ixszdf(iadd-1) + 1 == 1+i-ii     &
                   .and. cdbuf(1,iadd) == cdbuf(1,iadd-1) &
                   .and. cdbuf(2,iadd) == cdbuf(2,iadd-1) )  then
                 ixszdf2(isame) = ixszdf(iadd)
              else 
                 isame = iadd
              endif
           endif
        endif
     enddo lelem

     if(iadd.gt.0) then
        ! difference detected...

        idiff = idiff + 1
        if( lprt_all )      then
           if( varlist(jj)%steerable.ne.0 ) then
              csteer = '  '
           else
              csteer = '* '
              ierr = ierr + 1
           endif
           cix = '        '
           if( isize > 1 ) &
                write(cix,'(i3,a2,i3)') ixszdf(1), '..', ixszdf2(1)
           if( ixszdf2(1) == ixszdf(1) )  then
              cix(4:8) = '     '
              isame = 1
           else
              isame = ixszdf2(1)
           endif
!
           write(ilun,'(4x,4a,t35,a,2x,a,t70,a)')  &
                csteer, varlist(jj)%type,'changed ', &
                trim(varlist(jj)%name), cix, cdbuf(1,1),cdbuf(2,1)
           if(iadd .gt. 1 )                           then
              ix: do i=2,iadd
                 if (i <= isame )               cycle ix
                 write(cix,'(i3,a2,i3)') ixszdf(i), '..', ixszdf2(i)
                 if( ixszdf2(i) == ixszdf(i) ) cix(4:8) = '     '
                 isame = ixszdf2(i)
                 write(ilun,'(t35,a,2x,a,t70,a)') &
                      cix, cdbuf(1,i),cdbuf(2,i)
              enddo ix
              write(ilun,*) ' '
           endif
        else
           if( varlist(jj)%steerable.eq.0) then
              ierr = ierr + 1
              write(6,*) &
                   ' %splitn_diff: value change in non-TRDAT-steerable quantity: ',trim(varlist(jj)%name)
           endif
        endif

     endif

  enddo lvar

end subroutine splitn_diff
