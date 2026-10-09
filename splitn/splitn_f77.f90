subroutine splitn_f77(zfile,ios)

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: zfile   ! filename of file to write
  integer, intent(out) :: ios          ! status code, 0=OK

  !  write f77-style fortran-readable TRANSP namelist (TR.ZDA file)
  !  only namelist elements with non-default values are written.  All
  !  namelists are written but some may be empty.

  !----------------------------
  integer i,inaml,ivar,jj,ilen,ilenp,ilname,isize,irank,inrank
  integer ilina0,ilina1,ilina,ilin,iw1,iw2
  integer indx(maxrank)
  character*10 ibuf
  character*40 indxbuf
  character*40 zname
  !----------------------------

  open(unit=lun,file=zfile,status='new',iostat=ios)
  if(ios.ne.0) return

  do inaml=1,nnamls
     ! loop over individual namelists

     cur_naml=all_namls(inaml)
     write(lun,'(1x,"&",a)') cur_naml(1:len_trim(cur_naml))

     do ivar = 1,nvars
        jj = var_order(ivar)  ! in alphabetic order

        ! loop over all namelist variables, selecting those from the
        ! current namelist

        if(cur_naml.eq.varlist(jj)%naml) then

           ! got one.  Compute size & initial indexing
           zname = varlist(jj)%name
           ilname = len_trim(zname)

           inrank=varlist(jj)%rank
           isize=1
           indx=0
           do irank=1,inrank
              isize= isize* &
                   (varlist(jj)%dims(2,irank)-varlist(jj)%dims(1,irank)+1)
              indx(irank)=varlist(jj)%dims(1,irank)
           enddo
           if(inrank.gt.0) indx(1)=indx(1)-1

           ilina0=varlist(jj)%nlinadr
           ilina1=ilina0+isize-1

           do ilina=ilina0,ilina1
              !  loop over elements (indexing; bounds check was done
              !  when namelist was read).

              if(inrank.gt.0) then
                 irank=1
                 do
                    indx(irank)=indx(irank)+1
                    if(indx(irank).le.varlist(jj)%dims(2,irank)) then
                       exit
                    else
                       indx(irank)=varlist(jj)%dims(1,irank)
                       irank=irank+1   ! carry to next dimension
                    endif
                 enddo
              endif
              !  loop over elements (line references)
              if(ilines(ilina).gt.0) then
                 ilin=ordl(ilines(ilina))
                 if(ivrange(1,ilina).gt.0) then
                    iw1=ivrange(1,ilina)
                    iw2=ivrange(2,ilina)
                    ! we have a non-default value to write...

                    if(inrank.eq.0) then
                       ! scalar write

                       write(lun,'(2x,a,"=",a)') &
                            zname(1:ilname),textnl(ilin)(iw1:iw2)

                    else
                       ! array element write: generate index;
                       ! no imbedded blanks.
                       indxbuf='('  ! start of index expression in char form

                       ilen=1
                       do irank=1,inrank
                          ibuf=' '
                          write(ibuf,'(I10)') indx(irank) ! right justified
                          do i=10,1,-1
                             if(ibuf(i:i).eq.' ') then
                                ilenp=ilen
                                ilen=ilenp+(10-i)
                                indxbuf(ilenp+1:ilen)=ibuf(i+1:10) ! copy in
                                exit
                             endif
                          enddo
                          ilen=ilen+1
                          if(irank.lt.inrank) then
                             indxbuf(ilen:ilen)=','
                          else
                             indxbuf(ilen:ilen)=')'
                          endif
                       enddo

                       write(lun,'(1x,a,a,"=",a)') &
                            zname(1:ilname),indxbuf(1:ilen), &
                            textnl(ilin)(iw1:iw2)
                    endif  ! scalar/array
                 endif  ! value in file
              endif  ! element reference in file

           enddo  ! loop over elements of namelist item

        endif  ! item is in current namelist
     enddo  ! loop over all  namelist items

     write(lun,'(1x,"/")')

  enddo  ! loop over all namelists

  close(unit=lun)

end subroutine splitn_f77
