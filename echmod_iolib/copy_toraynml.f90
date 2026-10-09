!----------------------------------------------------------------------
!  COPY_TORAYNML -- subroutine to copy Toray template namelist
!
!  calling arguments only, trcom & TRDATA not used...
!
subroutine copy_toraynml(ctok,runid,toray_template, &
     kmdsplus,mds_server,mds_tree, &
     ztarget, &
     nonlin,ier)

  implicit none
  
  !-----------------------------------------------------------------------
  !  subroutine arguments:
  
  character(len=*), intent(in) :: ctok   ! tokamak ID
  character(len=*), intent(in) :: runid  ! TRANSP run ID
  character(len=*), intent(in) :: toray_template  ! template name
  
  integer, intent(inout) :: kmdsplus  ! flag, =1 or 2 for MDS+ access
  !  if =2, assume MDS+ tree is open;
  !    returnED VALUE is set to 0 if MDSplus access failed.
  
  character(len=*), intent(in) :: mds_server  ! MDS+ server
  character(len=*), intent(in) :: mds_tree    ! MDS+ tree name
  
  character(len=*), intent(in) :: ztarget   ! where to copy the file
  
  integer, intent(in) :: nonlin       ! I/O unit for messages
  
  integer, intent(out) :: ier ! status return code, 0=OK
  
  character(len=17) :: subr = 'COPY_TORAYNML: '
  !
  ! MDSplus Fortran interface
  !
  integer :: Mds_Value,idescr_long,idescr_cstringarr,idescr_cstring
  integer :: isize,status,idims(2)
  
  !  buffer for text
  
  integer :: iline,inlines
  character(len=150), dimension(:), allocatable :: nltext
  
  !-----------------------------------------------------------------------
  !  local:
  
  integer :: istat,iltok,ilnurd,ilun
  
  character(len=150) :: zval1,zval2,zval3,zval4
  character(len=200) :: ztmppath,ztag,zprefix,zufilespath
  character(len=300) :: zgentmplfile,zcmd
  character(len=3) ::   zcp

  !-----------------------------------------------------------------------
  !
  ier=0
  !
  !----------

  zcp='cp'

  iltok=ilnurd(ctok)

  !  Set the default template namelist

  zval1=' '
  zval2=' '
  zval3=' '
  zval4=' '
  call get_environment_variable('LOCAL',zval1)
  call get_environment_variable('TRANSP_LOCATION',zval2)
  call get_environment_variable('CODESYSDIR',zval3)
  call get_environment_variable('UFILES_ROOT',zval4)

  if(zval3.ne.' ') then
     zprefix = '$CODESYSDIR/tables'
  else if(zval1.ne.' ') then
     zprefix = '$LOCAL/tables'
  else if(zval2.ne.' ') then
     zprefix = '$TRANSP_LOCATION'
  else
     zprefix = 'unknown'
  end if

  zgentmplfile=trim(zprefix)//'/rf_nml/toray.in'

  zufilespath=' '
  if(zval4.ne.' ') then
     zufilespath='$UFILES_ROOT/bin'
  endif
  !  Get the namelist input template filename

  ztmppath=trim(adjustl(toray_template))

  ! check if mdsplus

  kmdsplus = max(0,min(2,kmdsplus))

  if(kmdsplus.eq.2) then
     write(nonlin,*) &
          ' copy_toraynml: try direct read on MDS+ server...'
     status = Mds_Value('size(.NML:TORAY)', &
          idescr_long(inlines),isize)
     if (mod(status,2).ne.1) then
        call mdserr(nonlin, &
             'Mds_Value(''size(.NML:TORAY)'',..)', status)
        write(nonlin,*) ' %TORAY template not found in MDS+...'
        kmdsplus = 0
     else if(inlines.le.0) then
        write(nonlin,*) ' %TORAY template length is zero in MDS+'
        kmdsplus = 0
     else
        allocate(nltext(inlines)); nltext = ' '
        idims(1)=inlines
        idims(2)=0

        status = Mds_Value('.NML:TORAY', &
             idescr_cstringarr(NLTEXT,idims,1),isize)
        if (mod(status,2).ne.1) then
           call mdserr(nonlin, &
                'Mds_Value(''.NML:TORAY'',..)', status)
           write(nonlin,*) ' %TORAY template read error in MDS+...'
           kmdsplus = 0
        else
           ! OK...
           call find_io_unit(ilun)
           open(unit=ilun,file=ztarget,status='unknown')
           do iline=1,inlines
              write(ilun,'(A)') trim(nltext(iline))
           end do
           close(unit=ilun)
           deallocate(nltext)
           return  ! all done
        end if
        deallocate(nltext)
     end if
  end if

  if(kmdsplus.eq.1) then
     write(nonlin,*) ' copy_toraynml: template fetch via '// &
          '"mds_get_inf" subprocess...'
     zcmd=trim(adjustl(zufilespath))//'/mds_get_inf s '//trim(mds_server)// &
          ' t '//trim(mds_tree)// &
          ' r '//trim(runid)// &
          ' l .NML:TORAY q'
     write(nonlin,*) subr//trim(zcmd)
     call execute_command_line(trim(zcmd),exitstat=istat)
     ier=abs(istat)
     if (istat .ne. 0) then
        write(nonlin,*) subr, &
             ' failed to extract .NML:TORAY'
        write(nonlin,*) '   (now look for template file).'
        kmdsplus=0
     else
        zgentmplfile=trim(runid)//'..NML:TORAY'
        zcp='mv '
     end if
  end if

  if(kmdsplus.eq.0) then
     write(nonlin,*) subr,'namelist: TORAY_TEMPLATE= "', &
          trim(ztmppath),'"'
     if(ztmppath.NE.' ') then
        if((ztmppath(1:1).ne.'~').and.(ztmppath(1:1).ne.'/').and. &
             (ztmppath(1:1).ne.'.')) then
           write(nonlin,*) subr//' file tag: '//trim(ztmppath)
           ztag=ztmppath
           call uupper(ztag)
           zgentmplfile=trim(zprefix)//'/rf_nml/toray.in'
        else
           zgentmplfile=ztmppath
        end if
     else
        !     default template file
        write(nonlin,*) subr// &
             'TORAY_TEMPLATE is blank; Default used'
     end if
  end if

  write(nonlin,*) subr//'Toray template file= '
  write(nonlin,*) '    '//trim(zgentmplfile)

  zcmd=zcp//trim(adjustl(zgentmplfile))//' ' &
       //trim(adjustl(ztarget))
  write(nonlin,*) subr//'zcmd= ',trim(zcmd)

  call execute_command_line(trim(zcmd),exitstat=istat)
  ier=abs(istat)
  write(nonlin,*) subr//' file copy status (0 is normal): ',ier

  return
end subroutine copy_toraynml
