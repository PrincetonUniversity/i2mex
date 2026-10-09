C----------------------------------------------------------------------
C  COPY_GENNML_SUB -- SUBROUTINE to copy Genray template namelist 
C
C  calling arguments only, trcom & TRDATA not used...
C
      SUBROUTINE COPY_GENNML_SUB(ctok,runid,genray_ech_template, 
     >        kmdsplus,mds_server,mds_tree,
     >        ztarget,
     >        nonlin,ier)

      use execsystem 
      IMPLICIT NONE

C-----------------------------------------------------------------------
C  subroutine arguments:
      
      character*(*), intent(in) :: ctok   ! tokamak ID
      character*(*), intent(in) :: runid  ! TRANSP run ID
      character*(*), intent(in) :: genray_ech_template  ! template name

      integer, intent(inout) :: kmdsplus  ! flag, =1 or 2 for MDS+ access
      !  if =2, assume MDS+ tree is open;
      !    RETURNED VALUE is set to 0 if MDSplus access failed.

      character*(*), intent(in) :: mds_server  ! MDS+ server
      character*(*), intent(in) :: mds_tree    ! MDS+ tree name

      character*(*), intent(in) :: ztarget   ! where to copy the file

      integer, intent(in) :: nonlin       ! I/O unit for messages

      integer, intent(out) :: ier ! status return code, 0=OK
      
      character*17 :: subr = 'COPY_GENNML_SUB: '
!
! MDSplus Fortran interface
!
      integer Mds_Value,idescr_long,idescr_cstringarr,idescr_cstring
      integer isize,status,idims(2)

!  buffer for text

      integer :: iline,inlines
      character*150, dimension(:), allocatable :: nltext

C-----------------------------------------------------------------------
C  local:

      INTEGER istat,iltok,ilnurd,ilun

      character*150 zval1,zval2,zval4
      character*200 ztmppath,ztag,zprefix,zufilespath
      character*300 zgentmplfile,zcmd
      character*3   zcp

C-----------------------------------------------------------------------
C
      ier=0
C
C----------

      zcp='cp'

      ILTOK=ILNURD(CTOK)

C  Set the default template namelist

      zval1=' '
      zval2=' '
      zval4=' '
      call get_environment_variable('LOCAL',zval1)
      call get_environment_variable('TRANSP_LOCATION',zval2)
      call get_environment_variable('UFILES_ROOT',zval4)

      if(zval1.ne.' ') then
         zprefix = '$LOCAL/tables'
      else if(zval2.ne.' ') then
         zprefix = '$TRANSP_LOCATION'
      else
         zprefix = 'unknown'
      endif

      zgentmplfile=trim(zprefix)//'/rf_nml/genray_ech.in'

      zufilespath=' '
      if(zval4.ne.' ') then
         zufilespath='$UFILES_ROOT/bin'
      endif
C  Get the namelist input template filename

      ztmppath=trim(adjustl(GENRAY_ECH_TEMPLATE))

! check if mdsplus 

      kmdsplus = max(0,min(2,kmdsplus))

      if(kmdsplus.eq.2) then
         write(nonlin,*) 
     >        ' copy_gennml_sub: try direct read on MDS+ server...'
         status = Mds_Value('size(.NML:GENRAY_ECH)',
     >        idescr_long(inlines),isize)
         if (mod(status,2).ne.1) then
            call mdserr(nonlin, 
     >           'Mds_Value(''size(.NML:GENRAY_ECH)'',..)', status)
            write(nonlin,*) ' %GENRAY template not found in MDS+...'
            kmdsplus = 0
         else if(inlines.le.0) then
            write(nonlin,*) ' %GENRAY template length is zero in MDS+'
            kmdsplus = 0
         else
            allocate(nltext(inlines)); nltext = ' '
            idims(1)=inlines
            idims(2)=0

            status = Mds_Value('.NML:GENRAY_ECH',
     >           idescr_cstringarr(NLTEXT,idims,1),isize)
            if (mod(status,2).ne.1) then
               call mdserr(nonlin, 
     >              'Mds_Value(''.NML:GENRAY_ECH'',..)', status)
               write(nonlin,*) ' %GENRAY template read error in MDS+...'
               kmdsplus = 0
            else
               ! OK...
               call find_io_unit(ilun)
               open(unit=ilun,file=ztarget,status='unknown')
               do iline=1,inlines
                  write(ilun,'(A)') trim(nltext(iline))
               enddo
               close(unit=ilun)
               deallocate(nltext)
               RETURN  ! all done
            endif
            deallocate(nltext)
         endif
      endif

      if(kmdsplus.eq.1) then
         write(nonlin,*) ' copy_gennml_sub: template fetch via '//
     >        '"mds_get_inf" subprocess...'
         zcmd=trim(adjustl(zufilespath))//'/mds_get_inf s '//
     >        trim(MDS_SERVER)//
     >        ' t '//trim(MDS_TREE)//
     >        ' r '//trim(runid)//
     >        ' l .NML:GENRAY_ECH q'
         write(nonlin,*) subr//trim(zcmd)
         istat=jsystem(trim(zcmd))
         if (istat .ne. 0) then
            write(nonlin,*) subr,
     >           ' failed to extract .NML:GENRAY_ECH'
            write(nonlin,*) '   (now look for template file).'
            kmdsplus=0
         else
            zgentmplfile=trim(runid)//'..NML:GENRAY_ECH'
            zcp='mv '
         endif
      endif

      if(kmdsplus.eq.0) then
         write(nonlin,*) subr,'namelist: GENRAY_ECH_TEMPLATE= "',
     >        trim(ztmppath),'"'
         IF(ztmppath.NE.' ') then
            IF((ztmppath(1:1).NE.'~').AND.(ztmppath(1:1).NE.'/').AND.
     >           (ztmppath(1:1).NE.'.')) then
               write(nonlin,*) subr//' file tag: '//trim(ztmppath)
               ztag=ztmppath
               call uupper(ztag)
               zgentmplfile=trim(zprefix)//'/rf_nml/'
     >              //'GENRAY_ECH_'//trim(adjustl(ztag))
     >              //'.IN'
            ELSE
               zgentmplfile=ztmppath
            ENDIF
         ELSE 
!     default template file
            write(nonlin,*) subr//
     >           'GENRAY_ECH_TEMPLATE is blank; Default used'
         ENDIF 
      endif

      write(nonlin,*) subr//'Genray template file= '
      write(nonlin,*) '    '//trim(zgentmplfile)

      zcmd=zcp//trim(adjustl(zgentmplfile))//' '
     >     //trim(adjustl(ztarget))
      write(nonlin,*) subr//'zcmd= ',trim(zcmd)

      istat=jsystem(trim(zcmd)) 
      ier=abs(istat)
      write(nonlin,*) subr//' file copy status (0 is normal): ',ier

      RETURN
      END
