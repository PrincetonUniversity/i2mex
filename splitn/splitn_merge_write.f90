subroutine splitn_merge_write(merge_input_filename, &
     backup_filename,output_filename, tagname, ios)

  ! save backup of current namelist, if backup_filename.ne.' '

  ! merge namelist text with text from a namelist section
  !   (in file merge_input_filename)

  ! write merged list to output_filename
  ! read merged file back in if no error so far.
  
  ! tag merged lines with "!$<tagname>".

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: merge_input_filename
  character*(*), intent(in) :: backup_filename
  character*(*), intent(in) :: output_filename
  character*(*), intent(in) :: tagname

  integer, intent(out) :: ios   ! status code on exit, 0=normal

  !----------------------------------
  !----------------------------------

  ios=0
  if(backup_filename.ne.' ') then

     call write_nl(backup_filename,ios)
     if(ios.ne.0) return

  endif

  call merge_read_nl(merge_input_filename,tagname,ios)
  if(ios.ne.0) return

  call write_nl(output_filename,ios)

  if(ios.ne.0) then

     call init_nl(.TRUE.)  ! leaving empty namelist -- error occurred.

  else

     call read_nl(output_filename,ios)

  endif

end subroutine splitn_merge_write
