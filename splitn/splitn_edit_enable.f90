subroutine splitn_edit_enable(program_name)

  use splitn_module
  implicit NONE

  ! activate editing capability

  character*(*) program_name

  !-------------------------------------

  if(edit_started) then
     write(6,*) ' %splitn_edit_enable: edit already enabled & started by: ', &
          trim(edit_program)
     edit_started = edit_program.eq.program_name
  endif

  edit_program = program_name
  edit_enabled = .TRUE.

end subroutine splitn_edit_enable
