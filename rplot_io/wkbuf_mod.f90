module wkbuf_mod

  ! utility storage module for RPLOT_IO & friends

  implicit NONE
  SAVE
  PUBLIC

  real, dimension(:), allocatable :: zdatbuf   ! dynamically allocated buffer

CONTAINS
  !
  ! --- need_zdatbuf ---
  ! insure zdatbuf is big enough and then some
  !
  subroutine need_zdatbuf(ineed)
    integer, intent(in) :: ineed  ! size of buffer needed

    integer :: idatbuf   ! size of existing zdatbuf
    integer :: isize     ! size of new buffer
    integer :: istat     ! allocation status

    real, dimension(:), allocatable :: tmp_copy  ! transfer old contents to new buffer

    idatbuf = 0
    if (allocated(zdatbuf)) idatbuf=size(zdatbuf)

    if (ineed<=idatbuf) return  ! zdatbuf already big enough

    isize = max(20000,ineed + ineed/2)         ! just a little bit more

    !print *, '%wkbuf_mod: allocating buffer size = ',isize

    if(idatbuf.eq.0) then
       allocate(zdatbuf(isize),stat=istat)     ! create zdatbuf for first time
       if (istat/=0) goto 100

       zdatbuf = 0.0
    else
       allocate(tmp_copy(idatbuf),stat=istat)  ! allocate copy -- how can this fail?
       if (istat/=0) call errmsg_exit("?wkbuf_mod::need_zdatbuf: error allocating temp buffer")
       tmp_copy = zdatbuf

       deallocate(zdatbuf)
       allocate(zdatbuf(isize),stat=istat)
       if (istat/=0) goto 100

       zdatbuf(1:idatbuf) = tmp_copy
       deallocate(tmp_copy)

       zdatbuf(idatbuf+1:isize) = 0.0
    endif

    return

    ! -- crash --
100 print *, "?wkbuf_mod::need_zdatbuf: error allocating buffer with size = ",isize
    call bad_exit
  end subroutine need_zdatbuf
end module wkbuf_mod
