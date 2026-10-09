subroutine trx_nb_sublist(ss_tr,ss_ref,inum_tr,inb_sublist,ilun,imatch,istat)

  ! compare reference neutral beam machine description to TRANSP description
  ! the TRANSP description is expected (hoped?) to be a subset with matching
  ! order

  use plasma_state_mod
  implicit NONE

  type (plasma_state) :: ss_tr   ! TRANSP state & machine description
  type (plasma_state) :: ss_ref  ! reference machine description

  integer :: inum_tr             ! number of beams, TRANSP machine description
  integer, intent(out) :: inb_sublist(inum_tr)  ! mapping tr -> ref

  integer, intent(in) :: ilun    ! LUN for messages
  integer, intent(out) :: imatch ! =1 if lists are identical, =0 otherwise
  integer, intent(out) :: istat  ! =0 if tr list matches or is ordered subset
                                 ! =1 otherwise

  !-----------------------------------
  ! local:

  integer :: inum_ref,ib_tr,ib_ref
  real*8, parameter :: ztol = 0.5d-3  ! 1/2 mm
  real*8 :: zdiff

  !-----------------------------------

  inb_sublist(:) = 0

  imatch=0
  istat=1

  inum_ref = ss_ref%nbeam

  if(inum_ref.lt.inum_tr) then
     write(ilun,*) ' %trx_nb_sublist: TRANSP machine description #beams = ',inum_tr
     write(ilun,*) '  reference machine description #beams = ',inum_ref
     write(ilun,*) '  TRANSP count exceeds reference count.'
     return
  endif

  !  OK total beam counts are consistent
  istat=0

  ib_ref = 0
  do ib_tr = 1,inum_tr
     do
        ib_ref = ib_ref + 1
        if(ib_ref.gt.inum_ref) then
           istat=1
           inb_sublist(:)=0
           write(ilun,*) ' %trx_nb_sublist: could not match TRANSP beam# ',ib_tr
           write(ilun,*) '  with beam from reference list.'
           exit
        endif

        zdiff = max(abs(ss_tr%sRtcen(ib_tr)-ss_ref%sRtcen(ib_ref)), &
             abs(ss_tr%Lbsctan(ib_tr)-ss_ref%Lbsctan(ib_ref)), &
             abs(ss_tr%Lbscap(ib_tr)-ss_ref%Lbscap(ib_ref)), &
             abs(ss_tr%Zbsc(ib_tr)-ss_ref%Zbsc(ib_ref)), &
             abs(ss_tr%Zbap(ib_tr)-ss_ref%Zbap(ib_ref)) )

        if(zdiff.gt.ztol) cycle

        ! match OK

        inb_sublist(ib_tr) = ib_ref
        exit
     enddo
     if(istat.eq.1) exit
  enddo

  if(istat.eq.0) then
     if(inum_tr.eq.inum_ref) imatch = 1
     if(imatch.eq.1) then
        write(ilun,*) ' %trx_nb_sublist:  TRANSP and reference beam lists match fully.'
     else
        write(ilun,*) ' %trx_nb_sublist:  TRANSP beams sublist is ordered subset of reference list.'
     endif
  endif

end subroutine trx_nb_sublist
