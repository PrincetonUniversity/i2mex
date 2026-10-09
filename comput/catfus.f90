subroutine catfus_num(ireact)
  use catfus_mod
  implicit none
 
  integer,intent(out) :: ireact
 
  ireact = nreact            ! number of known reaction categories

  return
end subroutine catfus_num


subroutine catfus_btindx(c3beam,c3target,itypbeam,indx)
  use catfus_mod
  implicit none

  character(len=*), intent(in) :: c3beam   ! "D" "T" or "HE3"
  character(len=*), intent(in) :: c3target ! "D:P" "D:N" "T" or "HE3"
  integer, intent(in) :: itypbeam     ! 1 for beam, 2 for fusion prod, 3 for
                                      ! RF tail; 0 for thermal

  integer, intent(out) :: indx        ! reaction index (or 0 if none found)

  !----------------------------------------------
  character(len=3) :: z3beam,z3target
  integer :: ityptarg = 0
  integer i
  !----------------------------------------------

  indx=0
  if(itypbeam.lt.0) return
  if(itypbeam.gt.3) return

  z3beam=c3beam;     call uupper(z3beam)

  if((z3beam.ne.'D  ').and.(z3beam.ne.'T  ').and. &
       (z3beam.ne.'HE3')) then
    return
  endif

  z3target=c3target; call uupper(z3target)
  if((z3target.ne.'D:P').and.(z3target.ne.'D:N').and.(z3target.ne.'D  ').and. &
       (z3target.ne.'T  ').and.(z3target.ne.'HE3')) then
    return
  endif

  do i=1,nreact
    if((r1(i).eq.z3beam).and.(r2(i).eq.z3target)) then
      if((typmin(i).le.itypbeam).and.(typmax(i).ge.itypbeam)) then
        if(ityptarg.le.typ2max(i)) then
          indx = i
          exit
        endif
      endif
    endif
  enddo

  return
end subroutine catfus_btindx
 
 
subroutine catfus_name(ireact,irtyp1,irtyp2,name)

  !  return a sanitized version of reaction name: first ":" => "_"; any other
  !  punctuation removed.

  !  for original name: use catfus_nam0, same arguments.

  !  mod DMC: if irtyp2=-1, return a "generic" name; pass 0 to catfus_nam0

  integer,intent(in) :: ireact    ! reaction index
  integer, intent(in) :: irtyp1   ! type of 1st reagent
  integer, intent(in) :: irtyp2   ! type of 2nd reagent

  character(len=*), intent(out) :: name   ! name (label string) for reaction

  !------------------------------------
  integer :: j,k,ilen,icolon,irtyp2a
  !------------------------------------

  irtyp2a=irtyp2
  if(irtyp2a.eq.-1) irtyp2a=0

  name=' '
  call catfus_nam0(ireact,irtyp1,irtyp2a,name)

  ilen=len(trim(name))

  icolon=0
  k=0

  do j=1,ilen
    if(name(j:j).eq.':') then
      icolon=icolon+1
      if(icolon.eq.1) then 
        k=k+1
        name(k:k)='_'
      endif

    else if(name(j:j).eq.'-') then
      continue  ! i.e. skip

    else if(name(j:j).eq.'>') then
      continue  ! i.e. skip

    else
      k=k+1
      if(k.lt.j) name(k:k)=name(j:j)
    endif
  enddo

  if(k.lt.ilen) name(k+1:ilen)=' '

  if(irtyp2a.ne.irtyp2) then
    if(name(1:3).eq.'BT_') then
      if(ireact.le.6) then
        name = name(4:)
      endif
    else if(name(1:3).eq.'FT_') then
      name(2:2)='P'
    else if(name(1:3).eq.'RT_') then
      name(2:2)='F'
    endif
  endif

  return
end subroutine catfus_name


subroutine catfus_nam0(ireact,irtyp1,irtyp2,name)

  use catfus_mod
  implicit none

  integer, intent(in) :: ireact   ! reaction index
  integer, intent(in) :: irtyp1   ! type of 1st reagent
  integer, intent(in) :: irtyp2   ! type of 2nd reagent

  !  rule:  irtyp1.ge.irtyp2 .and. irtyp2.le.1 .and. irtyp2.ge.0
  !  additional rules:  see the module data
  !     generally, reaction indices .gt.8 are for FP and RFI ions only.
  !
  !  irtyp1|2 = 0 -- thermal   irtyp1=0, irtyp2=0: thermonuclear reaction
  !
  !     in what follows, "target" means thermal target population.
  !
  !  irtyp1|2 = 1 -- beam ion  irtyp1=1, irtyp2=0: beam-target reaction
  !                            irtyp1=1, irtyp2=1: beam-beam reaction
  !
  !  irtyp1|2 = 2 -- fusion product (FP):  T or He3
  !                            irtyp1=2, irtyp2=0: FP-target reaction
  !                            irtyp1=2, irtyp2=1: FP-beam reaction
  !
  !    (D is not a fusion product; P and He4 are non-reacting fusion products)
  !
  !  irtyp1|2 = 3 -- RF tail ion (RFI)
  !                            irtyp1=3, irtyp2=0: RFI-target reaction
  !                            irtyp1=3, irtyp2=1: RFI-beam reaction
  !    (generally, T or He3 tail ions are of interest).
  !  D tail has been added , MG 08/19/19
  character(len=*), intent(out) :: name   ! name (label string) for reaction

  !----------------------------------------
  character(len=3) :: reagent1,reagent2
  character(len=8) :: combo
  integer ilr1,ilr2
  !----------------------------------------

  name='error'

  if(irtyp1.lt.irtyp2) return
  if(irtyp2.gt.1) return
  if(irtyp2.lt.0) return

  if((ireact.le.0).or.(ireact.gt.nreact)) return
  if((irtyp1.lt.typmin(ireact)).or.(irtyp1.gt.typmax(ireact))) return
  if(irtyp2.gt.typ2max(ireact)) return

  call catfus_reagent(ireact,reagent1,reagent2)
  ilr1=len_trim(reagent1)
  ilr2=len_trim(reagent2)

  if(irtyp1.eq.irtyp2) then

    combo=reagent1(1:ilr1)//'-'//reagent2(1:ilr2)

    if(irtyp1.eq.0) then
      name='TH:'//combo
    else if(irtyp1.eq.1) then
      name='BB:'//combo
    endif

  else

    combo=reagent1(1:ilr1)//'->'//reagent2(1:ilr2)

    if(irtyp1.eq.1) then
      name='BT:'//combo
    else if(irtyp1.eq.2) then
      if(irtyp2.eq.0) then
        name='FT:'//combo
      else if(irtyp2.eq.1) then
        name='FB:'//combo
      endif
    else if(irtyp1.eq.3) then
      if(irtyp2.eq.0) then
        name='RT:'//combo
      else if(irtyp2.eq.1) then
        name='RB:'//combo
      endif
    endif

  endif

  return
end subroutine catfus_nam0


subroutine catfus_mtype(ireact,itypfast)
 
  use catfus_mod
  implicit none
 
  integer,intent(in) :: ireact    ! reaction index
  integer,intent(out) :: itypfast ! fast-ion type
 
  !  generally, the categorization allows fast-ion -- fast-ion reactions
  !  to be considered, but the "second" fast-ion specie must be a beam
  !  ion specie; the "first" specie can be beam-ion (B),
  !  fusion-product-ion (FP), or RF-tail-ion (RFI).
  !
  !  categories for FP-FP, RFI-RFI, and FP-RFI are not supported;
  !  categories B-B, FP-B, and RFI-B *are* supported.
  !
  !----------------------------------------

  itypfast=0
  if((ireact.le.0).or.(ireact.gt.nreact)) return

  itypfast=typmax(ireact)   ! 1=beam, 2=fusion product, 3=RF tail

  return
end subroutine catfus_mtype


subroutine catfus_reagent(ireact,reagent1,reagent2)
  use catfus_mod
  implicit none

  integer,intent(in) :: ireact    ! reaction index

  character(len=*),intent(out) :: reagent1,reagent2    ! reagents

  !----------------------------------------

  reagent1='?'
  reagent2='?'

  if((ireact.le.0).or.(ireact.gt.nreact)) return

  reagent1 = r1(ireact)
  reagent2 = r2(ireact)

  return
end subroutine catfus_reagent
