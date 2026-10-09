logical function tdb_logchk_special(d,zitem,iwarn)
  ! dmc Apr 2005
  ! logical query -- presence (T) or absence (F) of "special" data items
  ! in trdat buffer.

  ! the items covered are those not handled by the code generator:
  !   "RPL" -- TF ripple vs. (R,Z)
  !   "RP2" -- TF ripple vs. (R,Z) 2nd component
  !   "NB2" -- beam power/voltage/etc data vs. time
  !   "RFP" -- ICRF antenna power vs. time
  !   "RFF" -- ICRF frequencies vs. time
  !   "LHP" -- LH antenna powers vs. time
  !   "ECP" -- ECH antenna powers vs. time
  !   "ECA" -- ECH poloidal aiming vs. time
  !   "ECB" -- ECH toroidal aiming vs. time
  !   "MMX" -- complete MHD equilibrium (inside plasma bdy) vs. time -- moments
  !   "RFS" -- complete MHD equilibrium (inside plasma bdy) vs. time -- spline
  !   "PSI" -- free boundary MHD equilibrium data Psi(R,Z,t)
  !   "LIM" -- EFIT style (R,Z) contour piecewise linear axisymmetric limiter
  !   "SAW" -- sawtooth event times
  !   "PEL" -- pellet event times
  !
  ! cf codesys/source/misc/trdatgen.spec -- all the above "trigraphs" are 
  ! associated with "special_handling" i.e. hand coded data channels.  Many,
  ! but not necessarily all, special handling channels are supported in this
  ! routine.

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: zitem   ! item queried
  integer, intent(out) :: iwarn        ! 0: OK; 1: item not recognized

  ! in case an item is unrecognized, this function prints out a warning,
  ! returns iwarn=1, and function value .FALSE.
  !----------------------------------------------
  integer lunmsg_tdb
  character*20 zitemc
  !----------------------------------------------

  iwarn=0
  tdb_logchk_special = .FALSE.

  zitemc=zitem
!-call uupper(trim(zitemc))   jec 22Dec2009  fails under debug compilation
  call uupper(     zitemc )

  if(zitemc.eq.'RPL') then
     tdb_logchk_special = ( d%lrpl .gt. 0 )

  else if(zitemc.eq.'RP2') then
     tdb_logchk_special = ( d%lrp2 .gt. 0 )

  else if(zitemc.eq.'NB2') then
     tdb_logchk_special = ( d%lpwrnb .gt. 0 )

  else if(zitemc.eq.'RFP') then
     tdb_logchk_special = ( d%lpwrrf .gt. 0 )

  else if(zitemc.eq.'RFF') then
     tdb_logchk_special = ( d%lfrqrff .gt. 0 )

  else if(zitemc.eq.'LHP') then
     tdb_logchk_special = ( d%lpwrlh .gt. 0 )

  else if(zitemc.eq.'ECP') then
     tdb_logchk_special = ( d%lpwrec .gt. 0 )

  else if(zitemc.eq.'ECA') then
     tdb_logchk_special = ( d%lfeca .gt. 0 )

  else if(zitemc.eq.'ECB') then
     tdb_logchk_special = ( d%lfecb .gt. 0 )

  else if(zitemc.eq.'MMX') then
     tdb_logchk_special = ( d%ldmmx .gt. 0 ) .and. ( d%nxmmx .gt. 1 )

  else if(zitemc.eq.'RFS') then
     tdb_logchk_special = ( d%lrfs .gt. 0 )

  else if(zitemc.eq.'PSI') then
     tdb_logchk_special = ( d%lfpsi .gt. 0 )

  else if(zitemc.eq.'LIM') then
     tdb_logchk_special = ( d%llim .gt. 0 )

  else if(zitemc.eq.'SAW') then
     tdb_logchk_special = ( d%ltsaw .gt. 0 )

  else if(zitemc.eq.'PEL') then
     tdb_logchk_special = ( d%npelda .gt. 0 )
     
  else
     write(lunmsg_tdb(0),*) &
          ' ?? trdatbuf_lib/tdb_logchk_special -- unrecognized data type: ', &
          trim(zitemc)
     iwarn=1
  endif

end function tdb_logchk_special

logical function tdb_logchk_nbi(d,zsubset,iwarn)
  ! dmc Apr 2005
  ! logical query -- presence (T) or absence (F) of NBI-related data items
  ! in trdat buffer.
  !
  ! if NB2 data is not in the buffer at all, always return F = FALSE.
  !
  ! known subset names:
  !   "PWR" -- beam power
  !   "VLT" -- beam voltage (full energy injected ptcls)
  !   "FUL" -- full energy fraction (#/sec)/(total #/sec)
  !   "HLF" -- half energy fraction (#/sec)/(total #/sec)
  !
  ! this code is matched to code in trdatusub which acquires NBI-related
  ! time dependent data...

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: zsubset ! subset item queried
  integer, intent(out) :: iwarn        ! 0: OK; 1: item not recognized

  ! in case an item is unrecognized, this function prints out a warning,
  ! returns iwarn=1, and function value .FALSE.
  !----------------------------------------------
  integer lunmsg_tdb
  character*20 zitemc
  !----------------------------------------------

  zitemc = zsubset
!-call uupper(trim(zitemc))   jec 22Dec2009  fails under debug compilation
  call uupper(     zitemc )

  iwarn = 0
  tdb_logchk_nbi=.FALSE.

  ! if there is no power data, there is no NBI data at all...

  if(d%lpwrnb.eq.0) then
     return
  endif

  ! there is some data...

  if(zitemc.eq.'PWR') then
     tdb_logchk_nbi = ( d%lpwrnb .gt. 0 )

  else if(zitemc.eq.'VLT') then
     tdb_logchk_nbi = ( d%lvltnb .gt. 0 )

  else if(zitemc.eq.'FUL') then
     tdb_logchk_nbi = ( d%lfulnb .gt. 0 )

  else if(zitemc.eq.'HLF') then
     tdb_logchk_nbi = ( d%lhlfnb .gt. 0 )

  else
     write(lunmsg_tdb(0),*) &
          ' ?? trdatbuf_lib/tdb_logchk_nbi -- unrecognized subset type: ', &
          trim(zitemc)
     iwarn=1
  endif

end function tdb_logchk_nbi
