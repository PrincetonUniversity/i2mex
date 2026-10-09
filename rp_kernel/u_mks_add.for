      subroutine u_mks_add(zuns_orig,zuns_mks,zconv)
c
      use cplotr_mod
c
c  add an entry to the MKS units conversion table...
c
      character*(*) zuns_orig           ! original units (non-MKS perhaps)
      character*(*) zuns_mks            ! MKS units string
      real zconv                        ! conversion factor orig->mks
c
c  example:
c     call u_mks_add('n/cm3','/m^3',1.0e6)
c
c  if the target units string is not also in the original units table,
c  it is inserted with a conversion factor of 1.0
c
c--------------------------------
      character*32 ztmp  ! units -> uppercase
c
      integer lunzer
c
c--------------------------------
c
c  check if original units entry is known.  If so, issue a warning, this
c  is not expected.
c
      ztmp=zuns_orig
      call uupper(ztmp)
      ifnd=ifind_ordr(units_orig,iordru,nu_mks,ztmp)
      if(ifnd.eq.0) then
         nu_mks=nu_mks+1
         units_orig(nu_mks)=ztmp
         units_mks(nu_mks)=zuns_mks
         conv_units(nu_mks)=zconv
         call aordr_add(units_orig,iordru,nu_mks)
      else
         write(lunzer(0),*)
     >      ' %u_mks_add:  duplicate:  ',units_orig(ifnd)
         write(lunzer(0),*) ' passed: ',
     >      '   zuns_orig = ',zuns_orig,' zuns_mks = ',zuns_mks,
     >      ' conversion factor = ',zconv
         write(lunzer(0),*) ' stored: ',
     >      '   units_orig = ',units_orig(ifnd),' units_mks = ',
     >      units_mks(ifnd),' conversion factor = ',conv_units(ifnd)
      endif
c
c  check if mks units are known; if not, add them with conv. factor = 1.0
c
      ztmp=zuns_mks
      call uupper(ztmp)
      ifnd=ifind_ordr(units_orig,iordru,nu_mks,ztmp)
      if(ifnd.eq.0) then
         nu_mks=nu_mks+1
         units_orig(nu_mks)=ztmp
         units_mks(nu_mks)=zuns_mks
         conv_units(nu_mks)=1.0
         call aordr_add(units_orig,iordru,nu_mks)
      endif
c
      return
      end
