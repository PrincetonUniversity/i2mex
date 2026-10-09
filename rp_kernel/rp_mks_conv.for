      subroutine rp_mks_conv(zitem,zuns_in,zuns_mks,zconv,iwarn)
c
      use cplotr_mod
c
c  find MKS units conversion factor for item "zitem" which prior
c  to conversion has units "zuns_in".
c
c  input:
      character*(*) zitem               ! item name -- for message if needed
      character*(*) zuns_in             ! item units (to be converted to MKS)
c
c  output:
      character*(*) zuns_mks            ! MKS units
      real zconv                        ! conversion factor, old units -> MKS
c
      integer iwarn                     ! iwarn=0 if conversion successful
c
c----------------
c
c  if the conversion fails, then a message is written on lunzer(0),
c  and, zuns_mks=zuns_in, zconv=1.0, iwarn=1
c
c  if the conversion succeeds zuns_mks contains the MKS units string,
c  and zconv the conversion factor, and, iwarn=0
c
c----------------
c
      character*32 ztmp  ! tmporary units string ^uppercase
c
c----------------
c
      ztmp=zuns_in
      call uupper(ztmp)
c
      i=ifind_ordr(units_orig,iordru,nu_mks,ztmp)
      if(i.eq.0) then
c
c  no match
c
         write(lunzer(0),*)
     >      ' %rp_mks_conv:  no MKS conversion available.'
         write(lunzer(0),*)
     >      '  item = "',zitem,'"  units = "',zuns_in,'"'
c
         zuns_mks=zuns_in
         zconv=1.0
         iwarn=1
c
      else
c
         zuns_mks=units_mks(i)
         zconv=conv_units(i)
         iwarn=0
c
      endif
c
      return
      end
