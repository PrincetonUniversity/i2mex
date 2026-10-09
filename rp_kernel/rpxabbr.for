      subroutine rpxabbr(abr_in,abr_out)
c
      implicit NONE
c
      character*(*), intent(in) :: abr_in ! <abbrev> or <abbrev>$<runid>
      character*(*), intent(out) :: abr_out ! <abbrev> or <runid>
c
c  if the label abr_in contains no "$" character, just copy it;
c  if it does contain "$" and it is the last character, copy replacing
c  the "$" with an "&"; if it contains "$" and it is not the last character,
c  copy the last part (after the "$") only.
c
c  abr_in is expected to be longer than abr_out
c
      integer ilen,idollr,ilo
c
      ilen=len_trim(abr_in)
      idollr=index(abr_in,'$')
      ilo=len(abr_out)
c
      if(idollr.eq.0) then
         abr_out = abr_in(1:ilo)
      else if(idollr.eq.ilen) then
         abr_out = abr_in(1:ilo)
         if(idollr.le.ilo) abr_out(idollr:idollr)='&'
      else
         abr_out = abr_in(idollr+1:min(ilen,idollr+ilo))
      endif
c
      return
      end
