c-----------------------------------
      logical function ck_rplist(abbrev,substr)
c
      character*(*) abbrev              ! abbreviation to test (input)
      character*(*) substr              ! substring (input) (blank if none)
c
c  output:  function value:  TRUE if abbrev qualifies; FALSE otherwise
c
      character*21 zabbrev
c
c----------------
c
      if(abbrev(1:1).eq.'%') then
         ck_rplist=.FALSE.
         return
      endif
c
      if(substr.eq.' ') then
         ck_rplist=.TRUE.
         return
      endif
c
c  substring test required
c
      zabbrev=abbrev
      call trcaps(zabbrev)
      ils=len_trim(substr)
c
      if(index(zabbrev,substr(1:ils)).gt.0) then
         ck_rplist=.TRUE.
      else
         ck_rplist=.FALSE.
      endif
c
      return
      end
