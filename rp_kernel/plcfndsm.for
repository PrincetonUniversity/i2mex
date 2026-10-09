C                                                                PLCFNDSM.FOR
 
c-------------------------------------------------------------
      subroutine PLCfndsm(symbol, symtab, symidx)          ! TBT 8/90
c
 
c This subroutine searches for the occurrence of an exact match symbol in
c a table of recognized names passed in the character array symtab.  the
c parameter symidx is returned with index of the matching name or is set
c to -1 if no match is found.
c
c  mod dmc -- using integer function str_cbc (VAXONLY library) instead
c  of VAX RTL routine STR$CASE_BLIND_COMPARE.  This improves portability
c  to unix workstations
c
      character*(*) symbol, symtab(*)
      integer symidx
      integer, external :: str_cbc
 
c
c  test for an exact match
      symidx=-1
      i=1
      do while (symtab(i)(1:1) .ne. ';')
       if (str_cbc(symbol,symtab(i)).eq.0) then
          symidx=i
          return
       end if
       i=i+1
      end do
      return
      end
