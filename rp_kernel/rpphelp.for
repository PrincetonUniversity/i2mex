      subroutine rpphelp(ilunt,zpre)
c
c  dmc 10 Aug 1999 -- write out help to describe auxilliary commands
c   of rplot calculator
c
      use rpcalc_mod

      integer ilunt                     ! i/o unit to write on
      character*(*) zpre                ! prefix character(s) for commands
c
c-----------------------------
c
      character*75 zbuf
      character*1 ztrl
c-----------------------------
c
      zbuf=' '
      ilpre=max(1,len_trim(zpre))
      ilzb=len(zbuf)
c
      write(ilunt,999)
 999  format(/
     >   ' RPLOT Calculator Auxilliary Commands:'/
     >   ' ====================================='//
     >   ' Auxilliary commands are used for functions that cannot'/
     >   ' readily be incorporated into the calculator expression'/
     >   ' evaluation language, either because the command performs'/
     >   ' an i/o function rather than a computation, or, because the'/
     >   ' computation involved does not fit the implicit data type'/
     >   ' assumptions of the main evaluation language.')
 
      write(ilunt,1000)
 1000 format(/
     >' There are two general syntax options for auxilliary commands:'/
     >' Positional argument syntax:'/
     >'   <command>(<arg1>,<arg2>,...,<argN>)'/
     >'   where <arg#> is either the desired value or empty if the ',
     >    'default'/
     >'   value is to be used.  Trailing ",,)" may be abbreviated as ',
     >    '")".'/
     >' Keyword argument syntax:'/
     >'   <command>(<key1>=<arg1>,...,<keyN>=<argN>)'/
     >'   where <key#> identifies the argument and <arg#> gives the',
     >    ' value.'/
     >'   Keyword/argument pairs can be given in any order.  Missing'/
     >'   keyword arguments get their default value.'/
     >' Most commands have an argument with keyword "EXPR" which ',
     >    'identifies'/
     >' data to which the command applies.  Sometimes there are ',
     >    'constraints'/
     >' as to the type of data.  The default value, EXPR="$", refers'/
     >' to the data computed by the most recent preceding calculator'/
     >' command or expression evaluation.')
 
      write(ilunt,1001) zpre(1:ilpre),zpre(1:ilpre),zpre(1:ilpre),
     >   zpre(1:ilpre),zpre(1:ilpre),zpre(1:ilpre)
 1001 format(/' Examples of Auxilliary Commands:'//
     >'   ',a,'SMOOTH(0.1,0.05,,,NE*MAX(TE,TI)) -or-'/
     >'   ',a,'SMOOTH(DELTA_T=0.1,DELTA_X=0.05,EXPR=NE*MAX(TE,TI))'/
     >' --smooth "NE*MAX(TE,TI)" using weighting functions with'/
     >' delta(t) of +/- 0.1 seconds, delta(x)= +/- 0.05, here "x"'/
     >' is the normalized flux coordinate of the NE, TE, TI profiles.'/
     >' Parameters EPS_X and EPS_T are defaulted.'//
     >'   ',a,'SAVE(PE,Electron Pressure,Pascals) -or-'/
     >'   ',a,'SAVE(ABBREV=PE,LABEL=Electron Pressure,UNITS=Pascals)'/
     >' --EXPR is defaulted, so, save the current calculator ',
     >   'accumulator'/
     >' contents as a named function "PE" with the given label and '/
     >' units label.'//
     >'   ',a,'DELETE(*) -or- ',a,'DELETE(FCN=*)'/
     >' --delete all user defined functions and reclaim the associated'/
     >' memory.')
 
      write(ilunt,1002)
 1002 format(/
     >' A note on command parsing.  Inside the parentheses, the two'/
     >' characters "=" (equals sign) and "," (comma) are significant.'/
     >' If a label needs to contain these characters than it must be'/
     >' enclosed in quotes.  A comma in an expression (EXPR argument)'/
     >' is safe without quotes because the comma will be appearing at'/
     >' a lower level of nesting inside parentheses.')
C
      write(ilunt,1003)
 1003 format(/' Complete List of RPLOT Calculator Auxilliary Commands:')
c
      do i=1,ncmdrpp
         write(ilunt,'(1x)')
         ilc=len_trim(rppcmds(i))
         zbuf=zpre(1:ilpre)//rppcmds(i)(1:ilc)//'('
         il=len_trim(zbuf)
         ztrl=','
         do j=1,ncmdargs(i)
            if(j.eq.ncmdargs(i)) ztrl=')'
            ilk=len_trim(rppkeys(j,i))
            if((il+ilk+1).gt.ilzb) then
               write(ilunt,'(1x,A)') zbuf(1:il)
               il=6
               zbuf=' '
            endif
            zbuf(il+1:il+ilk+1)=rppkeys(j,i)(1:ilk)//ztrl
            il=il+ilk+1
         enddo
         write(ilunt,'(1x,A)') zbuf(1:il)
         write(ilunt,1005) exprtype(i)
 1005    format(' ...EXPR argument type:  ',a/' ...Description:')
         do j=1,3
            if(rppdescr(j,i).ne.' ') then
               ild=len_trim(rppdescr(j,i))
               write(ilunt,'(6x,a)') rppdescr(j,i)(1:ild)
            endif
         enddo
         write(ilunt,1007)
 1007    format(' ...argument defaults:')
         do j=1,ncmdargs(i)
            ilk=len_trim(rppkeys(j,i))
            ild=len_trim(rppadfs(j,i))
            if(rppadfs(j,i).ne.' ') then
               write(ilunt,1010) rppkeys(j,i)(1:ilk),rppadfs(j,i)(1:ild)
            else
               write(ilunt,1011) rppkeys(j,i)(1:ilk)
            endif
 1010       format(6x,a,t18,'= ',a)
 1011       format(6x,'no default for ',a,', value must be specified.')
         enddo
         if(.not.lrepacc(i)) write(6,1013)
 1013    format(
     >     ' This command has no effect on the calculator accumulator.')
      enddo
c
      return
      end
