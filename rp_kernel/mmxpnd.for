      subroutine mmxpnd(zxpres,zwk)
C
C  expand MIN(a,b,c) -> MIN(a,MIN(b,c))
C  expand MAX(a,b,c) -> MAX(a,MAX(b,c))
C
      character*(*) zxpres              ! input expression, in uppercase
      character*(*) zwk                 ! workspace
C
C-------------------------------
C
      character*38 zlgl
      character*1 zchar,zquot,zquot0
C
C
      data zlgl/'1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZ_$'/
C
C------------------------------
C
      lunt=lunzer(0)
C
 10   continue
C
C  look for MIN or MAX operator tokens
C
      ilen=len_trim(zxpres)
      zquot0=' '
      do ic=1,ilen
         ifnd=0
         if(zquot0.ne.' ') then
            if(zxpres(ic:ic).eq.zquot0) zquot0=' '
         else if((zxpres(ic:ic).eq.'"').or.(zxpres(ic:ic).eq.'''')) then
            zquot0=zxpres(ic:ic)
         else if(zxpres(ic:ic).eq.'M') then
            iok=0
            if(ic.eq.1) then
               iok=1
            else if(ic.gt.ilen-4) then
               iok=0
            else
               zchar=zxpres(ic-1:ic-1)
               iok=1
               if(index(zlgl,zchar).gt.0) iok=0
            endif
            if(iok.eq.1) then
C  have an "M" which is the start of a token
               if((zxpres(ic:ic+2).ne.'MIN').and.
     >            (zxpres(ic:ic+2).ne.'MAX')) then
                  iok=0
               endif
            endif
            if(iok.eq.1) then
C  have a token which starts with MIN or MAX; look for following paren
               do ic2=ic+3,ilen
                  if((zxpres(ic2:ic2).ne.' ').and.
     >               (zxpres(ic2:ic2).ne.char(9))) then
                     if(zxpres(ic2:ic2).ne.'(') iok=0
                     go to 20
                  endif
               enddo
               iok=0
 20            continue
            endif
            if(iok.eq.1) then
C  have a MIN or MAX with following "(" at position ic2
C  look for occurrence of more than two arguments
C  output of this loop:  iok=1 if there are more than two arguments, and
C    iarg2 = position just before start of 2nd argument;
C    iarg3 = position of end of last argument
C
C  method:  count commas at current parentheses level; find closing ")"
               iparen=1                 ! current level =1
               icommas=0
               zquot=' '
               do ic3=ic2+1,ilen
                  if(zquot.ne.' ') then
                     if(zxpres(ic3:ic3).eq.zquot) then
                        zquot=' '       ! closing quote
                     endif
                  else
                     if(zxpres(ic3:ic3).eq.'(') then
                        iparen=iparen+1
                     else if(zxpres(ic3:ic3).eq.')') then
                        iparen=iparen-1
                        iarg3=ic3-1
                        if(iparen.eq.0) go to 50 ! closing paren
                     else if((zxpres(ic3:ic3).eq.'''').or.
     >                       (zxpres(ic3:ic3).eq.'"')) then
                        zquot=zxpres(ic3:ic3) ! opening quote
                     else
                        if(zxpres(ic3:ic3).eq.',') then
C  count comma if not in quote string and at correct paren level
                           if((zquot.eq.' ').and.(iparen.eq.1)) then
                              icommas=icommas+1
                              if(icommas.eq.1) iarg2=ic3
                           endif
                        endif
                     endif
                  endif
               enddo
               iok=0
 50            continue
               if(icommas.gt.1) then
                  iok=1
               else
                  iok=0
               endif
            endif
            if(iok.eq.1) then
C
C  expand min/max and re-enter loop at very top
C
C                                     |here is replecated MIN( or MAX(
               zwk = zxpres(1:iarg2)//zxpres(ic:ic+2)//'('//
     >            zxpres(iarg2+1:iarg3)//')'//zxpres(iarg3+1:ilen)
C                                         ^here is new closing paren
               zxpres=zwk
               go to 10
            endif
         endif
      enddo
C
      return
      end
