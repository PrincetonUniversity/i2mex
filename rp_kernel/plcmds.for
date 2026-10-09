      subroutine plcmds(zchar)
C
C  print out available calculator "commands"...
C
      use rpcalc_mod
C
      character*80 zbuf
C
      character*1 zchar
C
C-----------------------------------------------------------------
C
      lunt=lunzer(0)
C
      write(lunt,1001) zchar
 1001 format(/' RPLOT CALCULATOR commands:'/
     >     ' ** all start with the "',a1,'" character **'/)
C
      do ii=1,ncmdrpp
         i=icmdord(ii)
         zbuf=rppcmds(i)
         ilz=len_trim(zbuf)
         ilz=ilz+1
         zbuf(ilz:ilz)='('
         do j=1,nargrpp
            if(rppkeys(j,i).ne.' ') then
               ila=len_trim(rppkeys(j,i))
               zbuf(ilz+1:ilz+ila+1)=rppkeys(j,i)(1:ila)//','
               ilz=ilz+ila+1
            endif
         enddo
         zbuf(ilz:ilz)=')'
C
         write(lunt,'(/1x,a1,a)') zchar,zbuf(1:ilz)
C
         write(lunt,'('' argument defaults:'')')
         do j=1,nargrpp
            if(rppkeys(j,i).ne.' ') then
               if(rppadfs(j,i).eq.' ') then
                  zbuf='** no default, required argument'
                  ilz=len_trim(zbuf)
                  write(lunt,'(5x,a,'' =  '',a)')
     >               rppkeys(j,i),zbuf(1:ilz)
               else
                  zbuf=rppadfs(j,i)
                  ilz=len_trim(zbuf)
                  write(lunt,'(5x,a,'' = "'',a,''"'')')
     >               rppkeys(j,i),zbuf(1:ilz)
               endif
            endif
         enddo
C
         write(lunt,'('' summary description:'')')
         do k=1,3
            if(rppdescr(k,i).ne.' ') then
               zbuf=rppdescr(k,i)
               ilz=len_trim(zbuf)
               write(lunt,'(5x,a)') zbuf(1:ilz)
            endif
         enddo
C
      enddo
C
      write(lunt,'(//)')
C
      return
      end
