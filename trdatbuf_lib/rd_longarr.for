      subroutine rd_longarr(arg,array,numa,ier)
C
C  read data at expression (arg) into array (array)
C  return error code ier, 0 = normal
C
C---------------------------------
C History:
C 11/18/99 CAL: new calling arguments for mds_value
C
C------------------------------
      use UFILES
C
      character*(*) arg
      integer array(numa)
      integer ier
C
      character*128 zerrmsg
C
      integer     idescr_longarr, iretl
C
C----------------------------------
C
      ier=0
      ila=max(1,len_trim(arg))
C
      istat=mds_value(arg(1:ila),idescr_longarr(array,numa,1),iretl)
      if(mod(istat,2).ne.1) then
         call ufmds_err('UFMDS_READ',istat,ier,zerrmsg,ilm)
         write(ufelun,9510) arg(1:ila),zerrmsg(1:ilm)
 9510    format('  reading MDS+ signal data: ',a/2x,a)
      endif
C
      return
      end
