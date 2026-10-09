      subroutine plnames(zrunid)
C
      use cplotr_mod
C
      character*(*) zrunid      ! TRANSP runid
C
C----------------
C
      ilr=index(runid,' ')-1
      if(ilr.le.0) ilr=len(runid)
      tfiln=runid(1:ilr)//'TF.PLN'
      nfiln=runid(1:ilr)//'NF.PLN'
      mfiln=runid(1:ilr)//'MF.PLN'
C
      return
      end
