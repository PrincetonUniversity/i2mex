      subroutine set_rpluns(ilun_tf,ilun_nf,ilun_mf)
c
      use cplotr_mod
c
      if(ilun_tf.ne.0) lun_tf=ilun_tf
      if(ilun_nf.ne.0) lun_nf=ilun_nf
      if(ilun_mf.ne.0) lun_mf=ilun_mf
c
      return
      end
c-------------------
      subroutine get_rpluns(ilun_tf,ilun_nf,ilun_mf)
c
      use cplotr_mod
c
      ilun_tf=lun_tf
      ilun_nf=lun_nf
      ilun_mf=lun_mf
c
      return
      end
