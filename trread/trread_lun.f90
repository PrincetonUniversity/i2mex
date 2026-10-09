subroutine trread_lun(ilun)

  ! return LUN_TF in the subroutine argument -- available LUN for tmp file i/o

  use cplotr_mod

  ilun = lun_tf

end subroutine trread_lun
