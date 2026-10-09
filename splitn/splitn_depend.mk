$(OBJDIR)/splitn_module.o : splitn_module.f90 
$(OBJDIR)/splitn_dcod.o : splitn_dcod.f90 
$(OBJDIR)/splitn_getsize.o : splitn_getsize.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_chkdp.o : splitn_chkdp.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/c_splitn.o : c_splitn.f90 
$(OBJDIR)/splitn_addwarn.o : splitn_addwarn.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_edit_enable.o : splitn_edit_enable.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_io_list.o : splitn_io_list.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_print_varinfo.o : splitn_print_varinfo.f90 
$(OBJDIR)/splitn_getd.o : splitn_getd.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_write.o : splitn_write.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_update_times.o : splitn_update_times.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_update_check.o : splitn_update_check.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_f77.o : splitn_f77.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_update_mods.o : splitn_update_mods.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_rhs.o : splitn_rhs.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_printwarn.o : splitn_printwarn.f90 
$(OBJDIR)/splitn_diff.o : splitn_diff.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_bdecode.o : splitn_bdecode.f90 
$(OBJDIR)/splitn_getf.o : splitn_getf.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_get.o : splitn_get.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_put.o : splitn_put.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_putw.o : splitn_putw.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_merge_write.o : splitn_merge_write.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_seek_update.o : splitn_seek_update.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_lunset.o : splitn_lunset.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_puta.o : splitn_puta.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_read.o : splitn_read.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_wtrans.o : splitn_wtrans.f90 $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_lhs.o : splitn_lhs.for $(OBJDIR)/splitn_module.o 
$(OBJDIR)/idcod_splitn.o : idcod_splitn.for 
$(OBJDIR)/splitn_nxgrp.o : splitn_nxgrp.for 
$(OBJDIR)/splitn.o : splitn.for $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_rhs0.o : splitn_rhs0.for $(OBJDIR)/splitn_module.o 
$(OBJDIR)/splitn_zeff_switches.o : splitn_zeff_switches.for 
