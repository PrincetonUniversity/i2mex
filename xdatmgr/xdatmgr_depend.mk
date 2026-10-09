$(OBJDIR)/datmgr_mod.o : datmgr_mod.f90 
$(OBJDIR)/dmgini_chk.o : dmgini_chk.f90 $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/datmgr_no_delete.o : datmgr_no_delete.f90 $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmg_wkxpand.o : dmg_wkxpand.f90 $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/datmgr_delete_ok.o : datmgr_delete_ok.f90 $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/datadd.o : datadd.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmgalo.o : dmgalo.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmdloc.o : dmdloc.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dminew.o : dminew.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/wkput.o : wkput.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/wkget.o : wkget.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmprin.o : dmprin.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmidel.o : dmidel.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmgrpl.o : dmgrpl.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/datref.o : datref.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmg_macc_incr.o : dmg_macc_incr.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/dmgbsf.o : dmgbsf.for $(OBJDIR)/datmgr_mod.o 
$(OBJDIR)/iget_naxxvr.o : iget_naxxvr.for 
$(OBJDIR)/dmista.o : dmista.for $(OBJDIR)/datmgr_mod.o 
