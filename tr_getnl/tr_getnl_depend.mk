$(OBJDIR)/tr_getnl.o : tr_getnl.f90 
$(OBJDIR)/tr_getnl_ready.o : tr_getnl_ready.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_r4vec.o : tr_getnl_r4vec.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_putnl_lines.o : tr_putnl_lines.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_r8vec.o : tr_getnl_r8vec.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_intvec.o : tr_getnl_intvec.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_clear.o : tr_getnl_clear.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_strvals.o : tr_getnl_strvals.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_add_replace.o : tr_add_replace.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_lines.o : tr_getnl_lines.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_logvec.o : tr_getnl_logvec.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_getnl_ftext.o : tr_getnl_ftext.f90 $(OBJDIR)/tr_getnl.o 
$(OBJDIR)/tr_putnl_ftext.o : tr_putnl_ftext.f90 $(OBJDIR)/tr_getnl.o 
