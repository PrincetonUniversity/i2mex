$(OBJDIR)/fluxav.o : fluxav.f90 
$(OBJDIR)/fluxav_balloc.o : fluxav_balloc.f90 $(OBJDIR)/fluxav.o 
$(OBJDIR)/fluxav_brk_init.o : fluxav_brk_init.f90 $(OBJDIR)/fluxav.o 
$(OBJDIR)/fluxav_brk_add.o : fluxav_brk_add.f90 $(OBJDIR)/fluxav.o 
$(OBJDIR)/fluxav_opts.o : fluxav_opts.f90 $(OBJDIR)/fluxav.o 
$(OBJDIR)/fluxav_brk_add0.o : fluxav_brk_add0.f90 
$(OBJDIR)/dqng_frhoth.o : dqng_frhoth.for 
$(OBJDIR)/dqng_fth.o : dqng_fth.for 
