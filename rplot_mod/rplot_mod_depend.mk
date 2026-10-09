$(OBJDIR)/cplotr_mod.o : cplotr_mod.f90 
$(OBJDIR)/extrac2_mod.o : extrac2_mod.f90 
$(OBJDIR)/mfblok_mod.o : mfblok_mod.f90 $(OBJDIR)/cplotr_mod.o 
$(OBJDIR)/nltrdat_mod.o : nltrdat_mod.f90 
$(OBJDIR)/plfmpa_mod.o : plfmpa_mod.f90 $(OBJDIR)/cplotr_mod.o 
$(OBJDIR)/rpcalc_mod.o : rpcalc_mod.f90 
