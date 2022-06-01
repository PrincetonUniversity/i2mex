$(OBJDIR)/cont_mod.o : cont_mod.f90 
$(OBJDIR)/freeqbe_mod.o : freeqbe_mod.f90 
$(OBJDIR)/i2mex_mod.o : i2mex_mod.f90 
$(OBJDIR)/delta.o : delta.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/contour.o : contour.f90 $(OBJDIR)/cont_mod.o 
$(OBJDIR)/eqm_ball.o : eqm_ball.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getx.o : getx.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/metric.o : metric.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/poloidalremap.o : poloidalremap.f90 cubinterp.h 
$(OBJDIR)/input_geqdsk.o : input_geqdsk.f90 $(OBJDIR)/cont_mod.o $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/i2mex_free.o : i2mex_free.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getg.o : getg.f90 cubinterp.h $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/glasser.o : glasser.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/save.o : save.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getj.o : getj.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/integral.o : integral.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/scaleg.o : scaleg.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/i2mex_init.o : i2mex_init.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/parameters.o : parameters.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/scaleq.o : scaleq.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getboozer.o : getboozer.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/gets.o : gets.f90 cubinterp.h $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getdelstar.o : getdelstar.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/integrate1d.o : integrate1d.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/ratsurf.o : ratsurf.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/i2mex_decimate.o : i2mex_decimate.f90 
$(OBJDIR)/graphics.o : graphics.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/input.o : input.f90 cubinterp.h $(OBJDIR)/cont_mod.o $(OBJDIR)/freeqbe_mod.o $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/i2mex_error.o : i2mex_error.f90 
$(OBJDIR)/interp1d.o : interp1d.f90 
$(OBJDIR)/getz.o : getz.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getp.o : getp.f90 cubinterp.h $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/i2mex_theta_orient.o : i2mex_theta_orient.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/toaxis.o : toaxis.f90 cubinterp.h $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getq.o : getq.f90 cubinterp.h $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/getphi.o : getphi.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/gausspoints.o : gausspoints.f90 
$(OBJDIR)/output.o : output.f90 $(OBJDIR)/i2mex_mod.o 
$(OBJDIR)/original.o : original.f90 $(OBJDIR)/i2mex_mod.o 
