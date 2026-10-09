$(OBJDIR)/tconnect_mod.o : tconnect_mod.f90 
$(OBJDIR)/tconnect_close.o : tconnect_close.f90 $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/trgs2fetch.o : trgs2fetch.f90 
$(OBJDIR)/rd_species.o : rd_species.f90 
$(OBJDIR)/rd_nspecies.o : rd_nspecies.f90 
$(OBJDIR)/tr_getnl_text.o : tr_getnl_text.f90 
$(OBJDIR)/t1mhdeq.o : t1mhdeq.f90 
$(OBJDIR)/run_mg.o : run_mg.f90 
$(OBJDIR)/trread_lun.o : trread_lun.f90 
$(OBJDIR)/t1scalar.o : t1scalar.f90 
$(OBJDIR)/mds_socket_id.o : mds_socket_id.f90 
$(OBJDIR)/t1profil.o : t1profil.f90 
$(OBJDIR)/plcfget.o : plcfget.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/rplist.o : rplist.for 
$(OBJDIR)/rpcalcd.o : rpcalcd.for 
$(OBJDIR)/trprofx.o : trprofx.for 
$(OBJDIR)/plcfget_mds.o : plcfget_mds.for 
$(OBJDIR)/rpstats.o : rpstats.for 
$(OBJDIR)/tconnect.o : tconnect.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/trinfo.o : trinfo.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/rpnumx.o : rpnumx.for 
$(OBJDIR)/tmklbl.o : tmklbl.for 
$(OBJDIR)/rpscalar.o : rpscalar.for 
$(OBJDIR)/rpgetgc.o : rpgetgc.for 
$(OBJDIR)/trfunid.o : trfunid.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/rpmgcalc.o : rpmgcalc.for 
$(OBJDIR)/tgetscal.o : tgetscal.for 
$(OBJDIR)/trscalar.o : trscalar.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/rpcalc.o : rpcalc.for 
$(OBJDIR)/tgetpath.o : tgetpath.for 
$(OBJDIR)/trprofil.o : trprofil.for $(OBJDIR)/tconnect_mod.o 
$(OBJDIR)/rplabel.o : rplabel.for 
$(OBJDIR)/rptimav.o : rptimav.for 
$(OBJDIR)/rpdims.o : rpdims.for 
$(OBJDIR)/rptime_s.o : rptime_s.for 
$(OBJDIR)/trcaps.o : trcaps.for 
$(OBJDIR)/rpbufsiz.o : rpbufsiz.for 
$(OBJDIR)/tgetprof.o : tgetprof.for 
$(OBJDIR)/rdi_ckpdens.o : rdi_ckpdens.for 
$(OBJDIR)/trread_dd.o : trread_dd.for 
$(OBJDIR)/rprofile.o : rprofile.for 
$(OBJDIR)/rptime_p.o : rptime_p.for 
$(OBJDIR)/kconnect.o : kconnect.for 
$(OBJDIR)/mdscls.o : mdscls.for 
$(OBJDIR)/rpxname.o : rpxname.for 
$(OBJDIR)/trinfg2.o : trinfg2.for 
$(OBJDIR)/rpnlist.o : rpnlist.for 
$(OBJDIR)/plcexec.o : plcexec.for 
$(OBJDIR)/ck_rplist.o : ck_rplist.for 
$(OBJDIR)/rpfixuns.o : rpfixuns.for 
$(OBJDIR)/plcfget_dd.o : plcfget_dd.for 
$(OBJDIR)/tget_rlbl.o : tget_rlbl.for 
$(OBJDIR)/rpsetgc.o : rpsetgc.for 
$(OBJDIR)/tgethost.o : tgethost.for 
$(OBJDIR)/rpmulti.o : rpmulti.for 
