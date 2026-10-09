$(OBJDIR)/trdatbuf_obj.o : trdatbuf_obj.f90 
$(OBJDIR)/trdatbuf_aux.o : trdatbuf_aux.f90 
$(OBJDIR)/tdb_static.o : tdb_static.f90 
$(OBJDIR)/tdbsub_uts.o : tdbsub_uts.f90 trdatbuf_constants.incl 
$(OBJDIR)/trdatbuf_iface.o : trdatbuf_iface.f90 trdatbuf_constants.incl $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/trdatbuf_intmod.o : trdatbuf_intmod.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/trdatbuf_module.o : trdatbuf_module.f90 $(OBJDIR)/trdatbuf_iface.o $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/rd_trdatbuf.o : rd_trdatbuf.f90 $(OBJDIR)/trdatbuf_intmod.o $(OBJDIR)/trdatbuf_module.o 
$(OBJDIR)/str2real.o : str2real.f90 
$(OBJDIR)/tclwrite.o : tclwrite.f90 
$(OBJDIR)/tdb_ccwchk.o : tdb_ccwchk.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_ebmax.o : tdb_ebmax.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_fallback.o : tdb_fallback.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_find_nmoms.o : tdb_find_nmoms.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_free_bdy.o : tdb_free_bdy.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_getmmx.o : tdb_getmmx.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_psirz.o : tdb_get_psirz.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_rzbdy.o : tdb_get_rzbdy.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_getrzfs.o : tdb_getrzfs.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_rzgrids.o : tdb_get_rzgrids.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_rzsizes.o : tdb_get_rzsizes.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_impurities.o : tdb_impurities.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_izminmax.o : tdb_izminmax.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_logchk_special.o : tdb_logchk_special.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_log_misc.o : tdb_log_misc.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_lunmsg.o : tdb_lunmsg.f90 $(OBJDIR)/tdb_static.o 
$(OBJDIR)/tdb_merge_chk.o : tdb_merge_chk.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_merge_test.o : tdb_merge_test.f90 
$(OBJDIR)/tdb_mhd_bdy.o : tdb_mhd_bdy.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_nchan_find.o : tdb_nchan_find.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_ntimes.o : tdb_ntimes.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_onoff_times.o : tdb_onoff_times.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_pelda.o : tdb_pelda.f90 $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_pwrdata_avg.o : tdb_pwrdata_avg.f90 $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_ripple0.o : tdb_ripple0.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_rmp_bdy.o : tdb_rmp_bdy.f90 $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_sawtimes.o : tdb_sawtimes.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_sevtbl.o : tdb_sevtbl.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_species.o : tdb_species.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_steerxcep.o : tdb_steerxcep.f90 
$(OBJDIR)/tdbsub.o : tdbsub.f90 $(OBJDIR)/tdbsub_uts.o 
$(OBJDIR)/tdbsub_onoff.o : tdbsub_onoff.f90 $(OBJDIR)/tdbsub_uts.o 
$(OBJDIR)/tdbsub_slookup.o : tdbsub_slookup.f90 $(OBJDIR)/tdbsub_uts.o 
$(OBJDIR)/tdb_tlims.o : tdb_tlims.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_xisymp.o : tdb_xisymp.f90 $(OBJDIR)/trdatbuf_iface.o 
$(OBJDIR)/trdatbuf_intcopy.o : trdatbuf_intcopy.f90 
$(OBJDIR)/trdatbuf_intmod_sub1.o : trdatbuf_intmod_sub1.f90 
$(OBJDIR)/trdatbuf_intmod_sub2.o : trdatbuf_intmod_sub2.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/trdatbuf_intmod_sub3.o : trdatbuf_intmod_sub3.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/trdatbuf_io.o : trdatbuf_io.f90 $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/trdatbuf_iosub.o : trdatbuf_iosub.f90 
$(OBJDIR)/tree_close.o : tree_close.f90 
$(OBJDIR)/trnode_child.o : trnode_child.f90 
$(OBJDIR)/wr_cstringarr.o : wr_cstringarr.f90 
$(OBJDIR)/wr_cstring.o : wr_cstring.f90 
$(OBJDIR)/wr_doublearr.o : wr_doublearr.f90 
$(OBJDIR)/wr_longarr.o : wr_longarr.f90 
$(OBJDIR)/wr_trdatbuf.o : wr_trdatbuf.f90 $(OBJDIR)/trdatbuf_intmod.o $(OBJDIR)/trdatbuf_module.o 
$(OBJDIR)/iadci.o : iadci.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/int2d.o : int2d.for 
$(OBJDIR)/momind3.o : momind3.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/momind.o : momind.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/phdata_io.o : phdata_io.for 
$(OBJDIR)/rd_cstringarr.o : rd_cstringarr.for 
$(OBJDIR)/rd_doublearr.o : rd_doublearr.for 
$(OBJDIR)/rd_longarr.o : rd_longarr.for 
$(OBJDIR)/tdb_chkprof.o : tdb_chkprof.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_datalo.o : tdb_datalo.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_execsym.o : tdb_execsym.for $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_fintrp1.o : tdb_fintrp1.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_rpld2.o : tdb_get_rpld2.for $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_get_rpld.o : tdb_get_rpld.for $(OBJDIR)/tdbsub_uts.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_maksym.o : tdb_maksym.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_pfsetr.o : tdb_pfsetr.for $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_presym.o : tdb_presym.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_prfins.o : tdb_prfins.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_prfinv.o : tdb_prfinv.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_prflax.o : tdb_prflax.for 
$(OBJDIR)/tdb_prfnxp.o : tdb_prfnxp.for 
$(OBJDIR)/tdb_prfnxx.o : tdb_prfnxx.for 
$(OBJDIR)/tdb_profin.o : tdb_profin.for $(OBJDIR)/trdatbuf_iface.o $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_profli.o : tdb_profli.for $(OBJDIR)/trdatbuf_iface.o $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_rovera.o : tdb_rovera.for $(OBJDIR)/trdatbuf_aux.o 
$(OBJDIR)/tdbsub_ck2a.o : tdbsub_ck2a.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_symini.o : tdb_symini.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_symxmap.o : tdb_symxmap.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_unmap.o : tdb_unmap.for $(OBJDIR)/trdatbuf_aux.o $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_workalo.o : tdb_workalo.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_xygrid.o : tdb_xygrid.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_xyprof.o : tdb_xyprof.for $(OBJDIR)/trdatbuf_obj.o 
$(OBJDIR)/tdb_zyget.o : tdb_zyget.for $(OBJDIR)/trdatbuf_obj.o 
