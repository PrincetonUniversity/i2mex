$(OBJDIR)/trx_bxtr_options.o : trx_bxtr_options.f90 
$(OBJDIR)/trx_module.o : trx_module.f90 
$(OBJDIR)/trxplib_ps_options.o : trxplib_ps_options.f90 
$(OBJDIR)/trx_toray_in_search.o : trx_toray_in_search.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_spec_lbl.o : trx_spec_lbl.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_prof.o : trx_prof.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_ustore.o : trx_ustore.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_scal.o : trx_scal.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_genray_in_search.o : trx_genray_in_search.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trxplib_ps_write1.o : trxplib_ps_write1.f90 $(OBJDIR)/trxplib_ps_options.o 
$(OBJDIR)/trx_getlnlim.o : trx_getlnlim.f90 
$(OBJDIR)/trx_time.o : trx_time.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trxplib_ps_xplasma_ini.o : trxplib_ps_xplasma_ini.f90 $(OBJDIR)/trxplib_ps_options.o 
$(OBJDIR)/trx_ntimes.o : trx_ntimes.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_spec_edge.o : trx_spec_edge.f90 
$(OBJDIR)/trx_bxtr_setget.o : trx_bxtr_setget.f90 $(OBJDIR)/trx_bxtr_options.o 
$(OBJDIR)/trx_msgs.o : trx_msgs.f90 
$(OBJDIR)/trx_calc.o : trx_calc.f90 
$(OBJDIR)/trx_label.o : trx_label.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_extract_path.o : trx_extract_path.f90 
$(OBJDIR)/trx_wr_protimes.o : trx_wr_protimes.f90 $(OBJDIR)/trxplib_ps_options.o $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_getnlims.o : trx_getnlims.f90 
$(OBJDIR)/trxplib_ps_mdescr.o : trxplib_ps_mdescr.f90 $(OBJDIR)/trxplib_ps_options.o 
$(OBJDIR)/trx_ready.o : trx_ready.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_uget.o : trx_uget.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_nspec.o : trx_nspec.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trxplib_ps_connect.o : trxplib_ps_connect.f90 $(OBJDIR)/trxplib_ps_options.o 
$(OBJDIR)/trx_lintrans.o : trx_lintrans.f90 
$(OBJDIR)/trx_mhd.o : trx_mhd.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_nb_sublist.o : trx_nb_sublist.f90 
$(OBJDIR)/trxplib_ps_io.o : trxplib_ps_io.f90 $(OBJDIR)/trxplib_ps_options.o 
$(OBJDIR)/trx_psirz_load.o : trx_psirz_load.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_gen_state.o : trx_gen_state.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_rzbc.o : trx_rzbc.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_spec_prof.o : trx_spec_prof.f90 
$(OBJDIR)/trx_bxtr.o : trx_bxtr.f90 $(OBJDIR)/trx_bxtr_options.o $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_wall.o : trx_wall.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_psi0.o : trx_psi0.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_freebdy.o : trx_freebdy.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_trdatbuf_connect.o : trx_trdatbuf_connect.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_imas.o : trx_imas.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_gtime.o : trx_gtime.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_gen_ech_files.o : trx_gen_ech_files.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_getcrlim.o : trx_getcrlim.f90 
$(OBJDIR)/trx_tlims.o : trx_tlims.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_connect.o : trx_connect.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_chk_saw.o : trx_chk_saw.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_gafit_in_search.o : trx_gafit_in_search.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_init.o : trx_init.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_mks_conv.o : trx_mks_conv.f90 
$(OBJDIR)/trx_wr_stimes.o : trx_wr_stimes.f90 $(OBJDIR)/trx_module.o 
$(OBJDIR)/trx_gen_lh_files.o : trx_gen_lh_files.f90 $(OBJDIR)/trx_module.o 
