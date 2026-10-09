$(OBJDIR)/catfus_mod.o : catfus_mod.f90 
$(OBJDIR)/lincir_contour.o : lincir_contour.f90 
$(OBJDIR)/namerefs.o : namerefs.f90 
$(OBJDIR)/periodic_table.o : periodic_table.f90 
$(OBJDIR)/plasma_hash_code.o : plasma_hash_code.f90 
$(OBJDIR)/r8bsmoo_mod.o : r8bsmoo_mod.f90 
$(OBJDIR)/submsg.o : submsg.f90 
$(OBJDIR)/ode.o : ode.f90 
$(OBJDIR)/atoi.o : atoi.f90 
$(OBJDIR)/c8dft.o : c8dft.f90 
$(OBJDIR)/catfus.o : catfus.f90 $(OBJDIR)/catfus_mod.o 
$(OBJDIR)/catfuslis.o : catfuslis.f90 $(OBJDIR)/catfus_mod.o 
$(OBJDIR)/chk_extrapolate.o : chk_extrapolate.f90 
$(OBJDIR)/chk_interpolate.o : chk_interpolate.f90 
$(OBJDIR)/ck_sepseg.o : ck_sepseg.f90 
$(OBJDIR)/cloge_r8_fcn.o : cloge_r8_fcn.f90 
$(OBJDIR)/coulog.o : coulog.f90 
$(OBJDIR)/cperiodic_table.o : cperiodic_table.f90 $(OBJDIR)/periodic_table.o 
$(OBJDIR)/ech_angles.o : ech_angles.f90 
$(OBJDIR)/erf.o : erf.f90 
$(OBJDIR)/expi.o : expi.f90 
$(OBJDIR)/fhermitd.o : fhermitd.f90 
$(OBJDIR)/fhermite.o : fhermite.f90 
$(OBJDIR)/filfas.o : filfas.f90 
$(OBJDIR)/flin1.o : flin1.f90 
$(OBJDIR)/flint.o : flint.f90 
$(OBJDIR)/fpolar.o : fpolar.f90 
$(OBJDIR)/f_upwind0.o : f_upwind0.f90 
$(OBJDIR)/geotrn.o : geotrn.f90 
$(OBJDIR)/get_j_or_jdotb.o : get_j_or_jdotb.f90 
$(OBJDIR)/getlog2int.o : getlog2int.f90 
$(OBJDIR)/getmpa.o : getmpa.f90 
$(OBJDIR)/get_namerefs.o : get_namerefs.f90 $(OBJDIR)/namerefs.o 
$(OBJDIR)/hhb_resis.o : hhb_resis.f90 
$(OBJDIR)/integ_wts.o : integ_wts.f90 
$(OBJDIR)/irns.o : irns.f90 
$(OBJDIR)/jacevala.o : jacevala.f90 
$(OBJDIR)/jaceval.o : jaceval.f90 
$(OBJDIR)/lkupr.o : lkupr.f90 
$(OBJDIR)/lkupr_r8.o : lkupr_r8.f90 
$(OBJDIR)/lokjaca.o : lokjaca.f90 
$(OBJDIR)/lokjac.o : lokjac.f90 
$(OBJDIR)/makedens_gneg.o : makedens_gneg.f90 
$(OBJDIR)/mkquin2.o : mkquin2.f90 
$(OBJDIR)/namels.o : namels.f90 $(OBJDIR)/periodic_table.o 
$(OBJDIR)/namerefs_uts.o : namerefs_uts.f90 $(OBJDIR)/namerefs.o 
$(OBJDIR)/nnmod.o : nnmod.f90 
$(OBJDIR)/parsid.o : parsid.f90 
$(OBJDIR)/plcentr.o : plcentr.f90 
$(OBJDIR)/ps_namel_strip.o : ps_namel_strip.f90 
$(OBJDIR)/ps_namrd_ilist_chk.o : ps_namrd_ilist_chk.f90 
$(OBJDIR)/ps_namrd_slist_chk.o : ps_namrd_slist_chk.f90 
$(OBJDIR)/qhermitd.o : qhermitd.f90 
$(OBJDIR)/r4dft.o : r4dft.f90 
$(OBJDIR)/r4fftsc.o : r4fftsc.f90 
$(OBJDIR)/r4rfft_sub.o : r4rfft_sub.f90 
$(OBJDIR)/r8_blims_mmx.o : r8_blims_mmx.f90 
$(OBJDIR)/r8bsmoo.o : r8bsmoo.f90 $(OBJDIR)/r8bsmoo_mod.o 
$(OBJDIR)/r8bsmooh.o : r8bsmooh.f90 $(OBJDIR)/r8bsmoo_mod.o 
$(OBJDIR)/r8bsmoo_init.o : r8bsmoo_init.f90 $(OBJDIR)/r8bsmoo_mod.o 
$(OBJDIR)/r8_coulog.o : r8_coulog.f90 
$(OBJDIR)/r8fftsc.o : r8fftsc.f90 
$(OBJDIR)/r8_filfas.o : r8_filfas.f90 
$(OBJDIR)/r8_fpolar.o : r8_fpolar.f90 
$(OBJDIR)/r8_gausse.o : r8_gausse.f90 
$(OBJDIR)/r8_getmpa.o : r8_getmpa.f90 
$(OBJDIR)/r8jacevala.o : r8jacevala.f90 
$(OBJDIR)/r8jaceval.o : r8jaceval.f90 
$(OBJDIR)/r8_plcentr.o : r8_plcentr.f90 
$(OBJDIR)/r8rfft_sub.o : r8rfft_sub.f90 
$(OBJDIR)/r8ryevala.o : r8ryevala.f90 
$(OBJDIR)/r8ryeval.o : r8ryeval.f90 
$(OBJDIR)/r8_sevali.o : r8_sevali.f90 
$(OBJDIR)/r8_shafnl.o : r8_shafnl.f90 
$(OBJDIR)/r8sincos.o : r8sincos.f90 
$(OBJDIR)/r8splcrda.o : r8splcrda.f90 
$(OBJDIR)/r8splcrd.o : r8splcrd.f90 
$(OBJDIR)/r8_trirad.o : r8_trirad.f90 
$(OBJDIR)/r8_xintzb.o : r8_xintzb.f90 
$(OBJDIR)/r8_zeroin.o : r8_zeroin.f90 
$(OBJDIR)/rf_egridtr.o : rf_egridtr.f90 
$(OBJDIR)/rfx_antcen.o : rfx_antcen.f90 
$(OBJDIR)/ryevala.o : ryevala.f90 
$(OBJDIR)/ryeval.o : ryeval.f90 
$(OBJDIR)/sevali.o : sevali.f90 
$(OBJDIR)/shafnl.o : shafnl.f90 
$(OBJDIR)/sincos.o : sincos.f90 
$(OBJDIR)/splcrda.o : splcrda.f90 
$(OBJDIR)/splcrd.o : splcrd.f90 
$(OBJDIR)/split_imp.o : split_imp.f90 
$(OBJDIR)/surfgeo.o : surfgeo.f90 
$(OBJDIR)/tkbndr.o : tkbndr.f90 tkbparm.inc 
$(OBJDIR)/tkbrad.o : tkbrad.f90 tkbparm.inc 
$(OBJDIR)/trdate.o : trdate.f90 
$(OBJDIR)/tr_dsort.o : tr_dsort.f90 
$(OBJDIR)/trigdr.o : trigdr.f90 
$(OBJDIR)/tr_ir8sort.o : tr_ir8sort.f90 
$(OBJDIR)/tr_qelg.o : tr_qelg.f90 
$(OBJDIR)/tr_qpsrt.o : tr_qpsrt.f90 
$(OBJDIR)/tr_r8qags.o : tr_r8qags.f90 
$(OBJDIR)/tr_r8qelg.o : tr_r8qelg.f90 
$(OBJDIR)/tr_r8qpsrt.o : tr_r8qpsrt.f90 
$(OBJDIR)/tr_ssort.o : tr_ssort.f90 
$(OBJDIR)/unflatten.o : unflatten.f90 
$(OBJDIR)/vnewton.o : vnewton.f90 
$(OBJDIR)/xintrp.o : xintrp.f90 
$(OBJDIR)/xintz0.o : xintz0.f90 
$(OBJDIR)/xintzb.o : xintzb.f90 
$(OBJDIR)/xintzc.o : xintzc.f90 
$(OBJDIR)/zeroin.o : zeroin.f90 
$(OBJDIR)/zkbolt.o : zkbolt.f90 
$(OBJDIR)/zridderx.o : zridderx.f90 
$(OBJDIR)/zriddery.o : zriddery.f90 
