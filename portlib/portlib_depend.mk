$(OBJDIR)/mpi_proc_data.o : mpi_proc_data.f90 
$(OBJDIR)/mpi_env_mod.o : mpi_env_mod.f90 $(OBJDIR)/mpi_proc_data.o 
$(OBJDIR)/logmod.o : logmod.f90 
$(OBJDIR)/execsystem.o : execsystem.f90 $(OBJDIR)/mpi_proc_data.o $(OBJDIR)/logmod.o 
$(OBJDIR)/mmpi.o : mmpi.f90 $(OBJDIR)/logmod.o 
$(OBJDIR)/mpi_mkgroup.o : mpi_mkgroup.f90 
$(OBJDIR)/mpi_sset_env.o : mpi_sset_env.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/mpi_sget_env.o : mpi_sget_env.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/mpi_printenv.o : mpi_printenv.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/text_copy.o : text_copy.f90 
$(OBJDIR)/parcopy.o : parcopy.f90 
$(OBJDIR)/mpi_env_scomm.o : mpi_env_scomm.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/mpi_share_env.o : mpi_share_env.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/mkdir.o : mkdir.f90 
$(OBJDIR)/int_bcast.o : int_bcast.f90 
$(OBJDIR)/mpi_procname.o : mpi_procname.f90 
$(OBJDIR)/errset_mpi.o : errset_mpi.f90 $(OBJDIR)/mpi_env_mod.o 
$(OBJDIR)/system_call_echo.o : system_call_echo.f90 
$(OBJDIR)/find_io_unit.o : find_io_unit.f90 
$(OBJDIR)/parse.o : parse.f90 
$(OBJDIR)/wclock_diff.o : wclock_diff.f90 
$(OBJDIR)/is_ascii.o : is_ascii.f90 
$(OBJDIR)/c9date.o : c9date.f90 
$(OBJDIR)/wait_for_file.o : wait_for_file.f90 
$(OBJDIR)/get_fortran_const.o : get_fortran_const.f90 
$(OBJDIR)/f90_link1.o : f90_link1.for 
$(OBJDIR)/genopen.o : genopen.for 
$(OBJDIR)/str_pad.o : str_pad.for 
$(OBJDIR)/ierfnof.o : ierfnof.for 
$(OBJDIR)/bad_exit.o : bad_exit.for 
$(OBJDIR)/showfile.o : showfile.for 
$(OBJDIR)/errmsg_exit.o : errmsg_exit.for 
$(OBJDIR)/max_reclen.o : max_reclen.for 
$(OBJDIR)/sget_pid_str.o : sget_pid_str.for 
$(OBJDIR)/datcpu_u.o : datcpu_u.for 
$(OBJDIR)/fgetline.o : fgetline.for 
$(OBJDIR)/nblkfac.o : nblkfac.for 
$(OBJDIR)/sset_cwd.o : sset_cwd.for 
$(OBJDIR)/sset_env.o : sset_env.for 
$(OBJDIR)/good_exit.o : good_exit.for 
$(OBJDIR)/tmpfile.o : tmpfile.for 
$(OBJDIR)/isbatch.o : isbatch.for 
$(OBJDIR)/cftio.o : cftio.for 
$(OBJDIR)/fcopy.o : fcopy.for 
$(OBJDIR)/ufilnam.o : ufilnam.for 
$(OBJDIR)/frename.o : frename.for 
$(OBJDIR)/uupper.o : uupper.for 
$(OBJDIR)/cstring.o : cstring.for 
$(OBJDIR)/showdefl.o : showdefl.for 
$(OBJDIR)/lxasc.o : lxasc.for 
$(OBJDIR)/fdelete.o : fdelete.for 
$(OBJDIR)/sget_dsk.o : sget_dsk.for 
$(OBJDIR)/fclean_dir.o : fclean_dir.for 
$(OBJDIR)/utrnlog.o : utrnlog.for 
$(OBJDIR)/fwc_delete.o : fwc_delete.for 
$(OBJDIR)/term_out.o : term_out.for 
$(OBJDIR)/ulower.o : ulower.for 
$(OBJDIR)/tempfile.o : tempfile.for 
$(OBJDIR)/cvt_days.o : cvt_days.for 
$(OBJDIR)/str_length.o : str_length.for 
$(OBJDIR)/iopen.o : iopen.for 
$(OBJDIR)/wall_seconds.o : wall_seconds.for 
$(OBJDIR)/err_end.o : err_end.for 
$(OBJDIR)/gmkdir.o : gmkdir.for 
$(OBJDIR)/ugetlog.o : ugetlog.for 
$(OBJDIR)/term_input.o : term_input.for 
