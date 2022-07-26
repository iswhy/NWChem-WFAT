# $Id$
#
#      OBJ = dft_scf.o dft_s2_value.o dft_main0d.o
#      OBJ_OPTIMIZE = \
#	diis_driver.o	diis_bsolve.o	\
#	diis_core.o	diis_bld12.o diis_hamwgt.o 	\
#	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
#	dft_invdiag.o	dft_init_guess.o \
#	dft_mxovl.o  dft_mxspin_ovlp.o
#

	OBJ_OPTIMIZE = \
	dft_dencvg.o \
	dft_densm.o \
	dft_invdiag.o \
	dft_mxovl.o \
	dft_mxspin_ovlp.o \
	dft_s2_value.o \
	dft_scfcvg.o \
	diis_bld12.o \
	diis_bsolve.o \
	diis_core.o \
	diis_driver.o \
	diis_hamwgt.o \
	dft_fdist_init.o \
	cdft_util.o
	OBJ=\
	dft_swapab.o \
	dft_main0d.o \
	dft_scf.o \
	print_irrep_occp.o \
	dft_fockbld.o \
	dft_scfutil.o \
	dft_symbr.o \
	dft_pstat.o\
	fukui.o \
	dft_canorg.o \
	dft_diagn.o \
	dft_dump_info.o \
	dft_fant.o\
	dft_printinfo.o\
	dft_frac_mo.o \
	dft_nwmolden.o \
	wfat_input.o \
	wfat_process_input.o \
	wfat_task.o \
	wfat_prepare_task.o \
	wfat_init.o \
	wfat_init_test.o \
	wfat_calc_intgl.o \
	wfat_calc_structuref.o \
	wfat_calc_ionyield.o \
	wfat_calc_tmd.o \
	wfat_calc_ti.o \
	wfat_calc_td.o \
	wfat_field_factor.o \
	wfat_celf.o \
	wfat_chypgeom.o \
	wfat_density.o \
	wfat_eval_basis.o \
	wfat_grid.o \
	wfat_lebedev.o \
	wfat_omega.o \
	wfat_potential_calc.o \
	wfat_potential_calc_asy.o \
	wfat_potential.o \
	wfat_potential_asy.o \
	wfat_potential_rs.o \
	wfat_tools.o \
	wfat_par_tools.o \
	wfat_test.o \
	wfat_test2.o \
	wfat_numdens.o \
	wfat_orb_property.o \
	wfat_eval_mo.o \
	wfat_eval_laser.o \
	wfat_print.o \
	wfat_print_sf.o \
	wfat_print_rt.o \
	wfat_print_iy.o \
	wfat_print_tm.o \
	wfat_print_tmar.o \
	wfat_print_info.o \
	wfat_print_dens.o \
	wfat_util.o \
	wfat_orb_mulliken.o \
	wfat_opvec_coef.o \
	wfat_get_opmat.o \
	wfat_get_fock.o \
	wfat_perturb.o \
	wfat_perturb_bw.o \
	wfat_perturb_rs.o \
	wfat_perturb_analyse.o \
	wfat_perturb_print.o \
	wfat_partial_integ.o \
	mewfat_task.o \
	mewfat_init.o \
	mewfat_print_info.o \
	mewfat_need_input_v2.o \
	mewfat_channel_property.o \
	mewfat_calc_intgl.o \
	mewfat_io_intgl.o \
	mewfat_calc_detcoef.o \
	mewfat_density.o \
	mewfat_dfterg.o \
	mewfat_eval_dyson.o \
	mewfat_potential.o \
	mewfat_potential_calc.o \
	mewfat_dyson_print.o \
	imam_test.o

     LIBRARY = libnwdft.a
ifndef OLD_GA
_USE_SCALAPACK = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_SCALAPACK\ 1/ {print "Y"}')

ifdef _USE_SCALAPACK
ifeq ($(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk '/HAVE_ELPA\ 1/ {print "Y"}'),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
GOTGAPZHEEVD= $(shell [ -e ${NWCHEM_TOP}/src/tools/build/global/src/scalapack.o ] &&  nm ${NWCHEM_TOP}/src/tools/build/global/src/scalapack.o|egrep T\ ga_pzheevd_ |awk '/ga_pzheevd/ { print "Y"}')
ifeq ($(GOTGAPZHEEVD),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_PDSYEVR\ 1/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
endif

LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_mxovl.F diis_bld12.F diis_bsolve.F dft_densm.F dft_swapab.F dft_canorg.F cdft_util.F

HEADERS = dftps.fh dft_fdist.fh dft_conv.fh

include ../../config/makefile.h
include ../../config/makelib.h
