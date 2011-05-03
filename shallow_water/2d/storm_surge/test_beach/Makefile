#
# Local makefile for test_beach
#

# Path to storm surge specific source


# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = geoclaw
CLAW_setrun_file = setrun.py
CLAW_setplot_file = setplot.py
CLAW_OUTDIR = _output
CLAW_PLOTDIR = _plots

GEOCLAW=$(CLAW)/geoclaw/2d/lib
AMRCLAW=$(CLAW)/amrclaw/2d/lib
SINGLE_LAYER=../single_layer
HURRICANE_SRC=../src/

# Flagging
FC ?= gfortran
FFLAGS += -I$(GEOCLAW) -I$(AMRCLAW)
UNAME = $(shell uname)
ifeq ($(UNAME),Darwin)
	ifndef LFLAGS
		LFLAGS = -framework veclib
	else
		LFLAGS += -framework veclib
	endif
endif
ifeq ($(UNAME),Linux)
	ifndef LFLAGS
		LFLAGS = -llapack
	else
		LFLAGS += -llapack
	endif
endif

# ---------------------------------
# List of sources for this program:
# ---------------------------------
    
CLAW_SOURCES = \
    $(HURRICANE_SRC)/multilayer_module.f90 \
	$(HURRICANE_SRC)/hurricane.f90 \
	$(SINGLE_LAYER)/qinit_geo.f90 \
	$(HURRICANE_SRC)/setprob.f90 \
	$(HURRICANE_SRC)/b4step2_geo.f \
	$(HURRICANE_SRC)/valout_geo.f \
	$(HURRICANE_SRC)/setaux_geo.f \
	$(SINGLE_LAYER)/flag2refine_geo.f \
	$(SINGLE_LAYER)/src2_geo.f90 \
	$(SINGLE_LAYER)/src1d_geo.f90

CLAW_LIBSOURCES = \
	$(GEOCLAW)/geoclaw_mod.f90 \
	$(GEOCLAW)/rpn2ez_geo.f \
	$(GEOCLAW)/riemannsolvers_geo.f \
	$(GEOCLAW)/rpt2_geo.f \
	$(GEOCLAW)/topo_mod.f90 \
	$(GEOCLAW)/dtopo_mod.f90 \
	$(GEOCLAW)/readqinit_geo.f \
	$(GEOCLAW)/movetopo_geo.f \
	$(GEOCLAW)/cellgridintegrate_geo.f \
	$(GEOCLAW)/topointegral_geo.f \
	$(GEOCLAW)/bilinearintegral_geo.f \
	$(GEOCLAW)/stepgrid_geo.f \
	$(GEOCLAW)/fgridinterp_geo.f \
	$(GEOCLAW)/fgridout_geo.f \
	$(GEOCLAW)/step2_geo.f \
	$(GEOCLAW)/flux2fw_geo.f \
	$(GEOCLAW)/qad_geo.f \
	$(GEOCLAW)/filval_geo.f \
	$(GEOCLAW)/filpatch_geo.f \
	$(GEOCLAW)/bc2amr_geo.f \
	$(GEOCLAW)/update_geo.f \
	$(GEOCLAW)/amr2ez_geo.f \
	$(GEOCLAW)/setqinit_geo.f \
	$(GEOCLAW)/setregions_geo.f \
	$(GEOCLAW)/setgauges_geo.f \
	$(GEOCLAW)/setfixedgrids_geo.f \
	$(GEOCLAW)/allowflag_geo.f  \
	$(GEOCLAW)/dumpgauge_geo.f \
	$(GEOCLAW)/upbnd_geo.f  \
	$(GEOCLAW)/tick_geo.f \
	$(GEOCLAW)/setgrd_geo.f \
	$(GEOCLAW)/gfixup_geo.f \
	$(GEOCLAW)/ginit_geo.f \
	$(GEOCLAW)/getmaxspeed_geo.f \
	$(AMRCLAW)/prefilp.f \
	$(AMRCLAW)/quick_sort1.f \
	$(AMRCLAW)/advanc.f \
	$(AMRCLAW)/estdt.f \
	$(AMRCLAW)/bound.f \
	$(AMRCLAW)/check4nans.f \
	$(AMRCLAW)/spest.f  \
	$(AMRCLAW)/igetsp.f \
	$(AMRCLAW)/reclam.f \
	$(AMRCLAW)/birect.f \
	$(AMRCLAW)/cleanup.f \
	$(AMRCLAW)/colate.f \
	$(AMRCLAW)/errest.f \
	$(AMRCLAW)/bufnst.f \
	$(AMRCLAW)/errf1.f \
	$(AMRCLAW)/fixcapaq.f \
	$(AMRCLAW)/flglvl.f \
	$(AMRCLAW)/fluxad.f \
	$(AMRCLAW)/fluxsv.f \
	$(AMRCLAW)/grdfit.f \
	$(AMRCLAW)/intfil.f \
	$(AMRCLAW)/moment.f \
	$(AMRCLAW)/nestck.f \
	$(AMRCLAW)/prepc.f \
	$(AMRCLAW)/prepf.f \
	$(AMRCLAW)/projec.f \
	$(AMRCLAW)/signs.f \
	$(AMRCLAW)/findcut.f \
	$(AMRCLAW)/smartbis.f \
	$(AMRCLAW)/putnod.f \
	$(AMRCLAW)/putsp.f \
	$(AMRCLAW)/regrid.f \
	$(AMRCLAW)/setuse.f \
	$(AMRCLAW)/stst1.f \
	$(AMRCLAW)/trimbd.f \
	$(AMRCLAW)/nodget.f \
	$(AMRCLAW)/basic.f \
	$(AMRCLAW)/outval.f \
	$(AMRCLAW)/copysol.f \
	$(AMRCLAW)/outvar.f \
	$(AMRCLAW)/outmsh.f \
	$(AMRCLAW)/outtre.f \
	$(AMRCLAW)/domain.f  \
	$(AMRCLAW)/cellave.f \
	$(AMRCLAW)/fdisc.f \
	$(AMRCLAW)/fss.f \
	$(AMRCLAW)/zeroin.f \
	$(AMRCLAW)/setflags.f \
	$(AMRCLAW)/shiftset.f \
	$(AMRCLAW)/conck.f \
	$(AMRCLAW)/domshrink.f \
	$(AMRCLAW)/domprep.f \
	$(AMRCLAW)/domup.f \
	$(AMRCLAW)/domcopy.f \
	$(AMRCLAW)/coarsen.f \
	$(AMRCLAW)/auxcoarsen.f \
	$(AMRCLAW)/intcopy.f \
	$(AMRCLAW)/preintcopy.f \
	$(AMRCLAW)/icall.f \
	$(AMRCLAW)/preicall.f \
	$(AMRCLAW)/inlinelimiter.f \
	$(AMRCLAW)/cstore.f \
	$(AMRCLAW)/saveqc.f \
	$(AMRCLAW)/check.f \
	$(AMRCLAW)/restrt.f \
	$(AMRCLAW)/opendatafile.f \
	$(AMRCLAW)/init_alloc.f90 \
	$(AMRCLAW)/restrt_alloc.f90 \
	$(AMRCLAW)/resize_alloc.f90 \


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAW)/util/Makefile.common
	
### DO NOT remove this line - make depends on it ###
    
