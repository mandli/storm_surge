#
# Local makefile for test_beach
#

# Path to storm surge specific source
MULTILAYER_SRC = ../src

# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = geoclaw
CLAW_setrun_file = setrun.py
CLAW_setplot_file = setplot.py
CLAW_OUTDIR = $(DATA_PATH)/storm_surge/new_bathy/_output
CLAW_PLOTDIR = $(DATA_PATH)/storm_surge/new_bathy/_plots

# Clawpack libraries to be used:
CLAW_LIB = $(CLAW)/amrclaw/2d/lib
GEOLIB = $(CLAW)/geoclaw/2d/lib

# Flagging
FC = gfortran44
FFLAGS += -I$(GEOLIB)
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	ifndef LFLAGS
		LFLAGS = -framework veclib
	else
		LFLAGS += -framework veclib
	endif
endif
ifeq ($(UNAME), Linux)
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
	$(MULTILAYER_SRC)/multilayer_module.f90 \
	$(MULTILAYER_SRC)/hurricane.f90 \
	./setaux_geo.f \
	$(MULTILAYER_SRC)/b4step2_geo.f \
	$(MULTILAYER_SRC)/setprob.f90 \
	$(MULTILAYER_SRC)/src2_geo.f90 \
	$(MULTILAYER_SRC)/src1d_geo.f90 \
	$(MULTILAYER_SRC)/qinit.f90 \
	$(MULTILAYER_SRC)/valout_geo.f \
	$(MULTILAYER_SRC)/flux2fw.f \
	$(MULTILAYER_SRC)/flag2refine_geo.f \
	$(MULTILAYER_SRC)/rpn2_dry_multilayer.f90 \
	$(MULTILAYER_SRC)/rpt2_dry_multilayer.f90 \
	$(MULTILAYER_SRC)/multilayer_eigen.f90 \
	$(MULTILAYER_SRC)/tick.f \
	$(MULTILAYER_SRC)/qad_geo.f \
	$(MULTILAYER_SRC)/stepgrid_geo.f \
	$(CLAW_LIB)/filval.f \
	$(CLAW_LIB)/filpatch.f \
	$(CLAW_LIB)/update.f \
	$(CLAW_LIB)/upbnd.f

DATA_FILES = hurricane.data topo.data multilayer.data
	
CLAW_LIBSOURCES = \
	$(GEOLIB)/geoclaw_mod.f90 \
	$(GEOLIB)/riemannsolvers_geo.f \
	$(GEOLIB)/topo_mod.f90 \
	$(GEOLIB)/dtopo_mod.f90 \
	$(GEOLIB)/readqinit_geo.f \
	$(GEOLIB)/movetopo_geo.f \
    $(GEOLIB)/cellgridintegrate_geo.f \
	$(GEOLIB)/topointegral_geo.f \
	$(GEOLIB)/bilinearintegral_geo.f \
	$(GEOLIB)/fgridinterp_geo.f \
	$(GEOLIB)/fgridout_geo.f \
	$(GEOLIB)/step2_geo.f \
	$(GEOLIB)/bc2amr_geo.f \
	$(GEOLIB)/amr2ez_geo.f \
	$(GEOLIB)/setqinit_geo.f \
	$(GEOLIB)/setregions_geo.f \
	$(GEOLIB)/setgauges_geo.f \
	$(GEOLIB)/setfixedgrids_geo.f \
	$(GEOLIB)/allowflag_geo.f  \
	$(GEOLIB)/dumpgauge_geo.f \
	$(GEOLIB)/setgrd_geo.f \
	$(GEOLIB)/gfixup_geo.f \
	$(GEOLIB)/ginit_geo.f \
	$(GEOLIB)/getmaxspeed_geo.f \
	$(CLAW_LIB)/prefilp.f \
	$(CLAW_LIB)/quick_sort1.f \
	$(CLAW_LIB)/advanc.f \
	$(CLAW_LIB)/estdt.f \
	$(CLAW_LIB)/bound.f \
	$(CLAW_LIB)/check4nans.f \
	$(CLAW_LIB)/spest.f  \
	$(CLAW_LIB)/igetsp.f \
	$(CLAW_LIB)/reclam.f \
	$(CLAW_LIB)/birect.f \
	$(CLAW_LIB)/cleanup.f \
	$(CLAW_LIB)/colate.f \
	$(CLAW_LIB)/errest.f \
	$(CLAW_LIB)/bufnst.f \
	$(CLAW_LIB)/errf1.f \
	$(CLAW_LIB)/fixcapaq.f \
	$(CLAW_LIB)/flglvl.f \
	$(CLAW_LIB)/fluxad.f \
	$(CLAW_LIB)/fluxsv.f \
	$(CLAW_LIB)/grdfit.f \
	$(CLAW_LIB)/intfil.f \
	$(CLAW_LIB)/moment.f \
	$(CLAW_LIB)/nestck.f \
	$(CLAW_LIB)/prepc.f \
	$(CLAW_LIB)/prepf.f \
	$(CLAW_LIB)/projec.f \
	$(CLAW_LIB)/signs.f \
	$(CLAW_LIB)/findcut.f \
	$(CLAW_LIB)/smartbis.f \
	$(CLAW_LIB)/putnod.f \
	$(CLAW_LIB)/putsp.f \
	$(CLAW_LIB)/regrid.f \
	$(CLAW_LIB)/setuse.f \
	$(CLAW_LIB)/stst1.f \
	$(CLAW_LIB)/trimbd.f \
	$(CLAW_LIB)/nodget.f \
	$(CLAW_LIB)/basic.f \
	$(CLAW_LIB)/outval.f \
	$(CLAW_LIB)/copysol.f \
	$(CLAW_LIB)/outvar.f \
	$(CLAW_LIB)/outmsh.f \
	$(CLAW_LIB)/outtre.f \
	$(CLAW_LIB)/domain.f  \
	$(CLAW_LIB)/cellave.f \
	$(CLAW_LIB)/fdisc.f \
	$(CLAW_LIB)/fss.f \
	$(CLAW_LIB)/zeroin.f \
	$(CLAW_LIB)/setflags.f \
	$(CLAW_LIB)/shiftset.f \
	$(CLAW_LIB)/conck.f \
	$(CLAW_LIB)/domshrink.f \
	$(CLAW_LIB)/domprep.f \
	$(CLAW_LIB)/domup.f \
	$(CLAW_LIB)/domcopy.f \
	$(CLAW_LIB)/coarsen.f \
	$(CLAW_LIB)/auxcoarsen.f \
	$(CLAW_LIB)/intcopy.f \
	$(CLAW_LIB)/preintcopy.f \
	$(CLAW_LIB)/icall.f \
	$(CLAW_LIB)/preicall.f \
	$(CLAW_LIB)/inlinelimiter.f \
	$(CLAW_LIB)/cstore.f \
	$(CLAW_LIB)/saveqc.f \
	$(CLAW_LIB)/check.f \
	$(CLAW_LIB)/restrt.f \
	$(CLAW_LIB)/opendatafile.f \
	$(CLAW_LIB)/init_alloc.f90 \
	$(CLAW_LIB)/restrt_alloc.f90 \
	$(CLAW_LIB)/resize_alloc.f90 \
	
	# Library files being replaced
	# $(GEOLIB)/qinit_geo.f \
	# $(GEOLIB)/setprob_geo.f90 \
	# $(GEOLIB)/rpn2ez_geo.f \
	# $(GEOLIB)/rpt2_geo.f \
	# $(GEOLIB)/valout_geo.f \
	# $(GEOLIB)/flag2refine_geo.f  \
	# $(GEOLIB)/src2_geo.f \
	# $(GEOLIB)/src1d_geo.f \
	# $(GEOLIB)/setaux_geo.f \
	# $(GEOLIB)/b4step2_geo.f \
	
#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAW)/util/Makefile.common
	
### DO NOT remove this line - make depends on it ###
    
