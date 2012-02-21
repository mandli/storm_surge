#
# Local makefile for test_beach
#

# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = geoclaw
CLAW_setrun_file = setrun.py
CLAW_setplot_file = setplot.py
CLAW_OUTDIR = _output
CLAW_PLOTDIR = _plots

# Clawpack libraries to be used:
CLAW_LIB = $(AMRCLAW)/src/2d
GEOLIB = $(GEOCLAW)/src/2d

# Flagging
FC ?= gfortran
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
MODULES = \
    $(AMRCLAW)/src/2d/amr_module.f \
	$(ML_PATH)/src/2d/multilayer_module.f90 \
    $(GEOCLAW)/src/2d/geoclaw_module.f90 \
    $(GEOCLAW)/src/2d/topo_module.f90 \
    $(GEOCLAW)/src/2d/dtopo_module.f90 \
    $(GEOCLAW)/src/2d/regions_module.f90 \
    $(GEOCLAW)/src/2d/qinit_module.f90 \
    $(GEOCLAW)/src/2d/fixedgrids_module.f90 \
    $(GEOCLAW)/src/2d/gauges_module.f90

SOURCES = \
	$(ML_PATH)/src/2d/hurricane.f90 \
	$(ML_PATH)/src/2d/setaux.f \
	$(ML_PATH)/src/2d/b4step2.f \
	$(ML_PATH)/src/2d/setprob.f90 \
	$(ML_PATH)/src/2d/src2.f90 \
	$(ML_PATH)/src/2d/src1d.f90 \
	$(ML_PATH)/src/2d/qinit.f90 \
	$(ML_PATH)/src/2d/valout.f \
	$(ML_PATH)/src/2d/flag2refine.f \
	$(ML_PATH)/src/2d/rpn2_inundation_multilayer.f90 \
	$(ML_PATH)/src/2d/rpt2_dry_multilayer.f90 \
	$(ML_PATH)/src/2d/multilayer_eigen.f90 \
	$(ML_PATH)/src/2d/tick.f \
	$(ML_PATH)/src/2d/stepgrid.f \
	$(ML_PATH)/src/2d/flux2fw.f \
    $(ML_PATH)/src/2d/dumpgauge_ml.f \
    $(GEOCLAW)/src/2d/qad.f \
    $(GEOCLAW)/src/2d/readqinit.f \
    $(GEOCLAW)/src/2d/movetopo.f \
    $(GEOCLAW)/src/2d/cellgridintegrate.f \
    $(GEOCLAW)/src/2d/filval.f \
    $(GEOCLAW)/src/2d/topointegral.f \
    $(GEOCLAW)/src/2d/bilinearintegral.f \
    $(GEOCLAW)/src/2d/fgridinterp.f \
    $(GEOCLAW)/src/2d/fgridout.f \
    $(GEOCLAW)/src/2d/riemannsolvers.f \
    $(GEOCLAW)/src/2d/step2.f \
    $(GEOCLAW)/src/2d/filpatch.f \
    $(GEOCLAW)/src/2d/bc2amr.f \
    $(GEOCLAW)/src/2d/update.f \
    $(GEOCLAW)/src/2d/amr2ez.f \
    $(GEOCLAW)/src/2d/setqinit.f \
    $(GEOCLAW)/src/2d/setregions.f \
    $(GEOCLAW)/src/2d/setgauges.f \
    $(GEOCLAW)/src/2d/setfixedgrids.f \
    $(GEOCLAW)/src/2d/allowflag.f  \
    $(GEOCLAW)/src/2d/upbnd.f  \
    $(GEOCLAW)/src/2d/setgrd.f \
    $(GEOCLAW)/src/2d/gfixup.f \
    $(GEOCLAW)/src/2d/ginit.f \
    $(GEOCLAW)/src/2d/getmaxspeed.f \
    $(AMRCLAW)/src/2d/prefilp.f \
    $(AMRCLAW)/src/2d/quick_sort1.f \
    $(AMRCLAW)/src/2d/advanc.f \
    $(AMRCLAW)/src/2d/estdt.f \
    $(AMRCLAW)/src/2d/bound.f \
    $(AMRCLAW)/src/2d/check4nans.f \
    $(AMRCLAW)/src/2d/spest.f  \
    $(AMRCLAW)/src/2d/igetsp.f \
    $(AMRCLAW)/src/2d/reclam.f \
    $(AMRCLAW)/src/2d/birect.f \
    $(AMRCLAW)/src/2d/cleanup.f \
    $(AMRCLAW)/src/2d/colate.f \
    $(AMRCLAW)/src/2d/errest.f \
    $(AMRCLAW)/src/2d/bufnst.f \
    $(AMRCLAW)/src/2d/errf1.f \
    $(AMRCLAW)/src/2d/fixcapaq.f \
    $(AMRCLAW)/src/2d/flglvl.f \
    $(AMRCLAW)/src/2d/fluxad.f \
    $(AMRCLAW)/src/2d/fluxsv.f \
    $(AMRCLAW)/src/2d/grdfit.f \
    $(AMRCLAW)/src/2d/intfil.f \
    $(AMRCLAW)/src/2d/moment.f \
    $(AMRCLAW)/src/2d/nestck.f \
    $(AMRCLAW)/src/2d/prepc.f \
    $(AMRCLAW)/src/2d/prepf.f \
    $(AMRCLAW)/src/2d/projec.f \
    $(AMRCLAW)/src/2d/signs.f \
    $(AMRCLAW)/src/2d/findcut.f \
    $(AMRCLAW)/src/2d/smartbis.f \
    $(AMRCLAW)/src/2d/putnod.f \
    $(AMRCLAW)/src/2d/putsp.f \
    $(AMRCLAW)/src/2d/regrid.f \
    $(AMRCLAW)/src/2d/setuse.f \
    $(AMRCLAW)/src/2d/stst1.f \
    $(AMRCLAW)/src/2d/trimbd.f \
    $(AMRCLAW)/src/2d/nodget.f \
    $(AMRCLAW)/src/2d/basic.f \
    $(AMRCLAW)/src/2d/outval.f \
    $(AMRCLAW)/src/2d/copysol.f \
    $(AMRCLAW)/src/2d/outvar.f \
    $(AMRCLAW)/src/2d/outmsh.f \
    $(AMRCLAW)/src/2d/outtre.f \
    $(AMRCLAW)/src/2d/domain.f  \
    $(AMRCLAW)/src/2d/cellave.f \
    $(AMRCLAW)/src/2d/fdisc.f \
    $(AMRCLAW)/src/2d/fss.f \
    $(AMRCLAW)/src/2d/zeroin.f \
    $(AMRCLAW)/src/2d/setflags.f \
    $(AMRCLAW)/src/2d/shiftset.f \
    $(AMRCLAW)/src/2d/conck.f \
    $(AMRCLAW)/src/2d/domshrink.f \
    $(AMRCLAW)/src/2d/domprep.f \
    $(AMRCLAW)/src/2d/domup.f \
    $(AMRCLAW)/src/2d/domcopy.f \
    $(AMRCLAW)/src/2d/coarsen.f \
    $(AMRCLAW)/src/2d/auxcoarsen.f \
    $(AMRCLAW)/src/2d/intcopy.f \
    $(AMRCLAW)/src/2d/preintcopy.f \
    $(AMRCLAW)/src/2d/icall.f \
    $(AMRCLAW)/src/2d/preicall.f \
    $(AMRCLAW)/src/2d/inlinelimiter.f \
    $(AMRCLAW)/src/2d/cstore.f \
    $(AMRCLAW)/src/2d/saveqc.f \
    $(AMRCLAW)/src/2d/check.f \
    $(AMRCLAW)/src/2d/restrt.f \
    $(AMRCLAW)/src/2d/opendatafile.f \
    $(AMRCLAW)/src/2d/init_alloc.f90 \
    $(AMRCLAW)/src/2d/restrt_alloc.f90 \
    $(AMRCLAW)/src/2d/resize_alloc.f90
	
#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWUTIL)/src/Makefile.common
	
### DO NOT remove this line - make depends on it ###
    
