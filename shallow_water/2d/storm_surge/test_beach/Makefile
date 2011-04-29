#
# Local makefile for test_beach
#

# Path to storm surge specific source
HURRICANE_SRC = ../src

# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = geoclaw
LIB_DIM = 2d
SETRUN_FILE = setrun.py
SETPLOT_FILE = setplot.py
OUT_DIR = _output
PLOT_DIR = _plots
MOD_PATH = ./

# ---------------------------------
# List of sources for this program:
# ---------------------------------
SRC = \
	$(HURRICANE_SRC)/setprob.f90 \
	$(HURRICANE_SRC)/b4step2_geo.f \
	$(HURRICANE_SRC)/valout_geo.f \
	$(HURRICANE_SRC)/setaux_geo.f \
	$(HURRICANE_SRC)/flag2refine_geo.f \
	$(HURRICANE_SRC)/src2_geo.f90 \
	$(HURRICANE_SRC)/src1d_geo.f90

MOD_SRC = \
    $(HURRICANE_SRC)/hurricane.f90

MODULES = hurricane_module

# $(HURRICANE_SRC)/riemann_fwave_stormsurge.f90 \
# $(HURRICANE_SRC)/rpn2ez_geo.f
# $(HURRICANE_SRC)/qad_geo.f 
# $(HURRICANE_SRC)/filval_geo.f \
# $(HURRICANE_SRC)/filpatch_geo.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAW)/util/local_rules.mk
	
### DO NOT remove this line - make depends on it ###
    
