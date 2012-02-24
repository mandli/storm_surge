
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAWUTIL)/src/Makefile.common

# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Flagging
PPFLAGS ?= 
UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	LAPACK_LIB = -framework veclib
    # LFLAGS = -llapack -L/opt/local/lib
endif
ifeq ($(UNAME), Linux)
	LAPACK_LIB = -llapack
endif
ifeq ($(findstring gfortran,$(FC)),gfortran)
    PPFLAGS += -cpp
endif
ifeq ($(findstring ifort,$(FC)),ifort)
    PPFLAGS += -fpp
endif
MAKELEVEL ?= 0
ifeq ($(MAKELEVEL),0)
    FFLAGS += $(PPFLAGS) $(LAPACK_LIB)
    LFLAGS += $(LAPACK_LIB)
endif

# ---------------------------------
# List of sources for this program:
# ---------------------------------
MODULES = \
	$(ML_PATH)/src/1d/parameters_module.f90

SOURCES = \
	$(ML_PATH)/src/1d/driver.f \
	$(ML_PATH)/src/1d/qinit.f90 \
	$(ML_PATH)/src/1d/rp1_inundation.f90 \
	$(ML_PATH)/src/1d/out1.f \
	$(ML_PATH)/src/1d/b4step1.f90 \
	$(ML_PATH)/src/1d/bc1.f \
	$(ML_PATH)/src/1d/setaux.f90 \
	$(ML_PATH)/src/1d/src1.f90 \
	$(ML_PATH)/src/1d/setprob.f90 \
	$(CLAW_4)/clawpack/1d/lib/claw1ez.f \
	$(CLAW_4)/clawpack/1d/lib/claw1.f \
	$(CLAW_4)/clawpack/1d/lib/step1fw.f \
	$(CLAW_4)/clawpack/1d/lib/copyq1.f \
	$(CLAW_4)/clawpack/1d/lib/limiter.f \
	$(CLAW_4)/clawpack/1d/lib/philim.f \
	$(CLAW_4)/clawpack/1d/lib/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
	
### DO NOT remove this line - make depends on it ###
    
