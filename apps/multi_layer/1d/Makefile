
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# Adjust these variables if desired:
# ----------------------------------
CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots
# LFLAGS = -framework veclib
LFLAGS = -llapack -L/opt/local/lib

# ---------------------------------
# List of sources for this program:
# ---------------------------------
CLAW_SOURCES = \
	$(ML_SRC)/1d/setprob.f90 \
	$(ML_SRC)/1d/driver.f \
	$(ML_SRC)/1d/qinit.f90 \
	$(ML_SRC)/1d/rp1_inundation.f90 \
	$(ML_SRC)/1d/out1.f \
	$(ML_SRC)/1d/b4step1.f90 \
	$(ML_SRC)/1d/bc1.f \
	$(ML_SRC)/1d/setaux.f90 \
    $(ML_SRC)/1d/src1.f90
	
CLAW_LIB = $(CLAW)/clawpack/1d/lib
	
CLAW_LIBSOURCES = \
  $(CLAW_LIB)/claw1ez.f \
  $(CLAW_LIB)/claw1.f \
  $(CLAW_LIB)/step1fw.f \
  $(CLAW_LIB)/copyq1.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
	
### DO NOT remove this line - make depends on it ###
    
