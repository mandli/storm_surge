#!/usr/bin/env python
# encoding: utf-8
r"""
Runs tests for angles relative to a continental shelf

:Authors:
    Kyle Mandli (2011-04-20) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import subprocess

import numpy as np

from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

import setrun
import topo_data

def run_simulation(file_suffix):
    runclaw(xclawcmd='xclaw',outdir='_output%s' % file_suffix)
    plotclaw(outdir='_output%s' % file_suffix,plotdir='_plots%s' % file_suffix)


# Tests
tests = [{"velocity":5.0, "angle": 0.00 * np.pi, "eye":(0.0,0.0)},
         {"velocity":5.0, "angle": 0.25 * np.pi, "eye":(200e3,0.0)},
         {"velocity":5.0, "angle": 0.50 * np.pi, "eye":(400e3,0.0)},]

# Setup and run the tests
for (i,test) in enumerate(tests):    
    # Default data values
    # This is no longer necessary as it's in the setaux routine directly
    topo_data.write_topo_file('./topo.data',bathy_type='gulf_shelf',plot=False,force=False)
    rundata = setrun.setrun()
    hurricane_data = setrun.set_hurricane_data()
    multilayer_data = setrun.set_multilayer_data()
    
    # Set data particular to this test
    hurricane_data.hurricane_velocity = (test["velocity"] 
        * np.cos(test["angle"]),test["velocity"] * np.sin(test["angle"]))
    hurricane_data.R_eye_init = test["eye"]
    
    # Write out data
    rundata.write()
    hurricane_data.write()
    multilayer_data.write()
    
    # Run the simulation
    prefix = "_sl_%s" % i
    run_simulation(prefix)
    
    # Tar up the results
    cmd = "tar -cvzf ~/sl_%s_plots.tgz _plots_sl_%s" % (i,i)
    print cmd
    subprocess.Popen(cmd,shell=True).wait()
    
    
