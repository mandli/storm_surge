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
import sys
import os
import time

import numpy as np

from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

import setrun
import topo_data

# Parameters
if os.environ.has_key('DATA_PATH'):
    base_path = os.path.join(os.environ['DATA_PATH'],"storm_surge","hurricane")
else:
    base_path = os.getcwd()
parallel = True
process_queue = []
runclaw_cmd = "python $CLAW/python/pyclaw/runclaw.py"
plotclaw_cmd = "python $CLAW/python/pyclaw/plotters/plotclaw.py"

# Tests
tests = [{"name":"perpendicular","velocity":5.0, "angle": 0.00 * np.pi, "eye":(0.0,0.0)},
         {"name":"45angle","velocity":5.0, "angle": 0.25 * np.pi, "eye":(200e3,0.0)},
         {"name":"parallel","velocity":5.0, "angle": 0.50 * np.pi, "eye":(400e3,0.0)},]
         
if len(sys.argv) == 2: 
    tests = [tests[int(sys.argv[1])]]
print tests

# Setup and run the tests
for test in tests:    
    # Default data values
    # topo_data.write_topo_file('./topo.data',bathy_type='gulf_shelf',plot=False,force=False)
    rundata = setrun.setrun()
    hurricane_data = setrun.set_hurricane_data()
    multilayer_data = setrun.set_multilayer_data()

    # Set data particular to this test
    hurricane_data.hurricane_velocity = (test["velocity"] 
        * np.cos(test["angle"]),test["velocity"] * np.sin(test["angle"]))
    hurricane_data.R_eye_init = test["eye"]

    # Grid parameters
    factor = 3*2
    rundata.clawdata.mx = 70 * factor
    rundata.clawdata.my = 60 * factor
    rundata.clawdata.dt_initial = 0.6755e1
    rundata.clawdata.mxnest = -1

    # Write out data
    rundata.write()
    hurricane_data.write()
    multilayer_data.write()

    # Create output directory
    prefix = "ml_%s" % test['name']
    # tm = time.localtime(os.path.getmtime(outdir))
    tm = time.localtime()
    year = str(tm[0]).zfill(4)
    month = str(tm[1]).zfill(2)
    day = str(tm[2]).zfill(2)
    hour = str(tm[3]).zfill(2)
    minute = str(tm[4]).zfill(2)
    second = str(tm[5]).zfill(2)
    date = '_%s%s%s-%s%s%s' % (year,month,day,hour,minute,second)
    output_dirname = ''.join((prefix,date,"_output"))
    plots_dirname = ''.join((prefix,date,"_plots"))
    log_name = ''.join((prefix,date,"_log"))

    output_path = os.path.join(base_path,output_dirname)
    plots_path = os.path.join(base_path,plots_dirname)
    log_path = os.path.join(base_path,log_name)

    log_file = open(log_path,'w')

    # Run the simulation
    run_cmd = "%s xclaw %s" % (runclaw_cmd,output_path)
    plot_cmd = "%s %s %s" % (plotclaw_cmd,output_path,plots_path)
    tar_cmd = "tar -cvzf %s.tgz %s" % (plots_path,plots_path)
    cmd = ";".join((run_cmd,plot_cmd))
    print cmd

    # Run command
    if parallel:
        process_queue.append(subprocess.Popen(cmd,shell=True,
                                stdout=log_file,stderr=log_file))
    else:
        subprocess.Popen(cmd,shell=True,stdout=log_file,
		stderr=log_file).wait()
    
# if parallel:
#     for process in process_queue:
#         try:
#             process.wait()
#         except(KeyboardInterrupt,SystemExit):
#             print "Interrupt called and caught."
        
    

      
