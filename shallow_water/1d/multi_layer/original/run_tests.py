#!/usr/bin/env python
# encoding: utf-8
r"""
Run convergence tests for some of the multilayer test cases
"""

import shutil

import numpy as np

# Pyclaw functionality
from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

# Local imports
from setrun import setrun,MultilayerData
from setplot import setplot

def run_simulation(file_suffix):
    runclaw(xclawcmd='xclaw',outdir='_output%s' % file_suffix)
    plotclaw(outdir='_output%s' % file_suffix,plotdir='_plots%s' % file_suffix)

# Initialize basic rundata objects
run_data = setrun()

# ============================================================================
#  Hurricane run
# ============================================================================
# RAMP_UP_TIME = 12*60**2
# run_data.clawdata.xlower = -200e3 
# run_data.clawdata.xupper = 500e3
# run_data.clawdata.mx = 700
# num_hours = 40
# step = 0.25
# run_data.clawdata.t0 = -RAMP_UP_TIME
# run_data.clawdata.nout = int(num_hours / step) + int(np.ceil(RAMP_UP_TIME / (step*60**2)))
# run_data.clawdata.tfinal = num_hours * 60.0**2
# 
# prob_data = MultilayerData(run_data)
# prob_data.init_type = 0
# prob_data.depth = 100.0
# prob_data.wind_type = 2
# prob_data.total_depth = 100.0
# prob_data.A = 23.0
# prob_data.B = 1.5
# prob_data.Pn = 1005.0
# prob_data.Pc = 950.0
# prob_data.R_eye_init = 0.0
# prob_data.hurricane_velocity = 5.0
# prob_data.ramp_up_time = RAMP_UP_TIME

# ============================================================================
#  Oscillatory wind, blow up
# ============================================================================
# print "Removing previous output and plot directories"
# for var in ['A','N','omega']:
#     for dir in ['_output','_plots']:
#         [shutil.rmtree("%s_%s_%s" % (dir,var,x),True) for x in xrange(1,6)]

# prob_data = MultilayerData(run_data)
# prob_data.init_type = 0
# prob_data.wind_type = 3
# prob_data.A = 5.0
# prob_data.N = 2.0
# prob_data.omega = 2.0

# for A in xrange(1,6):
#     prob_data.A = A
#     # Write out config files
#     run_data.write()
#     prob_data.write(file='./problem.data',datasource='run_tests.py')
#     
#     # Run simulation
#     run_simulation("_A_%s" % str(int(A)))
# 
# prob_data.A = 2.0
# for N in xrange(1,6):
#     prob_data.N = N
#     # Write out config files
#     run_data.write()
#     prob_data.write(file='./problem.data',datasource='run_tests.py')
#     
#     # Run simulation
#     run_simulation("_N_%s" % str(int(N)))
# 
# prob_data.N = 2.0
# for omega in xrange(1,6):
#     prob_data.omega = omega
#     # Write out config files
#     run_data.write()
#     prob_data.write(file='./problem.data',datasource='run_tests.py')
#     
#     # Run simulation
#     run_simulation("_omega_%s" % str(int(omega)))
    
# ============================================================================
#  Lake run - not that interesting, needs mass exchange apparently
# ============================================================================
# run_data.clawdata.xlower = 0.0
# run_data.clawdata.max_steps = 15000
# run_data.clawdata.xupper = 12e3
# run_data.clawdata.mx = 120
# num_hours = 10.0
# step = 0.25
# run_data.clawdata.nout = int(num_hours / step) + int(np.ceil(0.0 / (step*60**2)))
# run_data.clawdata.tfinal = num_hours * 60.0**2
# 
# prob_data = MultilayerData(run_data)
# prob_data.init_type = 0
# prob_data.bathy_type = 1
# prob_data.depth = 55.0
# prob_data.wind_type = 1
# prob_data.total_depth = 55.0
# prob_data.A = 10.0

# ============================================================================
#  Dry state problem
# ============================================================================
run_data.clawdata.xlower = -1.0
run_data.clawdata.xupper = 1.0
run_data.clawdata.mthbc_xlower = 1
run_data.clawdata.mthbc_xupper = 1

prob_data = MultilayerData(run_data)
prob_data.init_type = 4
prob_data.wave_family = 1
prob_data.bathy_type = 2
prob_data.wind_type = 0
prob_data.depth = -1.0

prob_data.init_location = 0.0
# ============================================================================
#  Run simulation
# ============================================================================
# Write out config files
run_data.write()
prob_data.write(file='./problem.data',datasource='run_tests.py')

# Run simulation
runclaw(xclawcmd='xclaw',outdir='_output')
# Plot output
plotclaw(outdir='_output',plotdir='_plots')




